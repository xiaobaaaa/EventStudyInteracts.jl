param(
    [int]$Units = 100000,
    [int]$Periods = 20,
    [int]$Seed = 20260313,
    [int]$JuliaThreads = 0,
    [int]$FixestThreads = 0,
    [string]$JuliaExe = "$HOME\.julia\juliaup\julia-1.10.11+0.x64.w64.mingw32\bin\julia.exe",
    [string]$RscriptExe = "C:\Program Files\R\R-4.5.2\bin\Rscript.exe",
    [string]$StataExe = "C:\Program Files\Stata18\StataMP-64.exe",
    [string]$AsciiRepo = "F:\codex_eventstudyinteracts"
)

function Invoke-NativeCommand {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [string[]]$Arguments = @()
    )

    & $FilePath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $FilePath $($Arguments -join ' ')"
    }
}

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$BenchmarkRoot = Join-Path $RepoRoot 'benchmark'
$Artifacts = Join-Path $BenchmarkRoot 'artifacts'
New-Item -ItemType Directory -Path $Artifacts -Force | Out-Null

$HostThreads = [Environment]::ProcessorCount
$ResolvedJuliaThreads = if ($JuliaThreads -gt 0) { $JuliaThreads } else { $HostThreads }
$ResolvedFixestThreads = if ($FixestThreads -gt 0) { $FixestThreads } else { $HostThreads }
$env:JULIA_NUM_THREADS = [string]$ResolvedJuliaThreads

$DataPath = Join-Path $Artifacts 'benchmark_panel.csv'
$JuliaOut = Join-Path $Artifacts 'julia_results.csv'
$FixestOut = Join-Path $Artifacts 'fixest_results.csv'
$SvgOut = Join-Path $BenchmarkRoot 'speed_comparison.svg'
$CsvOut = Join-Path $BenchmarkRoot 'results_latest.csv'
$SummaryOut = Join-Path $BenchmarkRoot 'results_latest.md'

Invoke-NativeCommand -FilePath $JuliaExe -Arguments @("--project=$BenchmarkRoot", (Join-Path $BenchmarkRoot 'generate_benchmark_data.jl'), "--output=$DataPath", "--units=$Units", "--periods=$Periods", "--seed=$Seed")
Invoke-NativeCommand -FilePath $JuliaExe -Arguments @("--project=$BenchmarkRoot", (Join-Path $BenchmarkRoot 'run_julia_benchmark.jl'), "--input=$DataPath", "--output=$JuliaOut")
Invoke-NativeCommand -FilePath $RscriptExe -Arguments @((Join-Path $BenchmarkRoot 'run_fixest_benchmark.R'), "--input=$DataPath", "--output=$FixestOut", "--threads=$ResolvedFixestThreads")

$stataRepo = if ((Test-Path $AsciiRepo) -and ($RepoRoot -ne $AsciiRepo)) { $AsciiRepo } else { $RepoRoot }
$stataArtifacts = Join-Path $stataRepo 'benchmark\artifacts'
New-Item -ItemType Directory -Path $stataArtifacts -Force | Out-Null
$stataDo = Join-Path $stataArtifacts '.run_eventstudyinteract_benchmark.do'
$stataLog = Join-Path $stataArtifacts 'stata_benchmark.log'
$stataInput = Join-Path $stataRepo 'benchmark\artifacts\benchmark_panel.csv'
$stataOutput = Join-Path $stataRepo 'benchmark\artifacts\stata_results.csv'
@"
version 18.0
cap log close
log using `"$stataLog`", replace text
clear all
import delimited using `"$stataInput`", clear
scalar bench_nobs = _N
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.t)
timer clear 1
timer on 1
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.t)
timer off 1
timer list 1
scalar t1 = r(t1)
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.id1 i.id2 i.t)
timer clear 2
timer on 2
qui eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat) absorb(i.id i.id1 i.id2 i.t)
timer off 2
timer list 2
scalar t2 = r(t2)
clear
input str28 engine str24 spec double seconds long nobs
"Stata eventstudyinteract" "id + t" . .
"Stata eventstudyinteract" "id + id1 + id2 + t" . .
end
replace seconds = t1 in 1
replace seconds = t2 in 2
replace nobs = bench_nobs in 1/2
export delimited using `"$stataOutput`", replace
log close
exit, clear
"@ | Set-Content -Path $stataDo -Encoding ascii
try {
    $p = Start-Process -FilePath $StataExe -ArgumentList '/e', 'do', $stataDo -Wait -PassThru
    if ($p.ExitCode -ne 0) {
        throw "Stata benchmark failed with exit code $($p.ExitCode)."
    }
} finally {
    Remove-Item $stataDo -Force -ErrorAction SilentlyContinue
}

$StataOut = $stataOutput
Invoke-NativeCommand -FilePath 'python' -Arguments @((Join-Path $BenchmarkRoot 'render_speed_comparison.py'), "--inputs=$JuliaOut,$FixestOut,$StataOut", "--svg=$SvgOut", "--csv=$CsvOut", "--summary=$SummaryOut")

Write-Host "Benchmark complete."
Write-Host "Host logical threads: $HostThreads"
Write-Host "Julia threads: $ResolvedJuliaThreads"
Write-Host "Fixest threads: $ResolvedFixestThreads"
Write-Host "Results CSV: $CsvOut"
Write-Host "Figure: $SvgOut"
