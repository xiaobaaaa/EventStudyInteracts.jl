param(
    [switch]$Online,
    [string]$Julia = "$HOME\.julia\juliaup\julia-1.10.11+0.x64.w64.mingw32\bin\julia.exe"
)

$repoRoot = Split-Path -Parent $PSScriptRoot
$depotRoot = Join-Path $repoRoot '.julia_depot'
$registryRoot = Join-Path $depotRoot 'registries'
$userDepot = Join-Path $HOME '.julia'

function Ensure-Junction {
    param(
        [string]$Path,
        [string]$Target
    )

    if (Test-Path $Path) {
        return
    }

    New-Item -ItemType Junction -Path $Path -Target $Target | Out-Null
}

function Ensure-GeneralRegistry {
    $generalRoot = Join-Path $registryRoot 'General'
    $registryFile = Join-Path $generalRoot 'Registry.toml'
    if (Test-Path $registryFile) {
        return
    }

    New-Item -ItemType Directory -Force -Path $registryRoot | Out-Null
    $tempRoot = Join-Path $registryRoot '__extract__'
    if (Test-Path $tempRoot) {
        Remove-Item -Recurse -Force $tempRoot
    }
    New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
    tar -xzf (Join-Path $userDepot 'registries\General.tar.gz') -C $tempRoot
    Move-Item -Path $tempRoot -Destination $generalRoot
}

New-Item -ItemType Directory -Force -Path $depotRoot | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $depotRoot 'logs') | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $depotRoot 'compiled') | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $depotRoot 'environments\v1.10') | Out-Null
Ensure-GeneralRegistry
Ensure-Junction -Path (Join-Path $depotRoot 'packages') -Target (Join-Path $userDepot 'packages')
Ensure-Junction -Path (Join-Path $depotRoot 'artifacts') -Target (Join-Path $userDepot 'artifacts')

$env:JULIA_DEPOT_PATH = $depotRoot
$env:JULIA_PKG_OFFLINE = if ($Online) { 'false' } else { 'true' }

Push-Location $repoRoot
try {
    $manifest = Join-Path $repoRoot 'Manifest.toml'
    if (Test-Path $manifest) {
        Copy-Item -Force -Path $manifest -Destination (Join-Path $repoRoot 'Manifest.toml.bak')
        Remove-Item -Force $manifest
    }
    & $Julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
}
finally {
    Pop-Location
}
