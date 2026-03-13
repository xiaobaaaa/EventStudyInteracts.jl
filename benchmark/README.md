# Benchmark Workflow

This folder contains the reproducible speed-comparison workflow used in the repository README.

## What It Benchmarks

The workflow compares the same synthetic Sun-Abraham style event-study panel across three implementations:

- `Julia EventStudyInteracts.jl`
- `R fixest` via `feols(y ~ sunab(first_treat_fixest, t, ref.p = c(-1, -17)) | ...)`
- `Stata eventstudyinteract`

Two fixed-effect specifications are timed:

- `id + t`
- `id + id1 + id2 + t`

## Files

- `generate_benchmark_data.jl`: creates the synthetic panel once and writes `benchmark/artifacts/benchmark_panel.csv`
- `run_julia_benchmark.jl`: times `eventreg(...)` after a warmup run
- `run_fixest_benchmark.R`: times `fixest::feols(... sunab(...))` after a warmup run
- `run_eventstudyinteract_benchmark.do`: times Stata `eventstudyinteract` after a warmup run
- `render_speed_comparison.py`: merges raw timing CSVs and renders the README SVG
- `run_all.ps1`: orchestrates the full local benchmark

## Usage

From the repository root on Windows PowerShell:

```powershell
./benchmark/run_all.ps1
```

By default the runner uses all logical threads reported by the host for both Julia and `fixest`.
You can still override them manually if you want to profile a smaller setting:

```powershell
./benchmark/run_all.ps1 -Units 100000 -Periods 20 -JuliaThreads 24 -FixestThreads 24
```

The workflow writes:

- `benchmark/results_latest.csv`
- `benchmark/results_latest.md`
- `benchmark/speed_comparison.svg`

Raw intermediate files go under `benchmark/artifacts/` and are ignored by git.

## Notes

- Julia timings exclude the first-run JIT compilation cost by doing one warmup estimation before the timed run.
- The benchmark now uses all host logical threads by default instead of Julia's single-thread default on Windows terminals.
- The R side uses `fixest::sunab(...)` with `ref.p = c(-1, -17)` so that the omitted periods match the Julia/Stata event window exactly.
- If the repository lives under a non-ASCII Windows path, Stata should be run through the ASCII junction path `F:\codex_eventstudyinteracts`.
