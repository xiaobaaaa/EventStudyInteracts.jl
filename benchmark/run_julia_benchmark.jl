
using CSV
using DataFrames
using EventStudyInteracts
using FixedEffectModels
using StatsModels

function parse_cli_options(args)
    options = Dict{String, String}()
    for arg in args
        startswith(arg, "--") || continue
        key_value = split(arg[3:end], "="; limit = 2)
        length(key_value) == 2 || continue
        options[key_value[1]] = key_value[2]
    end
    return options
end

function run_spec(df::DataFrame, fes::Vector{Symbol}, rel_varlist::Vector{Symbol})
    formula = term(:Y) ~ sum(fe.(term.(fes)))
    eventreg(df, formula, rel_varlist, :never_treat, :first_treat, Vcov.simple(); progress_bar = false)
    GC.gc()
    elapsed = @elapsed eventreg(df, formula, rel_varlist, :never_treat, :first_treat, Vcov.simple(); progress_bar = false)
    return elapsed
end

options = parse_cli_options(ARGS)
input_path = get(options, "input", joinpath(@__DIR__, "artifacts", "benchmark_panel.csv"))
output_path = get(options, "output", joinpath(@__DIR__, "artifacts", "julia_results.csv"))
lead_max = parse(Int, get(options, "lead-max", "16"))
lag_max = parse(Int, get(options, "lag-max", "16"))

df = CSV.read(input_path, DataFrame)
rel_varlist = vcat(Symbol.("g_" .* string.(lead_max:-1:2)), Symbol.("g" .* string.(0:lag_max)))

specs = [
    ("id + t", [:id, :t]),
    ("id + id1 + id2 + t", [:id, :id1, :id2, :t]),
]

results = DataFrame(
    engine = String[],
    spec = String[],
    seconds = Float64[],
    nobs = Int[],
)

for (spec_name, fes) in specs
    elapsed = run_spec(df, fes, rel_varlist)
    push!(results, ("Julia EventStudyInteracts.jl", spec_name, elapsed, nrow(df)))
    println("Julia | $(spec_name) | $(round(elapsed, digits = 3)) seconds")
end

mkpath(dirname(output_path))
CSV.write(output_path, results)
