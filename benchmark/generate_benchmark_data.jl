
using CSV
using DataFrames
using Random
using Statistics

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

options = parse_cli_options(ARGS)
output_path = get(options, "output", joinpath(@__DIR__, "artifacts", "benchmark_panel.csv"))
units = parse(Int, get(options, "units", "100000"))
periods = parse(Int, get(options, "periods", "20"))
seed = parse(Int, get(options, "seed", "20260313"))
lead_max = parse(Int, get(options, "lead-max", "16"))
lag_max = parse(Int, get(options, "lag-max", "16"))

rng = MersenneTwister(seed)
obs = units * periods

id = repeat(1:units, inner = periods)
t = repeat(1:periods, outer = units)
id1_unit = rand(rng, 1:max(2000, units ÷ 10), units)
id2_unit = rand(rng, 1:100, units)
unit_fe = randn(rng, units)

treat_window = collect(4:max(4, periods - 2))
first_treat_unit = Vector{Union{Missing, Int}}(undef, units)
for i in eachindex(first_treat_unit)
    if rand(rng) < 0.25
        first_treat_unit[i] = missing
    else
        first_treat_unit[i] = rand(rng, treat_window)
    end
end

cohort_scale_unit = [ismissing(ft) ? 0.0 : 0.8 + 0.04 * ft for ft in first_treat_unit]

id1 = repeat(id1_unit, inner = periods)
id2 = repeat(id2_unit, inner = periods)
first_treat = repeat(first_treat_unit, inner = periods)
first_treat_fixest = Int.(coalesce.(first_treat, periods + 100))
never_treat = Int.(ismissing.(first_treat))
unit_fe_obs = repeat(unit_fe, inner = periods)
cohort_scale = repeat(cohort_scale_unit, inner = periods)

id1_effect = randn(rng, maximum(id1_unit)) .* 0.35
id2_effect = randn(rng, maximum(id2_unit)) .* 0.25

rel_time = Vector{Union{Missing, Int}}(undef, obs)
treated = zeros(Int, obs)
dynamic = zeros(Int, obs)

for i in eachindex(first_treat)
    ft = first_treat[i]
    if ismissing(ft)
        rel_time[i] = missing
        continue
    end

    rt = t[i] - ft
    rel_time[i] = rt
    if rt >= 0
        treated[i] = 1
        dynamic[i] = rt
    end
end

noise = randn(rng, obs) .* 0.8
trend = 0.2 .* t
heterogeneous_effect = treated .* (0.15 .* dynamic .+ cohort_scale .* log1p.(dynamic))
Y = 1.5 .+ unit_fe_obs .+ trend .+ id1_effect[id1] .+ id2_effect[id2] .+ heterogeneous_effect .+ noise

df = DataFrame(
    id = id,
    t = t,
    id1 = id1,
    id2 = id2,
    Y = Y,
    never_treat = never_treat,
    first_treat = first_treat,
    first_treat_fixest = first_treat_fixest,
    rel_time = rel_time,
)

for k in lead_max:-1:2
    df[!, Symbol("g_$(k)")] = Int.(coalesce.(df.rel_time .== -k, false))
end

for k in 0:lag_max
    df[!, Symbol("g$(k)")] = Int.(coalesce.(df.rel_time .== k, false))
end

mkpath(dirname(output_path))
CSV.write(output_path, df)

println("Wrote benchmark data to $(output_path)")
println("Rows: $(nrow(df))")
println("Columns: $(ncol(df))")
println("Never treated share: $(round(mean(df.never_treat), digits = 4))")
