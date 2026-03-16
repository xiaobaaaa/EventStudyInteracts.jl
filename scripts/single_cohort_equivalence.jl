using DataFrames
using EventStudyInteracts
using FixedEffectModels
using LinearAlgebra
using Random

const SINGLE_COHORT_REL_VARLIST = [:g_3, :g_2, :g0, :g1, :g2, :g3]
const SINGLE_COHORT_POST_VARLIST = [:g0, :g1, :g2, :g3]

function make_single_cohort_panel(; seed::Int = 20260316, n_ids::Int = 96, n_periods::Int = 10, treat_time::Int = 5)
    Random.seed!(seed)

    ids = repeat(1:n_ids, inner = n_periods)
    periods = repeat(1:n_periods, outer = n_ids)
    treated_ids = collect(1:(n_ids ÷ 2))
    treated = in.(ids, Ref(treated_ids))
    first_treat = ifelse.(treated, treat_time, missing)
    never_treat = .!treated
    rel_time = periods .- first_treat

    industry_by_id = mod1.(shuffle(1:n_ids), 6)
    industry = repeat(industry_by_id, inner = n_periods)

    unit_fe = randn(n_ids)
    time_fe = randn(n_periods)
    unit_slope = randn(n_ids)
    treated_time_shock = randn(n_periods)
    industry_time_shock = randn(maximum(industry_by_id), n_periods)
    id_noise = randn(n_ids, n_periods)
    iid_noise = randn(length(ids))

    dynamic_effect = Dict(-3 => 0.0, -2 => 0.0, 0 => 1.0, 1 => 1.4, 2 => 1.8, 3 => 2.1)
    effect = [get(dynamic_effect, Int(coalesce(rt, -999)), 0.0) for rt in rel_time]

    centered_t = periods .- ((n_periods + 1) / 2)
    residual = similar(iid_noise)
    for i in eachindex(ids)
        residual[i] = 0.35 * unit_slope[ids[i]] * centered_t[i] / n_periods +
            0.55 * treated[i] * treated_time_shock[periods[i]] +
            0.40 * industry_time_shock[industry[i], periods[i]] +
            0.25 * id_noise[ids[i], periods[i]] +
            0.20 * iid_noise[i]
    end

    y = similar(iid_noise)
    for i in eachindex(ids)
        y[i] = unit_fe[ids[i]] + time_fe[periods[i]] + effect[i] + residual[i]
    end

    df = DataFrame(
        id = ids,
        t = periods,
        industry = industry,
        first_treat = first_treat,
        never_treat = never_treat,
        rel_time = rel_time,
        y = y,
    )

    for k in 3:-1:2
        df[!, Symbol("g_$k")] = Int.(coalesce.(df.rel_time .== -k, false))
    end
    for k in 0:3
        df[!, Symbol("g$k")] = Int.(coalesce.(df.rel_time .== k, false))
    end

    return df
end

single_cohort_vcov_estimators() = [
    ("simple", Vcov.simple()),
    ("robust", Vcov.robust()),
    ("cluster(id)", Vcov.cluster(:id)),
    ("cluster(id, t)", Vcov.cluster(:id, :t)),
]

function fit_direct_single_cohort(df::DataFrame, vcov)
    formula = term(:y) ~ sum(term.(SINGLE_COHORT_REL_VARLIST)) + sum(fe.(term.([:id, :t])))
    return reg(df, formula, vcov; progress_bar = false)
end

function fit_iw_single_cohort(df::DataFrame, vcov)
    formula = term(:y) ~ sum(fe.(term.([:id, :t])))
    return eventreg(df, formula, SINGLE_COHORT_REL_VARLIST, :never_treat, :first_treat, vcov; progress_bar = false)
end

function _coef_lookup(model)
    return Dict(String(name) => value for (name, value) in zip(coefnames(model), coef(model)))
end

function _stderror_lookup(model)
    ses = sqrt.(diag(Matrix(vcov(model))))
    return Dict(String(name) => value for (name, value) in zip(coefnames(model), ses))
end

function average_effect(model, vars::Vector{Symbol} = collect(SINGLE_COHORT_POST_VARLIST))
    cn = String.(coefnames(model))
    idx = [findfirst(==(string(var)), cn) for var in vars]
    any(isnothing, idx) && error("Missing coefficients in average effect calculation.")
    weights = zeros(length(cn))
    weights[collect(idx)] .= 1 / length(vars)
    estimate = dot(weights, coef(model))
    variance = dot(weights, Matrix(vcov(model)) * weights)
    return (estimate = estimate, se = sqrt(variance))
end

function compare_single_cohort_models(direct, iw)
    direct_coef = _coef_lookup(direct)
    iw_coef = _coef_lookup(iw)
    direct_se = _stderror_lookup(direct)
    iw_se = _stderror_lookup(iw)

    period_rows = NamedTuple[]
    max_coef_diff = 0.0
    max_se_diff = 0.0
    for name in string.(SINGLE_COHORT_REL_VARLIST)
        coef_diff = abs(direct_coef[name] - iw_coef[name])
        se_diff = abs(direct_se[name] - iw_se[name])
        max_coef_diff = max(max_coef_diff, coef_diff)
        max_se_diff = max(max_se_diff, se_diff)
        push!(period_rows, (
            name = name,
            direct_estimate = direct_coef[name],
            iw_estimate = iw_coef[name],
            direct_se = direct_se[name],
            iw_se = iw_se[name],
            estimate_diff = coef_diff,
            se_diff = se_diff,
        ))
    end

    ate_direct = average_effect(direct)
    ate_iw = average_effect(iw)
    vcov_diff = maximum(abs.(Matrix(vcov(direct)) .- Matrix(vcov(iw))))

    return (
        period_rows = period_rows,
        max_abs_coef_diff = max_coef_diff,
        max_abs_se_diff = max_se_diff,
        max_abs_vcov_diff = vcov_diff,
        ate_direct = ate_direct.estimate,
        ate_iw = ate_iw.estimate,
        ate_estimate_diff = abs(ate_direct.estimate - ate_iw.estimate),
        ate_direct_se = ate_direct.se,
        ate_iw_se = ate_iw.se,
        ate_se_diff = abs(ate_direct.se - ate_iw.se),
    )
end

function run_single_cohort_equivalence()
    df = make_single_cohort_panel()
    results = NamedTuple[]
    for (label, estimator) in single_cohort_vcov_estimators()
        direct = fit_direct_single_cohort(df, estimator)
        iw = fit_iw_single_cohort(df, estimator)
        comparison = compare_single_cohort_models(direct, iw)
        push!(results, (label = label, direct = direct, iw = iw, comparison = comparison))
    end
    return (data = df, results = results)
end

function _fmt(x)
    return string(round(x; sigdigits = 8))
end

function render_single_cohort_markdown(report = run_single_cohort_equivalence())
    io = IOBuffer()
    println(io, "## Single-Cohort Equivalence Check")
    println(io)
    println(io, "When all treated units start treatment in the same period, the IW estimator collapses to the ordinary event-study regression because there is only one treated cohort and all cohort-share weights are equal to one.")
    println(io)
    println(io, "The script [`scripts/single_cohort_equivalence.jl`](scripts/single_cohort_equivalence.jl) generates a deterministic panel and compares `eventreg(...)` against a direct `FixedEffectModels.reg(...)` event-study regression under four covariance estimators.")
    println(io)
    println(io, "| Vcov | max abs coef diff | max abs se diff | max abs vcov diff | abs ATE diff | abs ATE se diff |")
    println(io, "| --- | ---: | ---: | ---: | ---: | ---: |")
    for row in report.results
        cmp = row.comparison
        println(io, "| $(row.label) | $(_fmt(cmp.max_abs_coef_diff)) | $(_fmt(cmp.max_abs_se_diff)) | $(_fmt(cmp.max_abs_vcov_diff)) | $(_fmt(cmp.ate_estimate_diff)) | $(_fmt(cmp.ate_se_diff)) |")
    end
    println(io)
    detail = only(filter(x -> x.label == "cluster(id, t)", report.results))
    println(io, "For the two-way clustered case (`cluster(id, t)`), the period-by-period estimates and the average post-treatment effect over `g0:g3` are:")
    println(io)
    println(io, "| Term | Direct FEM estimate | eventreg estimate | Direct FEM std. error | eventreg std. error |")
    println(io, "| --- | ---: | ---: | ---: | ---: |")
    for period in detail.comparison.period_rows
        println(io, "| $(period.name) | $(_fmt(period.direct_estimate)) | $(_fmt(period.iw_estimate)) | $(_fmt(period.direct_se)) | $(_fmt(period.iw_se)) |")
    end
    println(io, "| ATE(g0:g3) | $(_fmt(detail.comparison.ate_direct)) | $(_fmt(detail.comparison.ate_iw)) | $(_fmt(detail.comparison.ate_direct_se)) | $(_fmt(detail.comparison.ate_iw_se)) |")
    return String(take!(io))
end

function main(args = ARGS)
    output = nothing
    for arg in args
        startswith(arg, "--output=") || continue
        output = split(arg, "="; limit = 2)[2]
    end
    markdown = render_single_cohort_markdown()
    println(markdown)
    if output !== nothing
        write(output, markdown)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

