using DataFrames
using Dates
using EventStudyInteracts
using FixedEffectModels
using LinearAlgebra
using ReadStatTables
using RegressionTables
using Test
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE_DIR = joinpath(REPO_ROOT, "reference")

include(joinpath(REPO_ROOT, "scripts", "single_cohort_equivalence.jl"))

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

function reference_report_path()
    options = parse_cli_options(ARGS)
    return get(options, "reference-report", nothing)
end

function make_test_panel()
    ids = repeat(1:6, inner = 5)
    periods = repeat(1:5, outer = 6)
    first_treat = repeat(Union{Missing, Int}[3, 3, 4, 4, missing, missing], inner = 5)

    df = DataFrame(id = ids, t = periods, first_treat = first_treat)
    df.never_treat = ismissing.(df.first_treat)
    df.D = Int.(coalesce.(df.t .>= df.first_treat, false))
    df.rel_time = df.t .- df.first_treat
    df.Y = 1.0 .* df.id .+ 0.5 .* df.t .+ 2.0 .* df.D

    for k in 2:-1:2
        df[!, Symbol("g_$k")] = Int.(coalesce.(df.rel_time .== -k, false))
    end

    for k in 0:2
        df[!, Symbol("g$k")] = Int.(coalesce.(df.rel_time .== k, false))
    end

    return df
end

function make_clustered_share_inputs()
    cluster_a = repeat(1:3, inner = 4)
    cluster_b = repeat(1:4, outer = 3)

    x1 = Float64[1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1]
    x2 = Float64[0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1]
    e1 = 0.5 .* cluster_a .+ 0.25 .* cluster_b .+ x1
    e2 = -0.3 .* cluster_a .+ 0.4 .* cluster_b .- x2

    df = DataFrame(cluster_a = cluster_a, cluster_b = cluster_b)
    X = hcat(x1, x2)
    e = hcat(Float64.(e1), Float64.(e2))
    return df, X, e
end

function manual_robust_share_vcov(X, e)
    N = size(X, 1)
    Sxxi = Matrix(-FixedEffectModels.invsym!(Symmetric(copy((X' * X) / N))))
    K = kron(Matrix{Float64}(I, size(e, 2), size(e, 2)), Sxxi)
    return Symmetric(K * Matrix(EventStudyInteracts.avar(e, X, true)) * K / N)
end

function min_skipmissing(x)
    vals = collect(skipmissing(x))
    return isempty(vals) ? missing : minimum(vals)
end

function available_reference_files()
    isdir(REFERENCE_DIR) || return String[]
    files = readdir(REFERENCE_DIR; join = true)
    return sort(filter(path -> endswith(path, ".toml") && !endswith(path, ".example.toml"), files))
end

function metric_value(model, name::String)
    if name == "adjr2"
        return adjr2(model)
    elseif name == "coef_count"
        return length(coef(model))
    elseif name == "dof"
        return dof(model)
    elseif name == "f_stat"
        return model.F
    elseif name == "iterations"
        return something(model.iterations, -1)
    elseif name == "nobs"
        return nobs(model)
    elseif name == "p_value"
        return model.p
    elseif name == "r2"
        return r2(model)
    elseif name == "r2_within"
        return something(model.r2_within, NaN)
    end

    error("Unsupported metric '$name' in reference file.")
end

function build_nlswork_model(spec::Dict{String, Any})
    dataset_path = joinpath(REPO_ROOT, spec["dataset"])
    df = DataFrame(readstat(dataset_path))

    df.union_year = ifelse.(coalesce.(df.union, false) .== 1, df.year, missing)
    transform!(groupby(df, :idcode), :union_year => min_skipmissing => :first_union)
    select!(df, Not(:union_year))

    df.ry = df.year - df.first_union
    df.never_union = ismissing.(df.first_union)

    relmin = abs(minimum(skipmissing(df.ry)))
    relmax = abs(maximum(skipmissing(df.ry)))

    for k in relmin:-1:2
        df[!, Symbol("g_$k")] = Int.(coalesce.(df.ry .== -k, false))
    end

    for k in 0:relmax
        df[!, Symbol("g$k")] = Int.(coalesce.(df.ry .== k, false))
    end

    for name in get(spec, "force_zero", String[])
        column = Symbol(name)
        if hasproperty(df, column)
            df[!, column] .= 0
        else
            df[!, column] = zeros(Int, nrow(df))
        end
    end

    absorb = Symbol.(spec["absorb"])
    formula = term(Symbol(spec["response"])) ~ term(:south) + sum(fe.(term.(absorb)))

    return eventreg(
        df,
        formula,
        Symbol.(spec["rel_varlist"]),
        Symbol(spec["control_cohort"]),
        Symbol(spec["cohort"]),
        Vcov.cluster(Symbol(spec["cluster"]));
        progress_bar = false,
    )
end

function check_reference(path::String)
    spec = TOML.parsefile(path)
    model = build_nlswork_model(spec)

    atol = Float64(get(spec, "atol", 1.0e-4))
    rtol = Float64(get(spec, "rtol", 1.0e-4))
    metric_atol = Float64(get(spec, "metric_atol", 5.0e-3))
    metric_rtol = Float64(get(spec, "metric_rtol", 5.0e-3))

    observed_coefficients = Dict(string(name) => value for (name, value) in zip(coefnames(model), coef(model)))

    missing_coefficients = String[]
    mismatched_coefficients = String[]
    max_abs_diff = 0.0

    for name in spec["rel_varlist"]
        if !haskey(observed_coefficients, name)
            push!(missing_coefficients, name)
            continue
        end

        observed = observed_coefficients[name]
        expected = Float64(spec["coefficients"][name])
        diff = abs(observed - expected)
        max_abs_diff = max(max_abs_diff, diff)

        if !isapprox(observed, expected; atol = atol, rtol = rtol)
            push!(mismatched_coefficients, name)
        end
    end

    missing_metrics = String[]
    mismatched_metrics = String[]
    max_metric_abs_diff = 0.0

    for (metric_name, expected_raw) in get(spec, "metrics", Dict{String, Any}())
        observed = metric_value(model, metric_name)
        expected = Float64(expected_raw)

        if observed isa Integer
            if observed != Int(expected)
                push!(mismatched_metrics, metric_name)
            end
            max_metric_abs_diff = max(max_metric_abs_diff, abs(Float64(observed) - expected))
            continue
        end

        diff = abs(Float64(observed) - expected)
        max_metric_abs_diff = max(max_metric_abs_diff, diff)

        if !isapprox(Float64(observed), expected; atol = metric_atol, rtol = metric_rtol)
            push!(mismatched_metrics, metric_name)
        end
    end

    passed = isempty(missing_coefficients) && isempty(mismatched_coefficients) &&
        isempty(missing_metrics) && isempty(mismatched_metrics)

    return Dict(
        "name" => get(spec, "name", splitext(basename(path))[1]),
        "status" => passed ? "pass" : "fail",
        "checked_coefficients" => length(spec["rel_varlist"]),
        "max_abs_diff" => max_abs_diff,
        "mismatched_coefficients" => mismatched_coefficients,
        "missing_coefficients" => missing_coefficients,
        "checked_metrics" => sort!(collect(keys(get(spec, "metrics", Dict{String, Any}())))),
        "max_metric_abs_diff" => max_metric_abs_diff,
        "mismatched_metrics" => mismatched_metrics,
        "missing_metrics" => missing_metrics,
    )
end

function write_reference_report(path::String, results::Vector{Dict{String, Any}})
    mkpath(dirname(path))

    by_name = Dict{String, Any}()
    for result in results
        name = result["name"]
        by_name[name] = Dict(
            "status" => result["status"],
            "checked_coefficients" => result["checked_coefficients"],
            "checked_metrics" => result["checked_metrics"],
            "max_abs_diff" => result["max_abs_diff"],
            "max_metric_abs_diff" => result["max_metric_abs_diff"],
            "missing_coefficients" => result["missing_coefficients"],
            "mismatched_coefficients" => result["mismatched_coefficients"],
            "missing_metrics" => result["missing_metrics"],
            "mismatched_metrics" => result["mismatched_metrics"],
        )
    end

    payload = Dict(
        "generated_at_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
        "references" => by_name,
    )

    open(path, "w") do io
        TOML.print(io, payload)
    end
end

reference_results = Dict{String, Any}[]
report_path = reference_report_path()

try
    @testset "EventStudyInteracts" begin
        @testset "Smoke test" begin
            df = make_test_panel()
            rel_varlist = [:g_2, :g0, :g1, :g2]
            absorb = [:id, :t]
            formula = term(:Y) ~ sum(fe.(term.(absorb)))

            model = eventreg(
                df,
                formula,
                rel_varlist,
                :never_treat,
                :first_treat,
                Vcov.simple();
                progress_bar = false,
            )

            @test model isa EventStudyInteract
            @test length(coef(model)) == length(rel_varlist)
            @test coefnames(model) == string.(rel_varlist)
            @test size(vcov(model)) == (length(rel_varlist), length(rel_varlist))
            @test nobs(model) == nrow(df)
            @test isfinite(r2(model))
            @test size(confint(model)) == (length(rel_varlist), 2)
            @test !isempty(coeftable(model).rownms)
        end

        @testset "RegressionTables integration" begin
            df = make_test_panel()
            rel_varlist = [:g_2, :g0, :g1, :g2]
            absorb = [:id, :t]
            formula = term(:Y) ~ sum(fe.(term.(absorb)))

            model = eventreg(
                df,
                formula,
                rel_varlist,
                :never_treat,
                :first_treat,
                Vcov.simple();
                progress_bar = false,
            )

            table = regtable(
                model;
                render = AsciiTable(),
                regression_statistics = [Nobs, R2],
            )
            rendered = sprint(print, table)

            @test occursin("Y", rendered)
            @test occursin("g0", rendered)
            @test occursin("R2", rendered)
        end

        @testset "Share covariance uses cluster structure" begin
            share_df, X, e = make_clustered_share_inputs()

            robust_sigma = EventStudyInteracts._share_vcov(share_df, e, X, Vcov.robust())
            cluster1_sigma = EventStudyInteracts._share_vcov(share_df, e, X, Vcov.cluster(:cluster_a))
            cluster2_sigma = EventStudyInteracts._share_vcov(share_df, e, X, Vcov.cluster(:cluster_a, :cluster_b))

            @test robust_sigma ≈ manual_robust_share_vcov(X, e) atol = 1.0e-10 rtol = 1.0e-10
            @test all(isfinite, diag(cluster1_sigma))
            @test all(isfinite, diag(cluster2_sigma))
            @test maximum(abs.(Matrix(cluster1_sigma) .- Matrix(robust_sigma))) > 1.0e-8
            @test maximum(abs.(Matrix(cluster2_sigma) .- Matrix(robust_sigma))) > 1.0e-8
            @test maximum(abs.(Matrix(cluster1_sigma) .- Matrix(cluster2_sigma))) > 1.0e-8
        end

        @testset "Multiway cluster eventreg" begin
            df = make_test_panel()
            df.industry = repeat([1, 1, 2, 2, 3, 3], inner = 5)

            rel_varlist = [:g_2, :g0, :g1, :g2]
            absorb = [:id, :t]
            formula = term(:Y) ~ sum(fe.(term.(absorb)))

            model = eventreg(
                df,
                formula,
                rel_varlist,
                :never_treat,
                :first_treat,
                Vcov.cluster(:id, :industry);
                progress_bar = false,
            )

            @test model isa EventStudyInteract
            @test size(vcov(model)) == (length(rel_varlist), length(rel_varlist))
            @test all(isfinite, diag(vcov(model)))
        end

        @testset "Zero relative-time columns do not poison clustered share vcov" begin
            df = make_test_panel()
            df.g_3 = zeros(Int, nrow(df))

            formula = term(:Y) ~ sum(fe.(term.([:id, :t])))
            baseline = eventreg(
                df,
                formula,
                [:g_2, :g0, :g1, :g2],
                :never_treat,
                :first_treat,
                Vcov.cluster(:id);
                progress_bar = false,
            )
            with_zero = eventreg(
                df,
                formula,
                [:g_3, :g_2, :g0, :g1, :g2],
                :never_treat,
                :first_treat,
                Vcov.cluster(:id);
                progress_bar = false,
            )

            @test coefnames(with_zero) == string.([:g_3, :g_2, :g0, :g1, :g2])
            @test coef(with_zero)[1] == 0.0
            @test isnan(diag(vcov(with_zero))[1])
            @test all(isnan, vcov(with_zero)[1, :])
            @test all(isnan, vcov(with_zero)[:, 1])
            @test coef(with_zero)[2:end] ≈ coef(baseline) atol = 1.0e-10 rtol = 1.0e-10
            @test diag(vcov(with_zero))[2:end] ≈ diag(vcov(baseline)) atol = 1.0e-10 rtol = 1.0e-10
        end

        @testset "Single-cohort equivalence to direct FE regression" begin
            report = run_single_cohort_equivalence()

            for row in report.results
                cmp = row.comparison
                @testset "$(row.label)" begin
                    @test cmp.max_abs_coef_diff < 1.0e-10
                    @test cmp.max_abs_se_diff < 1.0e-10
                    @test cmp.max_abs_vcov_diff < 1.0e-10
                    @test cmp.ate_estimate_diff < 1.0e-10
                    @test cmp.ate_se_diff < 1.0e-10
                end
            end
        end

        @testset "Reference parity" begin
            references = available_reference_files()
            @test !isempty(references)

            for reference_path in references
                result = check_reference(reference_path)
                push!(reference_results, result)

                @testset "$(result["name"])" begin
                    @test result["status"] == "pass"
                    @test isempty(result["missing_coefficients"])
                    @test isempty(result["mismatched_coefficients"])
                    @test isempty(result["missing_metrics"])
                    @test isempty(result["mismatched_metrics"])
                end
            end
        end
    end
finally
    if report_path !== nothing
        resolved_path = isabspath(report_path) ? report_path : joinpath(REPO_ROOT, report_path)
        write_reference_report(resolved_path, reference_results)
    end
end



