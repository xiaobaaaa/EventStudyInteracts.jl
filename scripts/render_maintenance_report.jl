using Dates
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

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

function normalize_status(status::AbstractString)
    return status == "success" ? "pass" : status
end

function compat_tracks_latest(compat::AbstractString, latest_tag::AbstractString)
    latest_tag == "unknown" && return false
    version_text = replace(latest_tag, 'v' => "")
    version = try
        VersionNumber(version_text)
    catch
        return false
    end

    major_minor = string(version.major, ".", version.minor)
    return occursin(major_minor, compat) || compat == string(version.major)
end

function load_reference_results(path::String)
    isfile(path) || return Dict{String, Any}()
    parsed = TOML.parsefile(path)
    return get(parsed, "references", Dict{String, Any}())
end

function compute_risk_level(lts_status::String, release_status::String, reference_results, compat_ok::Bool)
    reference_failed = any(get(result, "status", "fail") != "pass" for result in values(reference_results))
    if lts_status != "pass" || release_status != "pass" || reference_failed
        return "high"
    elseif !compat_ok
        return "medium"
    else
        return "low"
    end
end

function push_section(lines::Vector{String}, title::String)
    push!(lines, "## $title")
end

options = parse_cli_options(ARGS)
output_path = get(options, "output", nothing)
output_path === nothing && error("Missing required argument --output=...")
output_path = isabspath(output_path) ? output_path : joinpath(REPO_ROOT, output_path)

project = TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))
compat_fem = project["compat"]["FixedEffectModels"]

reference_report_path = get(options, "reference-report", "")
reference_report_path = isabspath(reference_report_path) ? reference_report_path : joinpath(REPO_ROOT, reference_report_path)
reference_results = load_reference_results(reference_report_path)

fem_latest_tag = get(options, "fem-latest-tag", "unknown")
compat_ok = compat_tracks_latest(compat_fem, fem_latest_tag)

lts_status = normalize_status(get(options, "lts-status", "unknown"))
release_status = normalize_status(get(options, "release-status", "unknown"))
git_branch = get(options, "git-branch", "unknown")
git_dirty = get(options, "git-dirty", "unknown")
compathelper_priv_present = lowercase(get(options, "compathelper-priv-present", "false")) == "true"
boss_reference_configured = isfile(joinpath(REPO_ROOT, "reference", "boss_reference.toml"))

workflow_files = [
    ("CI", isfile(joinpath(REPO_ROOT, ".github", "workflows", "ci.yml"))),
    ("CompatHelper", isfile(joinpath(REPO_ROOT, ".github", "workflows", "CompatHelper.yml"))),
    ("TagBot", isfile(joinpath(REPO_ROOT, ".github", "workflows", "TagBot.yml"))),
    ("Maintenance", isfile(joinpath(REPO_ROOT, ".github", "workflows", "maintenance.yml"))),
]

risk_level = compute_risk_level(lts_status, release_status, reference_results, compat_ok)

lines = String[]
push!(lines, "# EventStudyInteracts Maintenance Report")
push!(lines, "")
push!(lines, "Generated at: $(Dates.format(now(UTC), dateformat"yyyy-mm-dd HH:MM:SS")) UTC")
push!(lines, "")

push_section(lines, "Status Summary")
push!(lines, "")
push!(lines, "- Risk level: $risk_level")
push!(lines, "- Git branch: $git_branch")
push!(lines, "- Working tree dirty: $git_dirty")
push!(lines, "- Julia 1.10 LTS maintenance test run: $lts_status")
push!(lines, "- Julia stable import smoke test: $release_status")
push!(lines, "")

push_section(lines, "FixedEffectModels Drift")
push!(lines, "")
push!(lines, "- Current compat in `Project.toml`: `$compat_fem`")
push!(lines, "- Latest upstream release: `$fem_latest_tag`")
push!(lines, "- Compat covers latest release: $(compat_ok ? "yes" : "no")")
push!(lines, "")

push_section(lines, "Reference Parity")
push!(lines, "")
if isempty(reference_results)
    push!(lines, "- No machine-readable reference report was produced in this run.")
else
    for name in sort(collect(keys(reference_results)))
        result = reference_results[name]
        push!(lines, "- `$name`: $(result["status"]) | coeff max abs diff = $(result["max_abs_diff"]) | metric max abs diff = $(result["max_metric_abs_diff"])")
        if !isempty(result["mismatched_coefficients"])
            push!(lines, "- `$name` mismatched coefficients: $(join(result["mismatched_coefficients"], ", "))")
        end
        if !isempty(result["mismatched_metrics"])
            push!(lines, "- `$name` mismatched metrics: $(join(result["mismatched_metrics"], ", "))")
        end
    end
end
push!(lines, "- Boss reference configured: $(boss_reference_configured ? "yes" : "no")")
push!(lines, "")

push_section(lines, "Workflow Coverage")
push!(lines, "")
for (name, present) in workflow_files
    push!(lines, "- $name workflow present: $(present ? "yes" : "no")")
end
push!(lines, "- `COMPATHELPER_PRIV` available: $(compathelper_priv_present ? "yes" : "no")")
push!(lines, "")

push_section(lines, "Recommended Next Steps")
push!(lines, "")
next_steps = String[]
if !compat_ok
    push!(next_steps, "Start the EventStudyInteracts compatibility update for FixedEffectModels $fem_latest_tag.")
end
if lts_status != "pass"
    push!(next_steps, "Investigate the Julia 1.10 LTS test failure before changing package compat.")
end
if release_status != "pass"
    push!(next_steps, "Investigate the Julia stable import failure to keep forward-compatibility visible.")
end
if !boss_reference_configured
    push!(next_steps, "Add `reference/boss_reference.toml` if you want the maintenance workflow to compare against your boss's saved output.")
end
if !compathelper_priv_present
    push!(next_steps, "Add the `COMPATHELPER_PRIV` repository secret so CompatHelper can open update PRs automatically.")
end
if isempty(next_steps)
    push!(next_steps, "No immediate maintenance action is required.")
end

for step in next_steps
    push!(lines, "- $step")
end

mkpath(dirname(output_path))
open(output_path, "w") do io
    write(io, join(lines, "\n"))
    write(io, "\n")
end
