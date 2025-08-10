using Pkg
Pkg.add(["EventStudyInteracts", "DataFrames", "ReadStatTables", "CUDA", "FixedEffectModels", "Plots",  "RegressionTables"])


Pkg.add("Revise")  
using Revise


Pkg.develop(path=raw"C:\D_Drive\我的github项目\EventStudyInteracts")  # 注册本地包
using EventStudyInteracts
println(pathof(EventStudyInteracts))

using DataFrames
using FixedEffectModels
using EventStudyInteracts
using ReadStatTables
using RegressionTables
using CUDA
using Plots

# Load the 1968 extract of the National Longitudinal Survey of Young Women and Mature Women.
df =DataFrame(readstat(raw"C:\Users\94933\.julia\packages\EventStudyInteracts\5NSeP\dataset\nlswork.dta"))

# Code the cohort categorical variable based on when the individual first joined the union, which will be inputted in cohort(varname).

df.union_year = ifelse.(coalesce.(df.union, false) .== 1, df.year, missing)
# This is really frustrating. In Stata, the simple code bysort idcode: egen first_union = min(union_year) can achieve this. Why is Julia so complicated?
function min_skipmissing(x)
    v = collect(skipmissing(x))
    return isempty(v) ? missing : minimum(v)
end

transform!(groupby(df, :idcode), :union_year => min_skipmissing => :first_union)

select!(df, Not(:union_year))

# Code the relative time categorical variable.
df.ry = df.year - df.first_union

# For the first example, we take the control cohort to be individuals that never unionized.
df.never_union = ismissing.(df.first_union)

# Check if there is a sufficient number of treated units for each relative time. With very few units it might be better to bin the relative times and assume constant treatment effects within the bin.
function tab1(df::DataFrame, catevar::Symbol)
    gdf = groupby(df, catevar)
    result = combine(gdf, nrow => :freq)
    result.percent = result.freq / nrow(df) * 100
    sort!(result, catevar)
    result.cum = cumsum(result.percent)
    return result
end

tab1(df,:ry)
# We will consider the dynamic effect of union status on income. We first generate these relative time indicators, and leave out the distant leads due to few observations. Implicitly this assumes that effects outside the lead windows are zero.
relmin = abs(minimum(skipmissing(df.ry)))
relmax = abs(maximum(skipmissing(df.ry)))

for k in relmin:-1:2
    df[!, Symbol("g_$k")] = Int.(coalesce.(df.ry .== -k, false))
end

for k in 0:relmax
    df[!, Symbol("g$k")] = Int.(coalesce.(df.ry .== k, false))
end

absorb = [:idcode,:year]

formula1 = term(:ln_wage) ~ term(:south) + sum(fe.(term.(absorb)))

# To test whether the package can run successfully when there are collinear variables.
df.g_19 .= 0

# The rel_varlist must be defined as a Vector{Symbol}. In the example below, when defining an empty Vector, you cannot use rel_varlist1 = [].

rel_varlist1 = Symbol[:g_4, :g_5]

for i in relmin:-1:2
    push!(rel_varlist1, Symbol("g_"*string(i)))
end

for i in 0:relmax
    push!(rel_varlist1, Symbol("g"*string(i)))
end

control_cohort1 = :never_union

cohort1 = :first_union

vcov1 = Vcov.cluster(:idcode)

#  first_stage::Bool = true: Should the first-stage F-stat and p-value be computed?
pathof(EventStudyInteracts)


Pkg.develop(path=raw"C:\D_Drive\我的github项目\EventStudyInteracts")  # 注册本地包
using EventStudyInteracts
println(pathof(EventStudyInteracts))

m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1, vcov1)


