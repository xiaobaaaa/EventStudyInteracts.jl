# EventStudyInteracts.jl
This package is a Julia replication of the Stata package [`eventstudyinteract` ](https://github.com/lsun20/EventStudyInteract)provided by [Sun and Abraham (2021)](https://www.sciencedirect.com/science/article/abs/pii/S030440762030378X).

The estimation process is identical to that of eventstudyinteract, and this package used [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl) to implement the Stata command [`reghdfe`](https://github.com/sergiocorreia/reghdfe) functionality.

As introduced in the Stata package [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract) provided by [Sun and Abraham (2021)](https://www.sciencedirect.com/science/article/abs/pii/S030440762030378X), the function of this package is:

> [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract) implements the interaction weighted (IW) estimator and constructs pointwise confidence interval for the estimation of dynamic treatment effects. To estimate the dynamic effects of an absorbing treatment, researchers often use two-way fixed effects (TWFE) regressions that include leads and lags of the treatment (event study specification). Units are categorized into different cohorts based on their initial treatment timing. Sun and Abraham (2020) proposes this estimator as an alternative to the TWFE regression in the presence of treatment effects heterogeneous across cohorts. Under treatment effects heterogeneity, the TWFE regression can result in estimates with uninterpretable weights, which can be assessed by the Stata module eventstudyweights. The IW estimator is implemented in three steps. First, estimate the interacted regression with reghdfe, where the interactions are between relative time indicators and cohort indicators. Second, estimate the cohort shares underlying each relative time. Third, take the weighted average of estimates from the first step, with weights set to the estimated cohort shares.

## Introduction

- Julia is faster than Stata, and [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl) is much faster than reghdfe.
  
- I am a beginner in programming and the code is a complete translation of the Stata package [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract). I used newbing to help me with programming, and it is very helpful for beginners.
  
- I used [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl) directly for estimation and I am grateful to the author of this package. Other parts of this package, such as reporting of estimation results, also refer to [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl).
  
- The Stata package [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract) is very flexible, and one important feature is that it can compare event study estimates for subsamples, i.e. perform heterogeneous estimation. I did not find this feature in the existing Julia package [DiffinDiffs.jl](https://github.com/JuliaDiffinDiffs/DiffinDiffs.jl) or the  R packages [`fixest`](https://lrberge.github.io/fixest/). I think the design of the Stata package [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract) is more user-friendly.
  
- This package fully supports the features of [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl), such as using GPU for estimation. The syntax of this package is similar to that of [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl).
  

## Installation

The package is registered in the [`General`](https://github.com/JuliaRegistries/General) registry and so can be installed at the REPL with `] add EventStudyInteracts`.

## Syntax

```julia
result =  reg(df,    
    formula,
    rel_varlist::Vector{Symbol},
    control_cohort::Symbol,
    cohort::Symbol,
    vcov,
    contrasts = contrasts,
    weights = weights,
    save = save,
    method = method,
    nthreads = nthreads,
    double_precision = double_precision,
    tol = tol,
    maxiter = maxiter,
    drop_singletons = drop_singletons,
    progress_bar = progress_bar,
    dof_add = dof_add,
    subset = subset)
```

The syntax is similar to [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl) （You must use a version higher than v1.9.1），but users need to specify separately `rel_varlist::Vector{Symbol}`, `control_cohort::Symbol`, and `cohort::Symbol`. The `rel_varlist` is the list of relative time indicators as you would have included in the canonical two-way fixed effects regression. The `rel_varlist` must be defined as a `Vector{Symbol}`.

```
rel_varlist = Symbol[:g_4, :g_3, :g_2, :g0, :g1, :g2, :3]
```

See illustration for an example of generating these relative time indicators. The relative time indicators should take the value of zero for never treated units. Users should shape their dataset to a long format where each observation is at the unit-time level. See illustration for an example of specifying the syntax. The syntax is similar to [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl) in specifying fixed effects and the type of standard error reported.

Furthermore, it also requires the user to specify the cohort categories as well as which cohort is the control cohort. Note that [Sun and Abraham (2021)](https://www.sciencedirect.com/science/article/abs/pii/S030440762030378X) only establishes the validity of the IW estimators for balanced panel data without covariates. eventstudyinteract does not impose a balanced panel by dropping units with observations that are missing in any iod. Instead, it evaluates the IW estimators for unbalanced panel data by estimating the interacted regression using all observations.

The `control_cohort` must be defined as a `Symbol`, which is a binary variable that corresponds to the control cohort, which can be never-treated units or last-treated units. If using **last-treated** unit as control cohort, **exclude the time periods when the last cohort receives treatment**.

The `cohort` must be defined as a `Symbol`, which is a categorical variable that corresponds to the initial treatment timing of each unit. If there are units that receive multiple treatments, [Sun and Abraham (2021)](https://www.sciencedirect.com/science/article/abs/pii/S030440762030378X) defines the initial treatment timing to be based on the first treatment. This categorical variable should be set to be **missing** for **never treated** units.

## Examples

- You can download the `nlswork.dta` data from the repository, which is the example data used by [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract), and the method of generating variables is the same as that of [`eventstudyinteract`](https://github.com/lsun20/EventStudyInteract).

```julia
using DataFrames
using FixedEffectModels
using EventStudyInteract
using ReadStatTables
using RegressionTables
using CUDA
using Plots

# Load the 1968 extract of the National Longitudinal Survey of Young Women and Mature Women.
df =DataFrame(readstat(EventStudyInteracts\dataset\nlswork.dta"))

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

m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1, vcov1)

#                    EventStudyInteract IW estimator
# ========================================================================
# Number of obs:                27979  Degrees of freedom:             166
# R2:                           0.666  R2 Adjusted:                  0.664
# F-Stat:                     2.81691  p-value:                      0.000
# R2 within:                    0.030  Iterations:                      11
# ========================================================================
# ln_wage |   Estimate Std.Error   t value Pr(>|t|)   Lower 95%  Upper 95%
# ------------------------------------------------------------------------
# g_19    |        0.0       0.0       NaN      NaN         0.0        0.0
# g_18    | -0.0609283 0.0513893  -1.18562    0.236   -0.161679  0.0398221
# g_17    | -0.0489819 0.0432146  -1.13346    0.257   -0.133706  0.0357418
# g_16    | -0.0680484 0.0495311  -1.37385    0.170   -0.165156  0.0290589
# g_15    |  -0.135218 0.0448805  -3.01284    0.003   -0.223208 -0.0472281
# g_14    |  -0.103295 0.0443376  -2.32975    0.020   -0.190221   -0.01637
# g_13    | -0.0684827 0.0416284   -1.6451    0.100   -0.150097  0.0131311
# g_12    |  -0.132015 0.0385567  -3.42392    0.001   -0.207606 -0.0564232
# g_11    |   -0.15239 0.0376567  -4.04683    0.000   -0.226217  -0.078563
# g_10    |   -0.11842 0.0291013  -4.06925    0.000   -0.175474 -0.0613663
# g_9     |  -0.172996 0.0330077  -5.24109    0.000   -0.237709  -0.108284
# g_8     |  -0.116192  0.030757  -3.77773    0.000   -0.176492 -0.0558915
# g_7     |  -0.140938 0.0287493  -4.90233    0.000   -0.197302 -0.0845745
# g_6     |  -0.133223 0.0369878  -3.60181    0.000   -0.205739 -0.0607073
# g_5     |  -0.131416 0.0238826  -5.50261    0.000   -0.178239 -0.0845938
# g_4     |  -0.150598 0.0330246  -4.56016    0.000   -0.215344 -0.0858518
# g_3     | -0.0559142 0.0202356  -2.76316    0.006  -0.0955869 -0.0162416
# g_2     | -0.0565837 0.0191366  -2.95683    0.003  -0.0941016 -0.0190657
# g0      |  0.0515795 0.0146284   3.52598    0.000      0.0229  0.0802591
# g1      |  0.0961542 0.0194692   4.93879    0.000   0.0579842   0.134324
# g2      |  0.0677024 0.0182556   3.70858    0.000   0.0319116   0.103493
# g3      |  0.0663372  0.018431   3.59923    0.000   0.0302027   0.102472
# g4      |   0.132223 0.0276378   4.78414    0.000   0.0780381   0.186408
# g5      |  0.0565864 0.0182183   3.10603    0.002   0.0208689  0.0923039
# g6      |   0.117224 0.0215077    5.4503    0.000   0.0750571    0.15939
# g7      |   0.072123 0.0221207   3.26042    0.001   0.0287545   0.115491
# g8      |  0.0666257 0.0193959   3.43504    0.001   0.0285993   0.104652
# g9      |  0.0929835 0.0349948   2.65707    0.008    0.024375   0.161592
# g10     |  0.0675955 0.0237323   2.84825    0.004   0.0210676   0.114124
# g11     |   0.106953 0.0212674   5.02896    0.000   0.0652575   0.148649
# g12     |  0.0647985 0.0365021    1.7752    0.076 -0.00676515   0.136362
# g13     |   0.134787 0.0412699   3.26599    0.001   0.0538762   0.215698
# g14     |   0.155875 0.0501801   3.10631    0.002   0.0574951   0.254255
# g15     |  0.0665542 0.0426601   1.56011    0.119  -0.0170823   0.150191
# g16     |   0.181412 0.0511286   3.54815    0.000   0.0811727   0.281651
# g17     |  0.0233864 0.0438817  0.532943    0.594  -0.0626452   0.109418
# g18     | -0.0477967 0.0750543 -0.636828    0.524   -0.194943  0.0993498
# ========================================================================
```

- **Event study plots**. We may define a ``coefplot`` function for an event study plot.

```julia
function mycoefplot(m)
        n = coefnames(m)
        vals = coef(m)
        errors = stderror(m)
        Plots.scatter(
            n,
            vals,
            legend = false,
            yerror = 1.96 .* errors,
            #title = "Coefficient plot",
            xtickfont = font(10, "Times"),
            ytickfont = font(10, "Times"),
            titlefont = font(14, "Times"),
            # markersize = 2
    )
    Plots.vline!([3.5,],linestyle = :dash,linecolor=:red)
    Plots.hline!([0,],linestyle = :dash,linecolor=:black,linewidth=:2)
end

plot1 = mycoefplot(m1)
```

- **Binning**. Pre-treatment effects seem relatively constant, which might suggest binning the many leads.  TODO: current implementation of bins does not follow Sun and Abraham (2020) exactly due to coding challenge.  But it is valid if effects in the bin are constant for each cohort.

```julia
df.g_l4 = Int.(coalesce.(df.ry .<= -4, false))

rel_varlist2 = Symbol[:g_l4]

for i in 3:-1:2
    push!(rel_varlist2, Symbol("g_"*string(i)))
end

for i in 0:18
    push!(rel_varlist2, Symbol("g"*string(i)))
end

m2 = eventreg(df, formula1, rel_varlist2, control_cohort1, cohort1, vcov1)
```

- **Using the last treated as the control cohort**. Alternatively, we can take the control cohort to be individuals that were unionized last.

```julia
df.last_union = Int.(coalesce.(df.first_union .== 88, false))
control_cohort2 = :last_union
# If using the last-treated cohort as the control, be sure to restrict the analysis sample to exclude the never-treated (if any) and to be before the treated periods for the last-treated cohort.

m3 = eventreg(df[(.!ismissing.(df.first_union)) .& (df.year .< 88), :], formula1, rel_varlist2, control_cohort1, cohort1, vcov1)
```

- **Aggregating event study estimates**. In Stata one can use `lincom` looks for coefficients and variance covariance matrix stored in e(b) and e(V). There is no package like ``lincom`` in Julia. I write another package [`LinComs.jl`]https://github.com/xiaobaaaa/LinComs.jl) to do this. You can refer to the following code after estimating the results.

```julia
rel_varlist1[20:38]
expr = Expr(:call, :/, Expr(:call, :+, rel_varlist1...), length(rel_varlist1))
# Or,
expr = :((g_20 + g_19 + g_18 + g_17 + g_16 + g_15 + g_14 + g_13 + g_12 + g_11 + g_10 + g_9 + g_8 + g_7 + g_6 + g_5 + g_4 + g_3 + g_2 + g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + g16 + g17 + g18) / 38)
# And then,
lincom(m1,expr)
# Returns,
                                      lincom                                      
==================================================================================
ln_wage            |    Estimate Std.Error   t value Pr(>|t|)  Lower 95% Upper 95%
----------------------------------------------------------------------------------
Linear Combination | -0.00369359 0.0131072 -0.281798    0.778 -0.0293907 0.0220035
==================================================================================
```

The results obtained from LinComs.jl can also be used with RegressionTables.jl to generate Regression Tables.
  
- **Compare event study estimates for subsamples**. Suppose we want to compare the average effect over the first five years of joining the union between college graduates and non-college graduates. We can first estimate their separate effects by interacting the relative time indicators with icator of college graduates:
  

```julia
for k in 18:-1:2
    df[!, Symbol("g_", k, "_collgrad0")] = ifelse.(coalesce.(df.ry .== -k, false) .& coalesce.(df.collgrad .== 0, false), 1, 0)
end
for k in 0:18
    df[!, Symbol("g", k, "_collgrad0")] = ifelse.(coalesce.(df.ry .== k, false) .& coalesce.(df.collgrad .== 0, false), 1, 0)
end
df.g_l4_collgrad0 = ifelse.(coalesce.(df.ry .<= -4, false) .& coalesce.(df.collgrad .== 0, false), 1, 0)
for k in 18:-1:2
    df[!, Symbol("g_", k, "_collgrad1")] = ifelse.(coalesce.(df.ry .== -k, false) .& coalesce.(df.collgrad .== 1, false), 1, 0)
end
for k in 0:18
    df[!, Symbol("g", k, "_collgrad1")] = ifelse.(coalesce.(df.ry .== k, false) .& coalesce.(df.collgrad .== 1, false), 1, 0)
end
df.g_l4_collgrad1 = ifelse.(coalesce.(df.ry .<= -4, false) .& coalesce.(df.collgrad .== 1, false), 1, 0)

rel_varlist3 = Symbol[:g_l4_collgrad0]

for i in 3:-1:2
    push!(rel_varlist3, Symbol("g_"*string(i)*"_collgrad0"))
end

for i in 0:18
    push!(rel_varlist3, Symbol("g"*string(i)*"_collgrad0"))
end

push!(rel_varlist3, :g_l4_collgrad1)
for i in 3:-1:2
    push!(rel_varlist3, Symbol("g_"*string(i)*"_collgrad1"))
end

for i in 0:18
    push!(rel_varlist3, Symbol("g"*string(i)*"_collgrad1"))
end

m4 = eventreg(df, formula1, rel_varlist3, control_cohort1, cohort1, vcov1)
```

## Output

`eventreg` returns a light object like ``FixedEffectModels.reg``. It is composed of the vector of coefficients & the covariance matrix (use coef, coefnames, vcov on the output of reg) a boolean vector reporting rows used in the estimation a set of scalars (number of observations, the degree of freedoms, r2, etc).

We can use RegressionTables.jl to get the Regression Tables.

```julia
regtable(m1, m2, m3,m4; renderSettings = asciiOutput(), estim_decoration=make_estim_decorator([0.01, 0.05, 0.1]))

regtable(m1, m2, m3,m4; renderSettings = htmlOutput("base_reg.html"), estim_decoration=make_estim_decorator([0.01, 0.05, 0.1]))
```

## Use other features of FixedEffectModels.jl

We can also use other features of FixedEffectModels.jl, such as using GPU for estimation.

```julia
eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1, vcov1,  method = :gpu, double_precision = false)
```

## Another example

Refer to [Difference-in-Difference (DiD) | DiD](https://asjadnaqvi.github.io/DiD/).

```julia
# Generate sample data
using Random

units = 30
start = 1
finish = 60
nlvls = 5
time = finish - start + 1
obsv = units * time

df = DataFrame(id = repeat(1:units, inner=time), t = repeat(start:finish, outer=units))

Random.seed!(20211222)

df.Y .= 0 # outcome variable
df.D .= 0 # intervention variable
df.cohort .= repeat([rand(0:nlvls) for i in 1:units], inner=time) # treatment cohort

effect = [rand(2:10) for i in unique(df.cohort)]
first_treat = [rand(start:(finish + 20)) for i in unique(df.cohort)]

df.effect = similar(df.cohort)
df.first_treat = similar(df.cohort, Union{Int64, Missing})

for i in 0:nlvls
    df.effect[df.cohort .== i] .= effect[i+1]
    df.first_treat[df.cohort .== i] .= first_treat[i+1]
end

df.first_treat[df.first_treat .> finish] .= missing


df.D = ifelse.(coalesce.(df.t .>= df.first_treat, false), 1, 0)
df.rel_time .= df.t .- df.first_treat
df.Y = df.id .+ df.t .+ ifelse.(df.D .== 1, df.effect .* df.rel_time, 0) .+ randn(nrow(df))

plot(df.t, df.Y, label=false)

tab1(df,:cohort)

df.never_treat = ismissing.(df.first_treat)

tab1(df,:rel_time)

# We will consider the dynamic effect of union status on income. We first generate these relative time indicators, and leave out the distant leads due to few observations. Implicitly this assumes that effects outside the lead windows are zero.

relmin = abs(minimum(skipmissing(df.rel_time)))
relmax = abs(maximum(skipmissing(df.rel_time)))

for k in relmin:-1:2
    df[!, Symbol("g_$k")] = Int.(coalesce.(df.rel_time .== -k, false))
end

for k in 0:relmax
    df[!, Symbol("g$k")] = Int.(coalesce.(df.rel_time .== k, false))
end

absorb = [:id,:t]

formula1 = term(:Y) ~ sum(fe.(term.(absorb)))

rel_varlist1 = Symbol[]

for i in relmin:-1:2
    push!(rel_varlist1, Symbol("g_"*string(i)))
end

for i in 0:relmax
    push!(rel_varlist1, Symbol("g"*string(i)))
end

control_cohort1 = :never_treat

cohort1 = :first_treat

m5 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1)

plot2 = mycoefplot(m5)
```

## Performance
EventStudyInteracts.jl is **3.69** times faster than Stata when controlling only for id and year fixed effects with a sample size of 2 million observations. EventStudyInteracts.jl is **4.75** times faster than Stata when controlling for an additional two fixed effects. As the sample size increases and more fixed effects are controlled for, performance improvement should be even greater. Furthermore, EventStudyInteracts.jl can also utilize CUDA. Many thanks to FixedEffectModels.jl for its outstanding performance!

- Julia Codes
  
  ```julia
  
  using CSV, DataFrames, Random, Plots, FixedEffectModels, EventStudyInteracts
  
  units = 100000
  start = 1
  finish = 20
  nlvls = 5
  time = finish - start + 1
  obsv = units * time
  K=100
  id1 = rand(1:(obsv/K), obsv)
  id2 = rand(1:K, obsv)
  
  df = DataFrame(id = repeat(1:units, inner=time), t = repeat(start:finish, outer=units), id1 = id1, id2 = id2)
  
  Random.seed!(20211222)
  
  df.Y .= 0 # outcome variable
  df.D .= 0 # intervention variable
  df.cohort .= repeat([rand(0:nlvls) for i in 1:units], inner=time) # treatment cohort
  
  effect = [rand(2:10) for i in unique(df.cohort)]
  first_treat = [rand(start:(finish + 20)) for i in unique(df.cohort)]
  
  df.effect = similar(df.cohort)
  df.first_treat = similar(df.cohort, Union{Int64, Missing})
  
  for i in 0:nlvls
      df.effect[df.cohort .== i] .= effect[i+1]
      df.first_treat[df.cohort .== i] .= first_treat[i+1]
  end
  
  df.first_treat[df.first_treat .> finish] .= missing
  
  
  df.D = ifelse.(coalesce.(df.t .>= df.first_treat, false), 1, 0)
  df.rel_time .= df.t .- df.first_treat
  df.Y = df.t .+ ifelse.(df.D .== 1, df.effect .* df.rel_time, 0) .+ randn(nrow(df))
  
  plot(df.t, df.Y, label=false)
  
  
  df.never_treat = ismissing.(df.first_treat)
  
  # We will consider the dynamic effect of union status on income. We first generate these relative time indicators, and leave out the distant leads due to few observations. Implicitly this assumes that effects outside the lead windows are zero.
  
  relmin = abs(minimum(skipmissing(df.rel_time)))
  relmax = abs(maximum(skipmissing(df.rel_time)))
  
  for k in relmin:-1:2
      df[!, Symbol("g_$k")] = Int.(coalesce.(df.rel_time .== -k, false))
  end
  
  for k in 0:relmax
      df[!, Symbol("g$k")] = Int.(coalesce.(df.rel_time .== k, false))
  end
  
  absorb = [:id,:t]
  
  formula1 = term(:Y) ~ sum(fe.(term.(absorb)))
  
  rel_varlist1 = Symbol[]
  
  for i in relmin:-1:2
      push!(rel_varlist1, Symbol("g_"*string(i)))
  end
  
  for i in 0:relmax
      push!(rel_varlist1, Symbol("g"*string(i)))
  end
  
  control_cohort1 = :never_treat
  
  cohort1 = :first_treat
  
  CSV.write("bench.csv", df)
  
  m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1)
  @time m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1)
  # 76.965495 seconds (395.81 k allocations: 218.291 GiB, 23.78% gc time)
  absorb = [:id,:id1,:id2,:t]
  formula2 = term(:Y) ~ sum(fe.(term.(absorb)))
  @time m2 = eventreg(df, formula2, rel_varlist1, control_cohort1, cohort1)
  # 88.459680 seconds (3.09 M allocations: 218.534 GiB, 21.19% gc time, 0.86% compilation time)
  ```
  
- Stata Codes
  
  ```stata
  import delimited "~\EventStudyInteracts\benchmark\bench.csv", clear
  timer clear
  set rmsg on
  gen never_treat1 = never_treat=="true"
  eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat1) absorb(i.id i.t)
  r; t=284.17
  eventstudyinteract y g_* g0-g16, cohort(first_treat) control_cohort(never_treat1) absorb(i.id i.id1 i.id2 i.t)
  r; t=420.61
  ```
## References

[GitHub - FixedEffects/FixedEffectModels.jl: Fast Estimation of Linear Models with IV and High Dimensional Categorical Variables](https://github.com/FixedEffects/FixedEffectModels.jl)

[lsun20/EventStudyInteract (github.com)](https://github.com/lsun20/EventStudyInteract)

[Difference-in-Difference (DiD) | DiD](https://asjadnaqvi.github.io/DiD/)
