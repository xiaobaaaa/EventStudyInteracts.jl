cd(pwd())
Pkg.activate("./EventStudyInteracts")
using EventStudyInteracts

using DataFrames, Random, Plots, FixedEffectModels
using CSV
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


CSV.write("./EventStudyInteracts/benchmark/bench.csv", df)

m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1)
@time m1 = eventreg(df, formula1, rel_varlist1, control_cohort1, cohort1)
# 76.965495 seconds (395.81 k allocations: 218.291 GiB, 23.78% gc time)
absorb = [:id,:id1,:id2,:t]
formula2 = term(:Y) ~ sum(fe.(term.(absorb)))
@time m2 = eventreg(df, formula2, rel_varlist1, control_cohort1, cohort1)
# 88.459680 seconds (3.09 M allocations: 218.534 GiB, 21.19% gc time, 0.86% compilation time)


