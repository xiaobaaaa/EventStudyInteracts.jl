
module EventStudyInteracts

# slows down tss
#if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
#	@eval Base.Experimental.@optlevel 1
#end

using DataFrames
using FixedEffectModels
using LinearAlgebra
using Printf
using Reexport
using Statistics
using StatsAPI
using StatsBase
using StatsFuns
@reexport using StatsModels
using Tables
using Vcov

include("EventStudyInteract.jl")
include("eventfit.jl")
include("ee.jl")
include("avar.jl")

# Export from StatsBase
export coef, coefnames, coeftable, responsename, vcov, stderror, nobs, dof_residual, r2, adjr2, islinear, deviance, rss, mss, confint, predict, residuals, fit


export eventreg,
EventStudyInteract,
Vcov
end 
