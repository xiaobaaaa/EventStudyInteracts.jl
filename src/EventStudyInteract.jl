
##############################################################################
##
## Type EventStudyInteract
##
##############################################################################

struct EventStudyInteract <: RegressionModel
    coef::Vector{Float64}   # Vector of coefficients
    vcov::Matrix{Float64}   # Covariance matrix
    vcov_type::CovarianceEstimator
    nclusters::Union{NamedTuple, Nothing}

    esample::BitVector      # Is the row of the original dataframe part of the estimation sample?
    residuals::Union{AbstractVector, Nothing}
    fe::DataFrame
    fekeys::Vector{Symbol}


    coefnames::Vector       # Name of coefficients
    yname::Union{String, Symbol} # Name of dependent variable
    formula::FormulaTerm        # Original formula
    formula_schema::FormulaTerm # Schema for predict
    contrasts::Dict

    nobs::Int64             # Number of observations
    dof::Int64              # Number parameters estimated - has_intercept
    dof_residual::Int64        # dof used for t-test and F-stat. nobs - degrees of freedoms with simple std

    rss::Float64            # Sum of squared residuals
    tss::Float64            # Total sum of squares
    r2::Float64             # R squared
    adjr2::Float64          # R squared adjusted

    F::Float64              # F statistics
    p::Float64              # p value for the F statistics

    # for FE
    iterations::Union{Int, Nothing}         # Number of iterations
    converged::Union{Bool, Nothing}         # Has the demeaning algorithm converged?
    r2_within::Union{Float64, Nothing}      # within r2 (with fixed effect

    feresult::FixedEffectModel
end

StatsAPI.coef(m::EventStudyInteract) = m.coef
StatsAPI.coefnames(m::EventStudyInteract) = m.coefnames
StatsAPI.responsename(m::EventStudyInteract) = m.yname
StatsAPI.vcov(m::EventStudyInteract) = m.vcov
StatsAPI.nobs(m::EventStudyInteract) = m.nobs
StatsAPI.dof(m::EventStudyInteract) = m.dof
StatsAPI.dof_residual(m::EventStudyInteract) = m.dof_residual
StatsAPI.r2(m::EventStudyInteract) = m.r2
StatsAPI.adjr2(m::EventStudyInteract) = m.adjr2
StatsAPI.islinear(m::EventStudyInteract) = true
StatsAPI.deviance(m::EventStudyInteract) = m.tss
StatsAPI.rss(m::EventStudyInteract) = m.rss
StatsAPI.mss(m::EventStudyInteract) = deviance(m) - rss(m)
StatsAPI.residuals(m::EventStudyInteract) = m.residuals2

function StatsAPI.confint(m::EventStudyInteract; level::Real = 0.95)
    scale = tdistinvcdf(StatsAPI.dof_residual(m), 1 - (1 - level) / 2)
    se = stderror(m)
    hcat(m.coef -  scale * se, m.coef + scale * se)
end


function StatsAPI.coeftable(m::EventStudyInteract; level = 0.95)
    cc = coef(m)
    se = stderror(m)
    coefnms = coefnames(m)
    conf_int = confint(m; level = level)
    # put (intercept) last
    if !isempty(coefnms) && ((coefnms[1] == Symbol("(Intercept)")) || (coefnms[1] == "(Intercept)"))
        newindex = vcat(2:length(cc), 1)
        cc = cc[newindex]
        se = se[newindex]
        conf_int = conf_int[newindex, :]
        coefnms = coefnms[newindex]
    end
    tt = cc ./ se
    CoefTable(
        hcat(cc, se, tt, fdistccdf.(Ref(1), Ref(StatsAPI.dof_residual(m)), abs2.(tt)), conf_int[:, 1:2]),
        ["Estimate","Std.Error","t value", "Pr(>|t|)", "Lower 95%", "Upper 95%" ],
        ["$(coefnms[i])" for i = 1:length(cc)], 4)
end


##############################################################################
##
## Display Result
##
##############################################################################

function title(m::EventStudyInteract)
    return "EventStudyInteract IW estimator"
end

format_scientific(x) = @sprintf("%.3f", x)

function top(m::EventStudyInteract)
    out = [
            "Number of obs" sprint(show, nobs(m), context = :compact => true);
            "Degrees of freedom" sprint(show, dof(m), context = :compact => true);
            "R2" format_scientific(r2(m));
            "R2 Adjusted" format_scientific(adjr2(m));
            "F-Stat" sprint(show, m.F, context = :compact => true);
            "p-value" format_scientific(m.p);
            ]

    out = vcat(out, 
        ["R2 within" format_scientific(m.r2_within);
        "Iterations" sprint(show, m.iterations, context = :compact => true);
            ])

    return out
end


function Base.show(io::IO, m::EventStudyInteract)
    ctitle = title(m)
    ctop = top(m)
    cc = coef(m)
    se = stderror(m)
    yname = responsename(m)
    coefnms = coefnames(m)
    conf_int = confint(m)
    # put (intercept) last
    if !isempty(coefnms) && ((coefnms[1] == Symbol("(Intercept)")) || (coefnms[1] == "(Intercept)"))
        newindex = vcat(2:length(cc), 1)
        cc = cc[newindex]
        se = se[newindex]
        conf_int = conf_int[newindex, :]
        coefnms = coefnms[newindex]
    end
    tt = cc ./ se
    mat = hcat(cc, se, tt, fdistccdf.(Ref(1), Ref(StatsAPI.dof_residual(m)), abs2.(tt)), conf_int[:, 1:2])
    nr, nc = size(mat)
    colnms = ["Estimate","Std.Error","t value", "Pr(>|t|)", "Lower 95%", "Upper 95%"]
    rownms = ["$(coefnms[i])" for i = 1:length(cc)]
    pvc = 4
    # print
    if length(rownms) == 0
        rownms = AbstractString[lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    if length(rownms) > 0
        rnwidth = max(4, maximum(length(nm) for nm in rownms) + 2, length(yname) + 2)
        else
            # if only intercept, rownms is empty collection, so previous would return error
        rnwidth = 4
    end
    rownms = [rpad(nm,rnwidth-1) * "|" for nm in rownms]
    widths = [length(cn)::Int for cn in colnms]
    str = [sprint(show, mat[i,j]; context=:compact => true) for i in 1:nr, j in 1:nc]
    if pvc != 0                         # format the p-values column
        for i in 1:nr
            str[i, pvc] = format_scientific(mat[i, pvc])
        end
    end
    for j in 1:nc
        for i in 1:nr
            lij = length(str[i, j])
            if lij > widths[j]
                widths[j] = lij
            end
        end
    end
    widths .+= 1
    totalwidth = sum(widths) + rnwidth
    if length(ctitle) > 0
        halfwidth = div(totalwidth - length(ctitle), 2)
        println(io, " " ^ halfwidth * string(ctitle) * " " ^ halfwidth)
    end
    if length(ctop) > 0
        for i in 1:size(ctop, 1)
            ctop[i, 1] = ctop[i, 1] * ":"
        end
        println(io, "=" ^totalwidth)
        halfwidth = div(totalwidth, 2) - 1
        interwidth = 2 +  mod(totalwidth, 2)
        for i in 1:(div(size(ctop, 1) - 1, 2)+1)
            print(io, ctop[2*i-1, 1])
            print(io, lpad(ctop[2*i-1, 2], halfwidth - length(ctop[2*i-1, 1])))
            print(io, " " ^interwidth)
            if size(ctop, 1) >= 2*i
                print(io, ctop[2*i, 1])
                print(io, lpad(ctop[2*i, 2], halfwidth - length(ctop[2*i, 1])))
            end
            println(io)
        end
    end
    println(io,"=" ^totalwidth)
    println(io, rpad(string(yname), rnwidth-1) * "|" *
            join([lpad(string(colnms[i]), widths[i]) for i = 1:nc], ""))
    println(io,"-" ^totalwidth)
    for i in 1:nr
        print(io, rownms[i])
        for j in 1:nc
            print(io, lpad(str[i,j],widths[j]))
        end
        println(io)
    end
    println(io,"=" ^totalwidth)
end