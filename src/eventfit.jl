function eventreg(@nospecialize(df),     
    @nospecialize(formula::FormulaTerm),
    @nospecialize(rel_varlist::Vector{Symbol}),
    @nospecialize(control_cohort::Symbol),
    @nospecialize(cohort::Symbol),
    @nospecialize(vcov::CovarianceEstimator = Vcov.simple());
    @nospecialize(contrasts::Dict = Dict{Symbol, Any}()),
    @nospecialize(weights::Union{Symbol, Nothing} = nothing),
    @nospecialize(save::Union{Bool, Symbol} = :none),
    @nospecialize(method::Symbol = :cpu),
    @nospecialize(nthreads::Integer = method == :cpu ? Threads.nthreads() : 256),
    @nospecialize(double_precision::Bool = true),
    @nospecialize(tol::Real = 1e-6),
    @nospecialize(maxiter::Integer = 10000),
    @nospecialize(drop_singletons::Bool = true),
    @nospecialize(progress_bar::Bool = true),
    @nospecialize(subset::Union{Nothing, AbstractVector} = nothing))
    StatsAPI.fit(EventStudyInteract, formula, df, rel_varlist, control_cohort, cohort, vcov; contrasts = contrasts, weights = weights, save = save, method = method, nthreads = nthreads, double_precision = double_precision, tol = tol, maxiter = maxiter, drop_singletons = drop_singletons, progress_bar = progress_bar, subset = subset)
end
  
function StatsAPI.fit(::Type{EventStudyInteract},     
    @nospecialize(formula::FormulaTerm),
    @nospecialize(df),
    @nospecialize(rel_varlist::Vector{Symbol}),
    @nospecialize(control_cohort::Symbol),
    @nospecialize(cohort::Symbol),
    @nospecialize(vcov::CovarianceEstimator = Vcov.simple());
    @nospecialize(contrasts::Dict = Dict{Symbol, Any}()),
    @nospecialize(weights::Union{Symbol, Nothing} = nothing),
    @nospecialize(save::Union{Bool, Symbol} = :none),
    @nospecialize(method::Symbol = :cpu),
    @nospecialize(nthreads::Integer = method == :cpu ? Threads.nthreads() : 256),
    @nospecialize(double_precision::Bool = true),
    @nospecialize(tol::Real = 1e-6),
    @nospecialize(maxiter::Integer = 10000),
    @nospecialize(drop_singletons::Bool = true),
    @nospecialize(progress_bar::Bool = true),
    @nospecialize(subset::Union{Nothing, AbstractVector} = nothing))

    df = DataFrame(df; copycols = false)
    
    # Prepare the varlists
    nvarlist = Symbol[]
    dvarlist = Symbol[]
    for l in rel_varlist
        push!(dvarlist,l)
        n_l = Symbol("n",l)
        df[!,n_l] = df[!,l]
        df[!,n_l][df[!,control_cohort].==1] .= 0
        push!(nvarlist,n_l)
    end

    cohort_list = unique(df[!,cohort][df[!,control_cohort].==0])
    cohort_list = sort(cohort_list)
    cohort_list = Int.(cohort_list)
    nrel_times = length(nvarlist)
    ncohort = length(cohort_list)

    # Prepare interaction terms for the interacted regression

    cohort_rel_varlist = Symbol[]
    for l in nvarlist
        for yy in cohort_list
            n_l_yy = Symbol("n"*string(l)*"_"*string(yy))
            df[!, :D] .= ifelse.(coalesce.(df[!, cohort] .== yy, false), 1, 0)
            df[!, n_l_yy] = df[!, :D] .* df[!, l]
            push!(cohort_rel_varlist, n_l_yy)
        end
    end

    # bcohort_rel_varlist = Symbol[]
    # for l in rel_varlist
    #     for yy in cohort_list
    #         push!(bcohort_rel_varlist, Symbol(string(l)*"_x_"*string(yy)))
    #     end
    # end

    # Estimate the interacted regression
    formula1 = formula.lhs ~ sum(term.(cohort_rel_varlist)) + formula.rhs

    result =  reg(df,     
    formula1,
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
    subset = subset)

    #Caculate the weights.

    df = df[result.esample .== 1, :]


    ff_w = zeros(0,length(nvarlist))
    nresidlist = Symbol[]
    
    for yy in cohort_list
        cohort_ind = Symbol("cohort_ind")
        resid_yy = Symbol("resid", yy)
        df[!, cohort_ind] .= ifelse.(coalesce.(df[!,cohort] .== yy,false) , 1 , 0)
        formula2 = term(cohort_ind) ~ sum(term.(nvarlist)) + term(0)
        reg1 = reg(df[df[!,control_cohort] .== 0,:], formula2)
        bb = reg1.coef
        ff_w  = vcat(ff_w,bb')
        df[!,resid_yy] = residuals(reg1,df)
        push!(nresidlist,resid_yy)
    end
    
    X = hcat([df[df[!, control_cohort] .== 0, i] for i in nvarlist]...)
    XX = X' * X
    Sxx = XX*(1/size(X)[1])
    Sxxi = FixedEffectModels.invsym!(Sxx)

    e = hcat([df[df[!, control_cohort] .== 0, i] for i in nresidlist]...)

    S_robust = avar(e, X, true) # using the robust estimator

    KSxxi = kron(Matrix{Float64}(I, ncohort, ncohort), Sxxi)

    Sigma_ff = KSxxi * S_robust * KSxxi/ size(X)[1]

    b = result.coef

    V = diag(result.vcov)

    replace!(V, NaN => 0)

    # Convert the delta estimate vector to a matrix where each column is a relative time
    end_ = 0
    evt_bb = zeros(ncohort,0)
    evt_VV = zeros(ncohort,0)

    for i in 1:nrel_times
        start = end_ + 1
        end_ = start + ncohort - 1
        b_i = b[start:end_]
        evt_bb = hcat(evt_bb, b_i)
        V_i = V[start:end_]
        evt_VV = hcat(evt_VV, V_i)
    end

    # Take weighted average for IW estimators
    w = ff_w
    delta = evt_bb
    b_iw = sum(w .* delta, dims=1)
    nc, nr = size(w)

    # Ptwise variance from cohort share estimation and interacted regression
    VV = result.vcov # VCV from the interacted regression
    replace!(VV, NaN => 0)

    VV = VV[1:nr*nc, 1:nr*nc] # in case reports _cons

    wlong = w'.*repeat(ee(1, nr)', 1, nc)

    for i in 2:nrel_times      # loop from 2 to nrel_times
        wlong = [wlong (w' .* repeat(ee(i, nr)', 1, nc))]  # concatenate 
    end

    V_iw = wlong * VV * wlong'


    Vshare = Sigma_ff
    Sigma_l = zeros(0, 0)
    share_idx = 0:nr:(nc-1)*nr

    for i in 1:nrel_times
        for j in 1:i
            Vshare_evt = Vshare[share_idx .+ i, share_idx .+ j]
            V_iw[i,j] += delta[:,i]' * Vshare_evt * delta[:,j]
        end
    end

    coef_iw = vec(b_iw')

    vcov_iw = Symmetric(V_iw)
    vcov_iw_diag = diag(vcov_iw)



    ##############################################################################
    ##
    ## Test Statistics
    ##
    ##############################################################################
    

    # Compute standard error

    # Compute Fstat
    F = result.F
    p = result.p
    dof = result.dof
    dof_residual = result.dof_residual
    # Compute rss, tss, r2, r2 adjusted

    rss = result.rss
    r2 = result.r2
    r2_within = result.r2_within
    tss_total = rss /(1-r2)
    mss = tss_total - rss
    adjr2 = result.adjr2

    # Others
    nclusters = result.nclusters
    esample = result.esample
    fekeys = result.fekeys
    coef_names = [string(x) for x in dvarlist]
    # coef_names = coefnames
    response_name = responsename(result)
    # response_name, coef_names = coefnames(formula_schema)

    formula_origin = formula1
    formula_schema = result.formula_schema
    contrasts = result.contrasts
    nobs = result.nobs
    residuals2 = nothing
    if (save == :residuals) | (save == :all)
        residuals2 = result.residuals2
    end
    iterations = result.iterations
    fe2 = result.fe
    converged = result.converged
    feresult = result

    return EventStudyInteract(coef_iw, vcov_iw, vcov, nclusters, esample, residuals2, fe2, fekeys, coef_names, response_name, formula_origin, formula_schema, contrasts, nobs, dof , dof_residual, rss, tss_total, r2, adjr2, F, p, iterations, converged, r2_within, feresult)
end