function _nonzero_relmask(df::DataFrame, nvarlist::Vector{Symbol})
    isempty(nvarlist) && return Bool[]
    X = Matrix{Float64}(df[:, nvarlist])
    return vec(any(.!iszero.(X), dims = 1))
end

function _estimable_relmask(vcov_mat::AbstractMatrix, nrel_times::Integer, ncohort::Integer)
    nrel_times == 0 && return Bool[]
    interaction_dim = nrel_times * ncohort
    diag_v = diag(Matrix(vcov_mat)[1:interaction_dim, 1:interaction_dim])
    keep_rel = falses(nrel_times)
    for i in 1:nrel_times
        block = ((i - 1) * ncohort + 1):(i * ncohort)
        keep_rel[i] = any(.!isnan.(diag_v[block]))
    end
    return keep_rel
end

function _zero_omitted_rowscols(vcov_mat::AbstractMatrix)
    VV = Matrix{Float64}(vcov_mat)
    omitted = isnan.(diag(VV))
    kept = .!omitted

    if any(kept) && any(isnan.(VV[kept, kept]))
        throw(ArgumentError("Non-omitted entries of the interacted-regression covariance matrix contain NaN values."))
    end

    if any(omitted)
        VV[omitted, :] .= 0.0
        VV[:, omitted] .= 0.0
    end

    return VV
end

function _reinsert_rel_outputs(coef_reduced::AbstractVector, vcov_reduced::AbstractMatrix, keep_rel::AbstractVector{Bool}, full_size::Integer)
    coef_full = zeros(Float64, full_size)
    vcov_full = fill(NaN, full_size, full_size)

    keep_idx = findall(keep_rel)
    coef_full[keep_idx] .= coef_reduced
    vcov_full[keep_idx, keep_idx] .= Matrix(vcov_reduced)

    return coef_full, Symmetric(vcov_full)
end

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
    @nospecialize(dof_add::Integer = 0),
    @nospecialize(subset::Union{Nothing, AbstractVector} = nothing))
    StatsAPI.fit(EventStudyInteract, formula, df, rel_varlist, control_cohort, cohort, vcov; contrasts = contrasts, weights = weights, save = save, method = method, nthreads = nthreads, double_precision = double_precision, tol = tol, maxiter = maxiter, drop_singletons = drop_singletons, progress_bar = progress_bar, dof_add = dof_add, subset = subset)
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
    @nospecialize(dof_add::Integer = 0),
    @nospecialize(subset::Union{Nothing, AbstractVector} = nothing))

    df = DataFrame(df; copycols = false)
    if dof_add != 0
        @warn "The dof_add keyword is no longer forwarded to FixedEffectModels.reg on version 1.13 and is currently ignored."
    end    
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
    full_dvarlist = copy(dvarlist)

    prekeep_rel = _nonzero_relmask(df, nvarlist)
    if !all(prekeep_rel)
        nvarlist = nvarlist[prekeep_rel]
        dvarlist = dvarlist[prekeep_rel]
    end
    isempty(nvarlist) && error("No estimable relative-time indicators remain after removing unsupported columns.")

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

    postkeep_rel = _estimable_relmask(result.vcov, nrel_times, ncohort)
    if !all(postkeep_rel)
        nvarlist = nvarlist[postkeep_rel]
        dvarlist = dvarlist[postkeep_rel]
    end
    isempty(nvarlist) && error("No estimable relative-time indicators remain after fitting the interacted regression.")

    nrel_times = length(nvarlist)
    rel_keep_mask = falses(length(rel_varlist))
    rel_keep_mask[findall(prekeep_rel)] .= postkeep_rel

    #Caculate the weights.

    df = df[result.esample .== 1, :]

    if ncohort == 1
        ff_w = ones(1, nrel_times)
        Sigma_ff = zeros(nrel_times, nrel_times)
    else
        ff_w = zeros(0, length(nvarlist))
        nresidlist = Symbol[]
        stage2df = df[df[!, control_cohort] .== 0, :]

        for yy in cohort_list
            cohort_ind = Symbol("cohort_ind")
            resid_yy = Symbol("resid", yy)
            stage2df[!, cohort_ind] .= ifelse.(coalesce.(stage2df[!,cohort] .== yy,false) , 1 , 0)
            formula2 = term(cohort_ind) ~ sum(term.(nvarlist)) + term(0)
            reg1 = reg(stage2df, formula2)
            bb = reg1.coef
            ff_w  = vcat(ff_w,bb')
            stage2df[!,resid_yy] = residuals(reg1, stage2df)
            push!(nresidlist,resid_yy)
        end

        X = Matrix{Float64}(stage2df[:, nvarlist])
        e = Matrix{Float64}(stage2df[:, nresidlist])
        Sigma_ff = _share_vcov(stage2df, e, X, vcov)
    end

    interaction_keep = repeat(postkeep_rel, inner = ncohort)
    interaction_idx = findall(interaction_keep)
    b = result.coef[interaction_idx]

    # Convert the interacted coefficients to a matrix where each column is a relative time
    evt_bb = reshape(b, ncohort, nrel_times)

    # Take weighted average for IW estimators
    w = ff_w
    delta = evt_bb
    b_iw = sum(w .* delta, dims=1)
    nc, nr = size(w)

    # Ptwise variance from cohort share estimation and interacted regression
    VV_full = _zero_omitted_rowscols(result.vcov) # VCV from the interacted regression
    VV = VV_full[interaction_idx, interaction_idx]

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

    coef_iw_reduced = vec(b_iw')
    vcov_iw_reduced = Symmetric(V_iw)
    coef_iw, vcov_iw = _reinsert_rel_outputs(coef_iw_reduced, vcov_iw_reduced, rel_keep_mask, length(rel_varlist))



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
    r2 = StatsAPI.r2(result)
    r2_within = result.r2_within
    tss_total = rss /(1-r2)
    mss = tss_total - rss
    adjr2 = StatsAPI.adjr2(result)

    # Others
    nclusters = result.nclusters
    esample = result.esample
    fekeys = result.fekeys
    coef_names = [string(x) for x in full_dvarlist]
    # coef_names = coefnames
    response_name = responsename(result)
    # response_name, coef_names = coefnames(formula_schema)

    formula_origin = formula1
    formula_schema = result.formula_schema
    contrasts = result.contrasts
    nobs = result.nobs
    residuals2 = nothing
    if (save == :residuals) | (save == :all)
        residuals2 = result.residuals
    end
    iterations = result.iterations
    fe2 = result.fe
    converged = result.converged
    feresult = result

    return EventStudyInteract(coef_iw, vcov_iw, vcov, nclusters, esample, residuals2, fe2, fekeys, coef_names, response_name, formula_origin, formula_schema, contrasts, nobs, dof , dof_residual, rss, tss_total, r2, adjr2, F, p, iterations, converged, r2_within, feresult)
end






