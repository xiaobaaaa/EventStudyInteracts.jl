function avar(e::AbstractMatrix, Z::AbstractMatrix, robust::Bool = true)
    N = size(e, 1)
    L = size(Z, 2)
    p = size(e, 2)
    S = zeros(Float64, L * p, L * p)
    for j in 1:p
        rows_j = (j - 1) * L + 1:j * L
        for k in 1:p
            rows_k = (k - 1) * L + 1:k * L
            if robust
                S[rows_j, rows_k] = (Z' * Diagonal(e[:, j] .* e[:, k]) * Z) / N
            else
                S[rows_j, rows_k] = (Z' * Z) * cov(e[:, j], e[:, k]) / N
            end
        end
    end
    return Symmetric(S)
end

_ginvsym(X::AbstractMatrix) = Matrix(-FixedEffectModels.invsym!(Symmetric(copy(X))))

function _share_vcov_expand(reduced::AbstractMatrix, keep_cols::AbstractVector{Bool}, nresid::Integer)
    full_cols = length(keep_cols)
    keep_idx = findall(keep_cols)
    expanded_idx = Int[]

    for resid_idx in 1:nresid
        append!(expanded_idx, (resid_idx - 1) * full_cols .+ keep_idx)
    end

    full = zeros(Float64, full_cols * nresid, full_cols * nresid)
    full[expanded_idx, expanded_idx] .= Matrix(reduced)
    return Symmetric(full)
end

function _share_vcov_dispatch(sample::DataFrame, emat::Matrix{Float64}, Xmat::Matrix{Float64}, vcov::CovarianceEstimator)
    vcov isa Vcov.SimpleCovariance && return _share_vcov_simple(emat, Xmat)
    vcov isa Vcov.RobustCovariance && return _share_vcov_robust(emat, Xmat)
    vcov isa Vcov.ClusterCovariance && return _share_vcov_cluster(sample, emat, Xmat, vcov)

    throw(ArgumentError("Unsupported covariance estimator $(typeof(vcov)) for the cohort-share step."))
end

function _share_vcov(table, e::AbstractMatrix, X::AbstractMatrix, vcov::CovarianceEstimator)
    sample = DataFrame(table; copycols = false)
    Xmat = Matrix{Float64}(X)
    emat = Matrix{Float64}(e)

    keep = Vcov.completecases(sample, vcov)
    if !all(keep)
        sample = sample[keep, :]
        Xmat = Xmat[keep, :]
        emat = emat[keep, :]
    end

    keep_cols = vec(any(.!iszero.(Xmat), dims = 1))
    if all(keep_cols)
        return _share_vcov_dispatch(sample, emat, Xmat, vcov)
    end

    reduced = _share_vcov_dispatch(sample, emat, Xmat[:, keep_cols], vcov)
    return _share_vcov_expand(reduced, keep_cols, size(emat, 2))
end

function _share_vcov_simple(e::Matrix{Float64}, X::Matrix{Float64})
    return _share_vcov_from_avar(e, X, false)
end

function _share_vcov_robust(e::Matrix{Float64}, X::Matrix{Float64})
    return _share_vcov_from_avar(e, X, true)
end

function _share_vcov_from_avar(e::Matrix{Float64}, X::Matrix{Float64}, robust::Bool)
    N = size(X, 1)
    XX = X' * X
    Sxxi = _ginvsym(XX / N)
    K = kron(Matrix{Float64}(I, size(e, 2), size(e, 2)), Sxxi)
    return Symmetric(K * Matrix(avar(e, X, robust)) * K / N)
end

function _share_vcov_cluster(table::DataFrame, e::Matrix{Float64}, X::Matrix{Float64}, vcov::Vcov.ClusterCovariance)
    XX = X' * X
    XXi = _ginvsym(XX)
    K = kron(Matrix{Float64}(I, size(e, 2), size(e, 2)), XXi)
    S = _clustered_share_meat(table, X, e, vcov.clusternames)
    return Symmetric(K * S * K)
end

function _clustered_share_meat(table::DataFrame, X::Matrix{Float64}, e::Matrix{Float64}, cluster_names::Tuple{Vararg{Symbol}})
    scores = _stacked_scores(X, e)
    S = zeros(Float64, size(scores, 2), size(scores, 2))
    nways = length(cluster_names)

    for mask in 1:(1 << nways) - 1
        names_subset = Tuple(cluster_names[i] for i in 1:nways if ((mask >> (i - 1)) & 1) == 1)
        sign = isodd(count_ones(mask)) ? 1.0 : -1.0
        S .+= sign .* _cluster_crossprod(scores, table, names_subset)
    end

    G = minimum(length(unique(table[!, name])) for name in cluster_names)
    scale = (size(X, 1) - 1) / max(size(X, 1) - size(X, 2), 1)
    G > 1 && (scale *= G / (G - 1))
    return Symmetric(scale .* S)
end

function _stacked_scores(X::Matrix{Float64}, e::Matrix{Float64})
    scores = zeros(Float64, size(X, 1), size(X, 2) * size(e, 2))
    idx = 0
    for k in 1:size(e, 2)
        for j in 1:size(X, 2)
            idx += 1
            scores[:, idx] .= X[:, j] .* e[:, k]
        end
    end
    return scores
end

function _cluster_crossprod(scores::Matrix{Float64}, table::DataFrame, cluster_names::Tuple)
    grouped_scores = Dict{Any, Vector{Float64}}()
    columns = map(name -> table[!, name], cluster_names)

    for i in 1:size(scores, 1)
        key = length(columns) == 1 ? columns[1][i] : tuple((column[i] for column in columns)...)
        bucket = get!(grouped_scores, key) do
            zeros(Float64, size(scores, 2))
        end
        bucket .+= @view scores[i, :]
    end

    S = zeros(Float64, size(scores, 2), size(scores, 2))
    for bucket in values(grouped_scores)
        S .+= bucket * bucket'
    end
    return S
end
