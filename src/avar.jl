
function avar(e::Matrix, Z::Matrix, robust::Bool=true)
    N = size(e, 1)
    L = size(Z, 2)
    p = size(e, 2)
    S = zeros(L*p, L*p)
    for j in 1:p
        for k in 1:p
            if robust
                S[(j-1)*L+1:j*L,(k-1)*L+1:k*L] = (Z' * Diagonal(e[:,j] .* e[:,k]) * Z) / N
            else
                S[(j-1)*L+1:j*L,(k-1)*L+1:k*L] = (Z' * Z) * cov(e[:,j], e[:,k]) / N
            end
        end
    end
    return S
end