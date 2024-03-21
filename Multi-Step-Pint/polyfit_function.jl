function vander!(V::AbstractMatrix, x::AbstractVector, n=length(x))
    m = length(x)
    (m, n) == size(V) || throw(DimensionMismatch())
    for j = 1:m
        @inbounds V[j, 1] = one(x[j])
    end
    for i = 2:n, j = 1:m
        @inbounds V[j, i] = x[j] * V[j, i-1]
    end
    return V
end

function vander(x::AbstractVector, n=length(x))
    return vander!(Array{eltype(x)}(undef, length(x), n), x, n)
end

function polyfit!(x::AbstractVector, y::AbstractVector, n)
    return vander(x, n+1) \ y
end

function polyfit(x::AbstractVector, y::AbstractVector, n)
    return polyfit!(x, y, n)
end