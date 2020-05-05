using Statistics
using Distributions

struct MarcenkoPastur{T<:UnivariateDistribution}
    rho::Real
end

alpha(rho) = (1-√rho)^2
beta(rho) = (1+√rho)^2
pdf(d::UnivariateDistribution, x::Real) = √((beta(rho)-x)(x-alpha(rho)))//(2πx)