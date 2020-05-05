using Statistics
using Distributions

struct MarcenkoPastur{T<:UnivariateDistribution}
    rho::Real
end

# Indicator function 1(x)
# Parameters
#   x::Real     Real number of which we are evaluating 1(x) at
#   a::Real     Lower bound of our set
#   b::Real     Upper bound of our set
function indicator(x::Real, a::Real, b::Real)
    if a <= x <= b
        return 1
    else
        return 0
end

alpha(rho) = (1-√rho)^2
beta(rho) = (1+√rho)^2
function pdf(d::UnivariateDistribution, x::Real)
    √((beta(rho)-x)(x-alpha(rho)))//(2πx)*indicator(x, alpha(rho),beta(rho))
end