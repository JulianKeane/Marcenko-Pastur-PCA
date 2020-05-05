using Statistics
using Distributions
import Distributions: @check_args
import Distributions: @distr_support

struct MarcenkoPastur{T<:Real} <: ContinuousUnivariateDistribution
    ρ::T
    MarcenkoPastur{T}(ρ::T) where {T<:Real} = new{T}(ρ)
end

function MarcenkoPastur(ρ::T, check_args=true) where {T <: Real}
    check_args && @check_args(MarcenkoPastur, ρ > zero(ρ))
    return MarcenkoPastur{T}(ρ)
end

MarcenkoPastur(ρ::Integer) = MarcenkoPastur(float(ρ))

@distr_support MarcenkoPastur (1-√d.ρ)^2 (1+√d.ρ)^2

#### Conversions
function convert(::Type{MarcenkoPastur{T}}, ρ::Real) where T<:Real
    MarcenkoPastur(T(ρ))
end
function convert(::Type{MarcenkoPastur{T}},d::MarcenkoPastur{S}) where {T <: Real, S<: Real}
    MarcenkoPastur(T(d.ρ),check_args=false)
end

### Parameters
leftbound(d::MarcenkoPastur) = (1-√d.ρ)^2
rightbound(d::MarcenkoPastur) = (1+√d.ρ)^2
params(d::MarcenkoPastur) = (d.ρ)

### Statistics
mean(d::MarcenkoPastur) = d.ρ
# TODO
# median(d::MarcenkoPastur)
# mode(d::MarcenkoPastur) = 2*leftbound(d)*rightbound(d)/(leftbound(d)+rightbound(d))
# var(d::MarcenkoPastur)

# TODO : make a separate file for functions such as 1(x)
# Indicator function 1(x)
# Parameters
#   x::Real     Real number of which we are evaluating 1(x) at
#   a::Real     Lower bound of our set
#   b::Real     Upper bound of our set
#function indicator(x::Real, a::Real, b::Real)
#    if a <= x <= b
#        return 1
#    else
#        return 0
#    end
#end

α(ρ) = (1-√ρ)^2
β(ρ) = (1+√ρ)^2
function pdf(d::MarcenkoPastur, x::Real)
    if ρ < 1
        return pdf(d, d.ρ*x)
    √((β(d.ρ)-x)*(x-α(d.ρ)))/(2π*x)*indicator(x, α(d.ρ), β(d.ρ))
end