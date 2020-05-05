"""
    MarcenkoPastur()
The Marcenko Pastur Distrbitution (sometimes called the Marchenko-Pastur 
Distribution or Law) has a probability distribution

'''math
f(x) = \frac{√{β(ϱ) - x)(x - α(ϱ))}}{2πx} 1_{[α(ϱ),β(ϱ)]}(x)
'''
where 1_{[a,b]}(x) is the indicator function that returns 1 if x exists in
[a,b] and returns 0 otherwise

External Links

* https://en.wikipedia.org/wiki/Marchenko%E2%80%93Pastur_distribution

"""

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

α(ρ) = (1-√ρ)^2
β(ρ) = (1+√ρ)^2

### Parameters
leftbound(d::MarcenkoPastur) = α(d.ρ)
rightbound(d::MarcenkoPastur) = β(d.ρ)
params(d::MarcenkoPastur) = (d.ρ)

### Statistics
mean(d::MarcenkoPastur) = d.ρ
# TODO
# median(d::MarcenkoPastur)
# mode(d::MarcenkoPastur) = 2*α(d.ρ)*β(d.ρ)/(α(d.ρ)+β(d.ρ))
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

function pdf(d::MarcenkoPastur, x::Real)
    if ρ < 1
        return pdf(d, d.ρ*x)
    √((β(d.ρ)-x)*(x-α(d.ρ)))/(2π*x)*indicator(x, α(d.ρ), β(d.ρ))
end