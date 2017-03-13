## Copyright (c) 2015, Júlio Hoffimann Mendes <juliohm@stanford.edu>
##
## Permission to use, copy, modify, and/or distribute this software for any
## purpose with or without fee is hereby granted, provided that the above
## copyright notice and this permission notice appear in all copies.
##
## THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
## WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
## ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
## WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
## ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
## OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

abstract AbstractEstimator
@doc doc"""
    estimate(est, xₒ)

  Evaluate estimator `est` at location `xₒ`
""" ->
estimate(::AbstractEstimator, xₒ::AbstractVector) = nothing

@doc doc"""
  Simple Kriging

  *INPUTS*:

    * X  ∈ ℜ^(mxn) - matrix of data locations
    * z  ∈ ℜⁿ      - vector of observations for X
    * cov          - covariance model
    * μ  ∈ ℜ       - mean of z
  """ ->
type SimpleKriging{T<:Real} <: AbstractEstimator
  # input fields
  X::AbstractMatrix{T}
  z::AbstractVector{T}
  cov::CovarianceModel
  μ::T

  # state fields
  C::AbstractMatrix{T}

  function SimpleKriging(X, z, cov, μ)
    @assert size(X, 2) == length(z) "incorrect data configuration"
    C = pairwise(cov, X)
    new(X, z, cov, μ, C)
  end
end
SimpleKriging(X, z, cov, μ) = SimpleKriging{eltype(z)}(X, z, cov, μ)


@doc doc"""
  Ordinary Kriging

  *INPUTS*:

    * X  ∈ ℜ^(mxn) - matrix of data locations
    * z  ∈ ℜⁿ      - vector of observations for X
    * cov          - covariance model
  """ ->
type OrdinaryKriging{T<:Real} <: AbstractEstimator
  # input fields
  X::AbstractMatrix{T}
  z::AbstractVector{T}
  cov::CovarianceModel

<<<<<<< HEAD
        # estimate and variance
        μ + y⋅λ, cov(0) - c⋅λ, C
    else                        # Ordinary Kriging
        C = [C ones(n); ones(n)' 0]
        c = [c; 1]
        λ = C \ c

        # estimate and variance
        z⋅λ[1:n], cov(0) - c⋅λ, C
    end
=======
  # state fields
  C::AbstractMatrix{T}

  function OrdinaryKriging(X, z, cov)
    @assert size(X, 2) == length(z) "incorrect data configuration"
    C = pairwise(cov, X)
    new(X, z, cov, C)
  end
>>>>>>> upstream/master
end
OrdinaryKriging(X, z, cov) = OrdinaryKriging{eltype(z)}(X, z, cov)


@doc doc"""
  Universal Kriging (a.k.a. Kriging with drift)

  *INPUTS*:

    * X  ∈ ℜ^(mxn) - matrix of data locations
    * z  ∈ ℜⁿ      - vector of observations for X
    * cov          - covariance model
    * degree       - polynomial degree for the mean

  Ordinary Kriging is recovered for 0th degree polynomial.
  """ ->
type UniversalKriging{T<:Real} <: AbstractEstimator
  # input fields
  X:: AbstractMatrix{T}
  z::AbstractVector{T}
  cov::CovarianceModel
  degree::Integer

  # state fields
  Γ::AbstractMatrix{T}
  F::AbstractMatrix{T}
  exponents::AbstractMatrix{Float64}

  function UniversalKriging(X, z, cov, degree)
    @assert size(X, 2) == length(z) "incorrect data configuration"
    @assert degree ≥ 0

<<<<<<< HEAD
    # X = round(X, 6)
    # z = round(z, 6)
    # x₀ = round(x₀, 6)
    println("X[:,1]: $x₀")

    γ(h) = sill(cov) - cov(h)

    dim = length(x₀)
=======
    dim = size(X, 1)
    nobs = length(z)
>>>>>>> upstream/master

    γ(h) = cov(0) - cov(h)
    Γ = pairwise(γ, X)
<<<<<<< HEAD

    gv = Array{typeof(X[1,1])}(n)
    gvn = Array{typeof(X[1,1])}(n)
    for j=1:n
      gvn[j] = norm(X[:,j]-x₀)
      gv[j] = γ(gvn[j])
    end
    #g = round(gv, 39)
    g = gv
    #g = Float64[γ(norm(X[:,j]-x₀)) for j=1:n]
    println("typeof g $(typeof(g))")
    #g = Float64[gv]
=======
>>>>>>> upstream/master

    # multinomial expansion
    exponents = zeros(0, dim)
    for d=0:degree
      exponents = [exponents; multinom_exp(dim, d, sortdir="descend")]
    end
    exponents = exponents'

    nterms = size(exponents, 2)
    println(nterms)

    F = Float64[prod(X[:,i].^exponents[:,j]) for i=1:nobs, j=1:nterms]

    new(X, z, cov, degree, Γ, F, exponents)
  end
end
UniversalKriging(X, z, cov, degree) = UniversalKriging{eltype(z)}(X, z, cov, degree)

##########################
### ESTIMATION METHODS ###
##########################

function estimate{T<:Real}(estimator::SimpleKriging{T}, xₒ::AbstractVector{T})
  X = estimator.X; z = estimator.z
  cov = estimator.cov; μ = estimator.μ
  C = estimator.C

  @assert size(X, 1) == length(xₒ) "incorrect location dimension"

  # evaluate covariance at location
  c = Float64[cov(norm(X[:,j]-xₒ)) for j=1:size(X,2)]

  # solve linear system
  y = z - μ
  λ = C \ c

  # return estimate and variance
  μ + y⋅λ, cov(0) - c⋅λ
end


function estimate{T<:Real}(estimator::OrdinaryKriging{T}, xₒ::AbstractVector{T})
  X = estimator.X; z = estimator.z; cov = estimator.cov
  C = estimator.C

  @assert size(X, 1) == length(xₒ) "incorrect location dimension"

  nobs = length(z)

  # evaluate covariance at location
  c = Float64[cov(norm(X[:,j]-xₒ)) for j=1:nobs]

  # solve linear system
  A = [C ones(nobs); ones(nobs)' 0]
  b = [c; 1]
  λ = A \ b

  # return estimate and variance
  z⋅λ[1:nobs], cov(0) - b⋅λ
end


function estimate{T<:Real}(estimator::UniversalKriging{T}, xₒ::AbstractVector{T})
  X = estimator.X; z = estimator.z
  cov = estimator.cov; degree = estimator.degree
  Γ = estimator.Γ; F = estimator.F; exponents = estimator.exponents

  @assert size(X, 1) == length(xₒ) "incorrect location dimension"

  nobs = length(z)

  # evaluate variogram at location
  γ(h) = cov(0) - cov(h)
  g = Float64[γ(norm(X[:,j]-xₒ)) for j=1:nobs]

<<<<<<< HEAD
    #
    # F = Array{Float64,2}(n,nterms)
    # for i=1:n, j=1:nterms
    #   F[i,j] = prod(X[:,i].^exponents[:,j])
    # end
    # println("size of F: $(size(F))")
    #
    # f = Array{Float64,1}(nterms)
    # for j=1:nterms
    #   f[j] = prod(x₀.^exponents[:,j])
    # end

    A = [Γ F; F' zeros(nterms,nterms)]
    a = [g; f]

    λ = A \ a
    println("Condition number of A: $(cond(A))")
    println("Condition number of F: $(cond(F))")
    println("Condition number of Γ: $(cond(Γ))")

    # estimate and variance
    z⋅λ[1:n], a⋅λ, Γ
=======
  # evaluate multinomial at location
  nterms = size(exponents, 2)
  f = Float64[prod(xₒ.^exponents[:,j]) for j=1:nterms]

  # solve linear system
  A = [Γ F; F' zeros(nterms,nterms)]
  b = [g; f]
  λ = A \ b

  # estimate and variance
  z⋅λ[1:nobs], b⋅λ
>>>>>>> upstream/master
end
