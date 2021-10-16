struct Risk
    λ::Float64
    α::Float64
    offset::Union{Dict{Leaf,Float64},Nothing}
    bound::Union{Float64,Nothing}
    penalty::Union{Float64,Nothing}
end

"""
	RiskNeutral()

Create a risk-neutral risk measure.
"""
function RiskNeutral()
    return Risk(1.0, 1.0)
end

"""
	Risk(λ::Float64,
         α::Float64;
         offset::Union{Dict{Leaf,Float64},Nothing}=nothing,
         bound::Union{Float64,Nothing}=nothing,
         penalty::Union{Float64,Nothing} = nothing)

Define the CVaR risk measure to be applied to the accumulated profits at the leaf nodes.

### Required Arguments
`λ` is weighting applied for the risk measure (max sum of weightings should be 1.0),
if sum of weightings is less than 1.0, expected value will make up the rest.

`α` is the probability in the tail of the distribution

### Optional Arguments
`offset` applies a negative offset to each leaf node. This can be used to reorder
the outcomes prior to applying the risk measure.

`bound` if used, this will create a constraint on CVaR(α) with this as the upper bound.

`penalty` if a constraint on CVaR is applied, then the marginal cost of violating the
constraint is `penalty`; if set to `nothing` then no constraint violation is allowed.)
"""
function Risk(
    λ::Float64,
    α::Float64;
    offset::Union{Dict{Leaf,Float64},Nothing} = nothing,
    bound::Union{Float64,Nothing} = nothing,
    penalty::Union{Float64,Nothing} = nothing,
)
    if α <= 0.0 || α > 1.0
        error("α must be >0 and <=1")
    elseif λ < 0.0 || λ > 1.0
        error("λ must be >=0 and <=1")
    end
    return Risk(λ, α, offset, bound, penalty)
end

"""
	Risk(α::Float64;
         offset::Union{Dict{Leaf,Float64},Nothing}=nothing,
         bound::Union{Float64,Nothing}=nothing,
         penalty::Union{Float64,Nothing} = nothing)

Define the CVaR risk constraint to be applied to the accumulated profits at the leaf nodes.

### Required Arguments
`λ` is weighting applied for the risk measure (max sum of weightings should be 1.0),
if sum of weightings is less than 1.0, expected value will make up the rest.

`α` is the probability in the tail of the distribution

### Optional Arguments
`offset` applies a negative offset to each leaf node. This can be used to reorder
the outcomes prior to applying the risk measure.

`bound` if used, this will create a constraint on CVaR(α) with this as the upper bound.

`penalty` if a constraint on CVaR is applied, then the marginal cost of violating the
constraint is `penalty`; if set to `nothing` then no constraint violation is allowed.
"""
function Risk(
    α::Float64;
    offset::Union{Dict{Leaf,Float64},Nothing} = nothing,
    bound::Union{Float64,Nothing} = nothing,
    penalty::Union{Float64,Nothing} = nothing,
)
    if α <= 0.0 || α > 1.0
        error("α must be >0 and <=1")
    end
    if bound == nothing
        error("Risk constraint cannot be created without setting a bound")
    end
    return Risk(0.0, α, offset, bound, penalty)
end
