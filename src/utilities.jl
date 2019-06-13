# helper function which changes the objective function coeffecient of a particular variable
# this really should be a jump standard
function changeobjcoef!(var::JuMP.Variable,coef::Float64)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        m.obj += var*coef
    else
        m.obj.aff.coeffs[pos] = coef
    end
end

# helper function which fetches the current objective value coef for a given variable
function getcoef(var::JuMP.Variable)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        return 0.0
    else
        return m.obj.aff.coeffs[pos]
    end
end

# get values of a variable from a node
function getvalueDW(jmodel::JuDGEModel, node::Node, var::Symbol)
    if jmodel.isbuilt
        return getvalue(jmodel.subprob[node][var])
    else
        println("Decomposition model not built.")
    end
end

function getvalueDW(jmodel::JuDGEModel, indices::Array{Int64,1}, var::Symbol)
    if jmodel.isbuilt
        return getvalue(jmodel.subprob[getnode(jmodel.tree,indices)][var])
    else
        println("Decomposition model not built.")
    end
end
