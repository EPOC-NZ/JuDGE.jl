function copy_variable!(toModel, variable)
    map(variable) do x
        return JuMP.add_variable(
            toModel,
            JuMP.build_variable(error, (get_info(x))),
        )
    end
end

function copy_variable!(toModel, variable, f)
    map(variable) do x
        return JuMP.add_variable(
            toModel,
            JuMP.build_variable(error, f(get_info(x))),
        )
    end
end

function copy_variable!(toModel, variable::JuMP.VariableRef)
    return JuMP.add_variable(
        toModel,
        JuMP.build_variable(error, get_info(variable)),
    )
end

function copy_variable!(toModel, variable::JuMP.VariableRef, f)
    return JuMP.add_variable(
        toModel,
        JuMP.build_variable(error, f(get_info(variable))),
        #        "exp",
    )
end

# constuct variable info object for a single variable
function get_info(x::VariableRef)
    has_lb_local = false
    lb_local = NaN
    has_ub_local = false
    ub_local = NaN
    is_fixed_local = false
    fixed_value_local = NaN
    has_start_local = false
    start_value_local = NaN
    is_binary_local = false
    is_integer_local = false

    if has_lower_bound(x)
        has_lb_local = true
        lb_local = lower_bound(x)
    end

    if has_upper_bound(x)
        has_ub_local = true
        ub_local = upper_bound(x)
    end

    if is_fixed(x)
        is_fixed_local = true
        fixed_value_local = fix_value(x)
    end

    if start_value(x) != nothing
        has_start_local = true
        start_value_local = start_value(x)
    end

    if is_binary(x)
        is_binary_local = true
    end

    if is_integer(x)
        is_integer_local = true
    end

    return VariableInfo(
        has_lb_local,
        lb_local,
        has_ub_local,
        ub_local,
        is_fixed_local,
        fixed_value_local,
        has_start_local,
        start_value_local,
        is_binary_local,
        is_integer_local,
    )
end

function relaxbinary(x::VariableInfo)
    return VariableInfo(
        true,
        0.0,
        true,
        1.0,
        x.has_fix,
        x.fixed_value,
        x.has_start,
        x.start,
        false,
        false,
    )
end

function relaxinteger(x::VariableInfo)
    return VariableInfo(
        true,
        x.lower_bound,
        true,
        x.upper_bound,
        x.has_fix,
        x.fixed_value,
        x.has_start,
        x.start,
        false,
        false,
    )
end

function UnitIntervalInformation(; UB::Float64 = 1.0)
    if UB == 1.0
        return VariableInfo(
            true,
            0.0,
            false,
            NaN,
            false,
            NaN,
            false,
            NaN,
            false,
            false,
        )
    else
        return VariableInfo(
            true,
            0.0,
            true,
            UB,
            false,
            NaN,
            false,
            NaN,
            false,
            false,
        )
    end
end

function coef(aff, x::JuMP.VariableRef)
    if x in keys(aff)
        aff[x]
    else
        0.0
    end
end

function densekey_to_tuple(key::Any)
    if typeof(key) <: JuMP.Containers.DenseAxisArrayKey
        return length(key.I) == 1 ? key.I[1] : key.I
    else
        return key
    end
end

function unpack_expansions(a::Dict{AbstractTree,Dict{Symbol,Any}})
    assign = Vector{Expr}()
    found = Dict{Symbol,Bool}()
    for (i, j) in a
        for (k, l) in j
            if k in keys(found)
                push!(assign, :($k[$i] = $l))
            else
                if isdefined(Main, k)
                    m = Symbol("→" * string(k))
                    push!(assign, :($m = $k))
                end
                push!(assign, :($k = Dict{AbstractTree,Any}()))
                push!(assign, :($k[$i] = $l))
                found[k] = true
            end
        end
    end
    return assign
end

function clear_expansions(a::Dict{AbstractTree,Dict{Symbol,Any}})
    assign = Vector{Expr}()
    for (i, j) in a
        for (k, l) in j
            m = Symbol("→" * string(k))
            if isdefined(Main, m)
                push!(assign, :($k = $m))
                push!(assign, :($m = nothing))
            else
                push!(assign, :($k = nothing))
            end
        end
        break
    end
    return assign
end

function compute_objval(
    scenarios::Dict{Leaf,Float64},
    probabilities::Dict{AbstractTree,Float64},
    risk::Union{Risk,Vector{Risk}},
)
    scenario_objs = Vector{Tuple{Float64,Float64,Leaf}}()
    EV_weight = 1.0
    EV = 0.0
    for (leaf, val) in scenarios
        push!(scenario_objs, (probabilities[leaf], val, leaf))
        EV += probabilities[leaf] * val
    end

    obj = 0.0

    if typeof(risk) == Risk
        if risk.α == 1.0 || risk.λ == 0.0
            risk = []
        else
            risk = [risk]
        end
    end
    for i in 1:length(risk)
        EV_weight -= risk[i].λ
        so = scenario_objs
        if risk[i].offset != nothing
            for j in 1:length(so)
                so[j] =
                    (so[j][1], so[j][2] - risk[i].offset[so[j][3]], so[j][3])
            end
        end
        sort!(so, by = i -> i[2], rev = true)
        beta = risk[i].α
        for scen in so
            if scen[1] > beta
                obj += scen[2] * risk[i].λ * beta / risk[i].α
                beta = 0
            else
                obj += scen[2] * risk[i].λ * scen[1] / risk[i].α
                beta -= scen[1]
            end
        end
    end

    return obj + EV_weight * EV
end

"""
    solution_to_dictionary(jmodel::JuDGEModel; prefix::String = "")

Create a nested dictionary with the solution values for each node copied across to
a standardised structure.

### Required Arguments
`jmodel` is a solved JuDGEModel

### Optional Arguments
`prefix` is a string that will be prepended to each of the variable names (e.g. if comparing to versions of the same model)
"""
function solution_to_dictionary(jmodel::JuDGEModel; prefix::String = "")
    function helper(node::AbstractTree, solution::T where {T<:Dict})
        solution[node] = Dict{Symbol,Any}()
        var_groups = JuMP.object_dictionary(jmodel.sub_problems[node])
        for (var_group, vars) in var_groups
            sym = Symbol(prefix * string(var_group))
            if typeof(vars) == VariableRef
                solution[node][sym] = JuMP.value(vars)
            elseif typeof(vars) <: AbstractArray
                skip = false
                for v in vars
                    if typeof(v) != VariableRef
                        skip = true
                    end
                    break
                end
                if skip
                    continue
                end
                if sym ∉ keys(solution[node])
                    solution[node][sym] = Dict{String,Float64}()
                end
                vals = JuMP.value.(vars)
                if typeof(vars) <: JuMP.Containers.SparseAxisArray
                    vals = vals.data
                end

                for key in keys(vals)
                    temp = ""
                    if typeof(vals) <: Array
                        strkey = string(key)
                        strkey = replace(strkey, "CartesianIndex(" => "")
                        strkey = replace(strkey, ")" => "")
                        strkey = replace(strkey, ", " => ",")
                        temp = strkey
                    elseif typeof(vals) <: Dict
                        strkey = string(key)
                        strkey = replace(strkey, ")" => "")
                        strkey = replace(strkey, "(" => "")
                        strkey = replace(strkey, ", " => ",")
                        temp = strkey
                    else
                        strkey = string(key.I)
                        strkey = replace(strkey, ")" => "")
                        strkey = replace(strkey, "(" => "")
                        strkey = replace(strkey, "\"" => "")
                        strkey = replace(strkey, ", " => ",")
                        temp = strkey
                    end
                    solution[node][sym][temp] = vals[key]
                end
            end
        end

        for (x, var) in jmodel.master_problem.ext[:expansions][node]
            sym = Symbol(prefix * string(x) * "_master")
            if typeof(var) <: AbstractArray
                if sym ∉ keys(solution[node])
                    solution[node][sym] = Dict{String,Float64}()
                end
                val = JuMP.value.(var)
                if typeof(var) <: JuMP.Containers.SparseAxisArray
                    val = val.data
                end
                for key in keys(val)
                    temp = ""
                    if typeof(val) <: Array
                        strkey = string(key)
                        strkey = replace(strkey, "CartesianIndex(" => "")
                        strkey = replace(strkey, ")" => "")
                        strkey = replace(strkey, ", " => ",")
                        temp *= strkey
                    elseif typeof(val) <: Dict
                        strkey = string(key)
                        strkey = replace(strkey, ")" => "")
                        strkey = replace(strkey, "(" => "")
                        strkey = replace(strkey, ", " => ",")
                        temp *= strkey
                    else
                        for i in 1:length(val.axes)-1
                            temp *= string(key[i]) * ","
                        end
                        temp *= string(key[length(val.axes)])
                    end

                    solution[node][sym][temp] = val[key]
                end
            else
                solution[node][sym] = JuMP.value(var)
            end
        end

        if typeof(node) == Tree
            for child in node.children
                helper(child, solution)
            end
        else
            solution[node][Symbol(prefix * "scenario_obj")] =
                JuMP.value(jmodel.master_problem.ext[:scenprofit_var][node])
        end
        return solution
    end

    if termination_status(jmodel.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel.master_problem) != MOI.TIME_LIMIT &&
       termination_status(jmodel.master_problem) != MOI.LOCALLY_SOLVED
        error("You need to first solve the decomposed model.")
    end

    if prefix != ""
        prefix *= "_"
    end

    solution = Dict{AbstractTree,Dict{Symbol,Any}}()
    return helper(jmodel.tree, solution)
end

"""
    solution_to_dictionary(deteq::DetEqModel; prefix::String = "")

Create a nested dictionary with the solution values for each node copied across to
a standardised structure.

### Required Arguments
`deteq` is a solved DetEqModel

### Optional Arguments
`prefix` is a string that will be prepended to each of the variable names (e.g. if comparing to versions of the same model)
"""
function solution_to_dictionary(deteq::DetEqModel; prefix::String = "")
    function helper(node::AbstractTree, solution::T where {T<:Dict})
        solution[node] = Dict{Symbol,Any}()

        for (x, var) in deteq.problem.ext[:vars][node]
            ss = split(string(x), '[')
            if length(ss) == 1
                solution[node][Symbol(prefix * ss[1])] = JuMP.value(var)
            else
                if Symbol(prefix * ss[1]) ∉ keys(solution[node])
                    solution[node][Symbol(prefix * ss[1])] =
                        Dict{String,Float64}()
                end
                solution[node][Symbol(prefix * ss[1])][ss[2][1:end-1]] =
                    JuMP.value(var)
            end
        end

        for (x, var) in deteq.problem.ext[:master_vars][node]
            if typeof(var) == VariableRef
                ss = split(string(x), '[')
                ss[1] *= "_master"
                if length(ss) == 1
                    solution[node][Symbol(prefix * ss[1])] = JuMP.value(var)
                else
                    # if Symbol(prefix * ss[1]) ∉ keys(solution[node])
                    #     solution[node][Symbol(prefix * ss[1])] =
                    #         Dict{String,Float64}()
                    # end
                    # solution[node][Symbol(prefix * ss[1])][ss[2][1:end-1]] = JuMP.value(var)
                end
            elseif typeof(var) <: Dict
                for i in eachindex(var)
                    name = deteq.problem.ext[:master_names][node][x][i]
                    ss = split(string(name), '[')
                    ss[1] *= "_master"
                    if Symbol(prefix * ss[1]) ∉ keys(solution[node])
                        solution[node][Symbol(prefix * ss[1])] =
                            Dict{String,Float64}()
                    end
                    solution[node][Symbol(prefix * ss[1])][ss[2][1:end-1]] =
                        JuMP.value(var[i])
                end
            end
        end

        if typeof(node) == Tree
            for child in node.children
                helper(child, solution)
            end
        else
            solution[node][Symbol(prefix * "scenario_obj")] =
                JuMP.value(deteq.problem.ext[:scenario_obj][node])
        end
        return solution
    end

    if termination_status(deteq.problem) != MOI.OPTIMAL &&
       termination_status(deteq.problem) != MOI.TIME_LIMIT &&
       termination_status(deteq.problem) != MOI.INTERRUPTED
        error("You need to first solve the deterministic equivalent model.")
    end

    if prefix != ""
        prefix *= "_"
    end

    solution = Dict{AbstractTree,Dict{Symbol,Any}}()
    return helper(deteq.tree, solution)
end

"""
    set_starting_solution!(deteq::DetEqModel, jmodel::JuDGEModel)

Use the best solution from a JuDGEModel to warm start the deterministic equivalent.

### Required Arguments
`deteq` is an unsolved DetEqModel
`jmodel` is a solved JuDGEModel
"""
function set_starting_solution!(deteq::DetEqModel, jmodel::JuDGEModel)
    for node in collect(jmodel.tree)
        for (name, exps) in jmodel.master_problem.ext[:expansions][node]
            for index in
                keys(jmodel.master_problem.ext[:expansions][node][name])
                key = JuDGE.densekey_to_tuple(index)
                set_start_value(
                    deteq.problem.ext[:master_vars][node][name][key],
                    JuMP.value(
                        jmodel.master_problem.ext[:expansions][node][name][index],
                    ),
                )
            end
        end
    end

    for (node, sp) in jmodel.sub_problems
        temp = Dict{String,VariableRef}()
        for variable in all_variables(sp)
            temp[string(variable)] = variable
        end
        for variable in keys(deteq.problem.ext[:vars][node])
            set_start_value(
                deteq.problem.ext[:vars][node][variable],
                JuMP.value(temp[string(variable)]),
            )
        end
    end
end

"""
    add_to_dictionary!(
        original::Dict{AbstractTree,Dict{Symbol,Any}},
        add::T where {T<:Dict},
        sym::Union{Symbol,Vector{Symbol}})

Given a dictionary produced by the `solution_to_dictionary()` function, and a Symbol or Symbol[],
this function adds that/those Symbol(s) to the dictionary.

### Required Arguments
`original` is a dictionary produced by `solution_to_dictionary()`
`add` is the dictionary (with the keys being `AbstractTree` objects) that you wish to add to
the `original` dictionary
`sym` is a Symbol or Symbol vector of keys to add to each node within the dictionary
"""
function add_to_dictionary!(
    original::Dict{AbstractTree,Dict{Symbol,Any}},
    add::T where {T<:Dict},
    sym::Symbol,
)
    no_warn = false
    if length(keys(original)) == 0
        no_warn = true
    end

    for k in keys(add)
        if haskey(original, k)
            original[k][sym] = copy(add[k])
        else
            if !no_warn
                @warn("Key: " * string(k) * " not in original dictionary.")
            end
            original[k] = Dict{Symbol,Any}()
            original[k][sym] = copy(add[k])
        end
    end
end

"""
    remove_from_dictionary!(
        original::Dict{AbstractTree,Dict{Symbol,Any}},
        sym::Union{Symbol,Vector{Symbol}})

Given a dictionary produced by the `solution_to_dictionary()` function, and a Symbol or Symbol[],
this function removes that/those Symbol(s) from the dictionary.

### Required Arguments
`original` is a dictionary produced by `solution_to_dictionary()`
`sym` is a Symbol or Symbol vector of keys to remove from each node within the dictionary.
"""
function remove_from_dictionary!(
    original::Dict{AbstractTree,Dict{Symbol,Any}},
    sym::Union{Symbol,Vector{Symbol}},
)
    if typeof(sym) == Symbol
        sym = [sym]
    end

    for data in values(original)
        for s in sym
            if haskey(data, s)
                delete!(data, s)
            end
        end
    end
end

function combine_dictionaries(
    dict1::Dict{AbstractTree,Dict{Symbol,Any}},
    dict2::Dict{AbstractTree,Dict{Symbol,Any}},
)
    temp = Dict{AbstractTree,Dict{Symbol,Any}}()

    for node in keys(dict1)
        if !haskey(temp, node)
            temp[node] = Dict{Symbol,Any}()
        end
        for sym in keys(dict1[node])
            if haskey(temp[node], sym)
                error(
                    "Key " *
                    string(sym) *
                    " exists in both dictionaries for the same nodes.",
                )
            end
            temp[node][sym] = copy(dict1[node][sym])
        end
    end
    for node in keys(dict2)
        if !haskey(temp, node)
            temp[node] = Dict{Symbol,Any}()
        end
        for sym in keys(dict2[node])
            if haskey(temp[node], sym)
                error(
                    "Key " *
                    string(sym) *
                    " exists in both dictionaries for the same nodes.",
                )
            end
            temp[node][sym] = copy(dict2[node][sym])
        end
    end
    return temp
end

"""
	get_active_columns(jmodel::JuDGEModel; inttol = 10^-7)

Returns a list of tuples of non-zero columns along with their corresponding value.

### Required Arguments
`jmodel` is a solved JuDGEModel

### Optional Arguments
`inttol` is the minimum value for a variable to be deemed non-zero
"""
function get_active_columns(jmodel::JuDGEModel; inttol = 10^-7)
    active = Dict{AbstractTree,Vector{Any}}()
    for node in collect(jmodel.tree)
        active[node] = []
        for col in jmodel.master_problem.ext[:columns][node]
            if JuMP.value(col.var) > inttol
                push!(active[node], (col, JuMP.value(col.var)))
            end
        end
    end

    return active
end

function overprint(str)
    print("\e[2K")
    print("\e[1G")
    return print(str)
end

function printleft(str)
    print("\u1b[40G")
    print("\u1b[1K")
    print("\u1b[1G")
    return print(str)
end

function printright(str)
    print("\u1b[41G")
    print("\u1b[0K")
    return print(str)
end
