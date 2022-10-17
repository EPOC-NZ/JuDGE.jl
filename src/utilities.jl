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

    if start_value(x) !== nothing
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

function key_to_tuple(key::Any)
    if typeof(key) <: JuMP.Containers.DenseAxisArrayKey ||
       typeof(key) <: CartesianIndex
        return length(key.I) == 1 ? key.I[1] : key.I
    else
        return key
    end
end

function key_to_string(key::Union{Int,Tuple,Symbol,AbstractString})
    #return replace(string(key), ")" => "", "(" => "", "\"" => "", ", " => ",") # needs Julia 1.7
    str = replace(string(key), ")" => "")
    str = replace(str, "(" => "")
    str = replace(str, "\"" => "")
    return replace(str, ", " => ",")
end

function get_keys(var::Union{AbstractArray,Dict})
    return typeof(var) <: JuMP.Containers.SparseAxisArray ? keys(var.data) :
           keys(var)
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

"""
    compute_objval(
        scenarios::Dict{Leaf,Float64},
        probabilities::Dict{AbstractTree,Float64},
        risk::Union{Risk,Vector{Risk}},
    )

This function finds the risk-adjusted objective value, given scenario objectives, a probability distribution
and a risk measure.

### Required arguments

`scenarios` is a dictionary mapping leaf nodes to scenario objectives.

`probabilities` is a dictionary mapping each leaf nodes to the probability of reaching that node.

`risk` is a `JuDGE.Risk` object or a vector of such objects.
"""
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
    for i in eachindex(risk)
        EV_weight -= risk[i].λ
        so = scenario_objs
        if risk[i].offset !== nothing
            for j in eachindex(so)
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
    function get_risk_probs(
        model::Union{JuDGEModel,DetEqModel},
        mode::Symbol = :marginal,
    )
Returns a dictionary mapping the nodes in a tree to a probability (either marginal or conditional).

### Required arguments

`model` is a JuDGE or DetEq model for which we wish to find the risk-adjusted probabilities for each scenario.

`mode` is either `:marginal` or `:conditional`.
"""
function get_risk_probs(
    model::Union{JuDGEModel,DetEqModel},
    mode::Symbol = :marginal,
)
    m = typeof(model) == JuDGEModel ? model.master_problem : model.problem

    cons = all_constraints(m, AffExpr, MOI.EqualTo{Float64})

    if mode ∉ [:marginal, :conditional]
        error(
            "Invalid probability mode. mode must be :conditional or :marginal",
        )
    end

    if has_duals(m)
        π = Dict{Leaf,Float64}()
        leafnodes = JuDGE.get_leafnodes(model.tree)
        if typeof(model) == JuDGEModel
            for leaf in leafnodes
                π[leaf] = -dual(model.master_problem.ext[:scenprofit_con][leaf])
            end
        else
            for leaf in leafnodes
                π[leaf] = -dual(model.problem.ext[:scenario_con][leaf])
            end
        end

        pr = Dict{AbstractTree,Float64}()
        nodes = reverse(collect(model.tree))
        for node in nodes
            if typeof(node) == Leaf
                pr[node] = π[node]
            else
                pr[node] = sum(pr[n] for n in node.children)
            end
        end

        if mode == :conditional
            for node in nodes
                if node.parent !== nothing && pr[node.parent] > 0
                    pr[node] /= pr[node.parent]
                else
                    delete!(pr, node)
                end
            end
            return pr
        else
            return pr
        end
    else
        @warn(
            "Dual variables not present; computing approximate probabilities."
        )
    end

    return compute_risk_probs(
        get_scen_objs(model),
        model.probabilities,
        model.risk,
        mode,
    )
end

"""
    function compute_risk_probs(
        scenarios::Dict{Leaf,Float64},
        probabilities::Dict{AbstractTree,Float64},
        risk::Union{Risk,Vector{Risk}},
        mode::Symbol = :marginal,
    )
Returns a dictionary mapping the nodes in a tree to a probability (either marginal or conditional).

Note that this method is only approximate and may not be able to provide useful insights when there is degeneracy.

### Required arguments

`scenarios` is a dictionary mapping leaf nodes to scenario costs.

`probabilities` is a dictionary mapping nodes to marginal probabilities.

`risk` is a `JuDGE.Risk` object, or a vector of such objects.

`mode` is either `:marginal` or `:conditional`.
"""
function compute_risk_probs(
    scenarios::Dict{Leaf,Float64},
    probabilities::Dict{AbstractTree,Float64},
    risk::Union{Risk,Vector{Risk}},
    mode::Symbol = :marginal,
)
    scenario_objs = Vector{Tuple{Float64,Float64,Leaf}}()
    riskprob = Dict{Leaf,Float64}()

    for (leaf, val) in scenarios
        push!(scenario_objs, (probabilities[leaf], val, leaf))
        riskprob[leaf] = 0
    end
    EV_weight = 1.0

    if typeof(risk) == Risk
        if risk.α == 1.0 || risk.λ == 0.0
            risk = []
        else
            risk = [risk]
        end
    end
    for i in eachindex(risk)
        EV_weight -= risk[i].λ
        so = copy(scenario_objs)
        if risk[i].offset !== nothing
            for j in eachindex(so)
                so[j] =
                    (so[j][1], so[j][2] - risk[i].offset[so[j][3]], so[j][3])
            end
        end
        sort!(so, by = i -> i[2], rev = true)
        beta = risk[i].α
        for scen in so
            if scen[1] > beta
                riskprob[scen[3]] += risk[i].λ * beta / risk[i].α
                beta = 0
            else
                riskprob[scen[3]] += risk[i].λ * scen[1] / risk[i].α
                beta -= scen[1]
            end
        end
    end

    for i in keys(scenarios)
        riskprob[i] += EV_weight * probabilities[i]
    end

    node = collect(keys(probabilities))[1]
    while node.parent !== nothing
        node = node.parent
    end
    nodes = reverse(collect(node))

    pr = Dict{AbstractTree,Float64}()
    for node in nodes
        if typeof(node) == Leaf
            pr[node] = riskprob[node]
        else
            pr[node] = sum(pr[n] for n in node.children)
        end
    end

    if mode == :conditional
        for node in nodes
            if node.parent !== nothing && pr[node.parent] > 0.0
                pr[node] /= pr[node.parent]
            else
                delete!(pr, node)
            end
        end
        return pr
    elseif mode == :marginal
        return pr
    else
        error(
            "Invalid probability mode. mode must be :conditional or :marginal",
        )
    end
end

"""
    get_regret(
        model::Union{JuDGEModel,DetEqModel},
        baseline::Union{JuDGEModel,DetEqModel},
    )

The function takes a pair of solved `JuDGEModel` or `DetEqModel` objects that use the same `JuDGE.Tree`,
and compares the cost in each scenario. A dictionary of regret for each leaf node of the tree is returned.

### Required Arguments
`model` is the model whose solution is the policy that we may wish to implement.

`baseline` is the same model, but using a different policy (it was solved using a different risk measure).
"""
function get_regret(
    model::Union{JuDGEModel,DetEqModel},
    baseline::Union{JuDGEModel,DetEqModel},
)
    if model.tree != baseline.tree
        error(
            "This function (currently) requires both policies to use the same trees when computing regret.",
        )
    end

    objs1 = get_scen_objs(model)
    objs2 = get_scen_objs(baseline)

    regret = Dict{Leaf,Float64}()
    for leaf in get_leafnodes(model.tree)
        regret[leaf] = objs1[leaf] - objs2[leaf]
    end
    return regret
end

"""
    solution_to_dictionary(model::Union{JuDGEModel,DetEqModel}; prefix::String = "")

Create a nested dictionary with the solution values for each node copied across to
a standardised structure.

### Required Arguments
`model` is a solved JuDGEModel or DetEqModel

### Optional Arguments
`prefix` is a string that will be prepended to each of the variable names (e.g. if comparing to versions of the same model)
"""
function solution_to_dictionary(
    model::Union{JuDGEModel,DetEqModel};
    prefix::String = "",
)
    function helper(node::AbstractTree, solution::T where {T<:Dict})
        solution[node] = Dict{Symbol,Any}()

        var_groups =
            typeof(model) == DetEqModel ? model.problem.ext[:all_vars][node] :
            model.sub_problems[node].ext[:all_vars]

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

                for key in get_keys(vals)
                    strkey = key_to_string(key_to_tuple(key))
                    solution[node][sym][strkey] = vals[key]
                end
            elseif typeof(vars) <: Dict
                skip = false
                for v in values(vars)
                    if typeof(v) ∉ [VariableRef, AffExpr]
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

                for key in keys(vars)
                    strkey = key_to_string(key_to_tuple(key))
                    solution[node][sym][strkey] = value(vars[key])
                end
            end
        end

        if typeof(node) == Tree
            for child in node.children
                helper(child, solution)
            end
        elseif typeof(model) == DetEqModel
            solution[node][Symbol(prefix * "scenario_obj")] =
                JuMP.value(model.problem.ext[:scenario_obj][node])
        else
            solution[node][Symbol(prefix * "scenario_obj")] =
                JuMP.value(model.master_problem.ext[:scenprofit_var][node])
        end
        return solution
    end

    if typeof(model) == DetEqModel
        if termination_status(model.problem) != MOI.OPTIMAL &&
           termination_status(model.problem) != MOI.TIME_LIMIT &&
           termination_status(model.problem) != MOI.INTERRUPTED
            error("You need to first solve the deterministic equivalent model.")
        end
    else
        if termination_status(model.master_problem) != MOI.OPTIMAL &&
           termination_status(model.master_problem) != MOI.INTERRUPTED &&
           termination_status(model.master_problem) != MOI.TIME_LIMIT &&
           termination_status(model.master_problem) != MOI.LOCALLY_SOLVED
            error("You need to first solve the decomposed model.")
        end
    end

    if prefix != ""
        prefix *= "_"
    end

    solution = Dict{AbstractTree,Dict{Symbol,Any}}()
    return helper(model.tree, solution)
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
            if typeof(jmodel.master_problem.ext[:expansions][node][name]) ==
               VariableRef
                set_start_value(
                    deteq.problem.ext[:master_vars][node][name],
                    JuMP.value(
                        jmodel.master_problem.ext[:expansions][node][name],
                    ),
                )
            else
                for index in
                    get_keys(jmodel.master_problem.ext[:expansions][node][name])
                    key = key_to_tuple(index)
                    set_start_value(
                        deteq.problem.ext[:master_vars][node][name][key],
                        JuMP.value(
                            jmodel.master_problem.ext[:expansions][node][name][key],
                        ),
                    )
                end
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
        sym::Symbol,
    )

Given a dictionary produced by the `solution_to_dictionary()` function, and a Symbol or Symbol[],
this function adds that/those Symbol(s) to the dictionary.

### Required Arguments
`original` is a dictionary produced by `solution_to_dictionary()`

`add` is the dictionary (with the keys being `AbstractTree` objects) that you wish to add to
the `original` dictionary

`sym` is a Symbol to use as the key within the dictionary
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

"""
    combine_dictionaries(
        dict1::Dict{AbstractTree,Dict{Symbol,Any}},
        dict2::Dict{AbstractTree,Dict{Symbol,Any}},
    )

This function combines two dictionaries that have been created using `solution_to_dictionary(...)`.
Note that the `prefix` argument should be used so that different keys are created in each dictionary.

### Required Arguments
`dict1` is a dictionary created using `solution_to_dictionary(...)`.

`dict2` is a dictionary created using `solution_to_dictionary(...)`.
"""
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

"""
    get_scen_objs(jmodel::JuDGEModel)

This function creates a dictionary that contains the scenario profit for each leaf node.

### Required Arguments
    `jmodel` is a solved `JuDGEModel`.
"""
function get_scen_objs(jmodel::JuDGEModel)
    offset = Dict{Leaf,Float64}()
    for leaf in JuDGE.get_leafnodes(jmodel.tree)
        offset[leaf] = value(jmodel.master_problem.ext[:scenprofit_var][leaf])
    end
    return offset
end

"""
    get_scen_objs(deteq::DetEqModel)

This function creates a dictionary that contains the scenario profit for each leaf node.

### Required Arguments
    `deteq` is a solved `DetEqModel`.
"""
function get_scen_objs(deteq::DetEqModel)
    offset = Dict{Leaf,Float64}()
    for leaf in JuDGE.get_leafnodes(deteq.tree)
        offset[leaf] = value(deteq.problem.ext[:scenario_obj][leaf])
    end
    return offset
end

"""
    scenarios_CDF(model::Union{JuDGEModel,DetEqModel}; tol::Float64 = 1e-8)

Given a solved `JuDGEModel` or `DetEqModel`, this function will return a set of points
for a cumulative distribution function for the scenario costs.

### Required Arguments
`model` is either a solved `JuDGEModel` or `DetEqModel`.

### Optional Arguments
`tol` is the tolerance used to reduce the number of distinct objective values.
"""
function scenarios_CDF(model::Union{JuDGEModel,DetEqModel}; tol::Float64 = 1e-8)
    scenobj = Tuple{Float64,Float64}[]
    if typeof(model) == JuDGEModel
        for leaf in JuDGE.get_leafnodes(model.tree)
            push!(
                scenobj,
                (
                    value(model.master_problem.ext[:scenprofit_var][leaf]),
                    model.probabilities[leaf],
                ),
            )
        end
    else
        for leaf in JuDGE.get_leafnodes(model.tree)
            push!(
                scenobj,
                (
                    value(model.problem.ext[:scenario_obj][leaf]),
                    model.probabilities[leaf],
                ),
            )
        end
    end
    sort!(scenobj)

    total = 0.0
    for i in eachindex(scenobj)
        total += scenobj[i][2]
        scenobj[i] = (scenobj[i][1], total)
    end

    scenobj2 = Tuple{Float64,Float64}[]

    prev = 0.0
    prev_obj = -Inf
    for i in eachindex(scenobj)
        if i < length(scenobj) && scenobj[i+1][1] - scenobj[i][1] < tol
            continue
        end
        push!(scenobj2, (scenobj[i][1], prev))
        prev = scenobj[i][2]
        push!(scenobj2, (scenobj[i][1], prev))
    end

    return scenobj2
end

function Base.getindex(jmodel::JuDGEModel, node::AbstractTree)
    return jmodel.sub_problems[node].ext[:all_vars]
end

function Base.getindex(deteq::DetEqModel, node::AbstractTree)
    return deteq.problem.ext[:all_vars][node]
end

# function create_dash(
#     blocks::Union{Nothing,Vector{DashWrapper.Block}},
#     callbacks::Union{
#         Nothing,
#         DashWrapper.DashCallback,
#         Vector{DashWrapper.DashCallback},
#     };
#     path::String = pwd(),
#     title::String = "JuDGE Solution Dashboard",
#     data::Any = Dict(),
#     extra_files::Vector{String} = String[],
#     external_js::Vector{String} = String[],
#     external_css::Vector{String} = String[],
# )
#     extra_files = [
#         [
#             joinpath(dirname(@__DIR__), "dash", "tree.js"),
#             joinpath(dirname(@__DIR__), "dash", "judge-dash.js"),
#             joinpath(dirname(@__DIR__), "dash", "colors.js"),
#             joinpath(
#                 dirname(@__DIR__),
#                 "dash",
#                 "jsonview",
#                 "jsonview.bundle.css",
#             ),
#             joinpath(
#                 dirname(@__DIR__),
#                 "dash",
#                 "jsonview",
#                 "jsonview.bundle.js",
#             ),
#             joinpath(
#                 dirname(@__DIR__),
#                 "dash",
#                 "svg-pan-zoom",
#                 "svg-pan-zoom.min.js",
#             ),
#         ]
#         extra_files
#     ]
#
#     return app = DashWrapper.create_dash(
#         blocks,
#         callbacks,
#         path = path,
#         title = title,
#         data = data,
#         extra_files = extra_files,
#         external_js = external_js,
#         external_css = external_css,
#     )
# end
#
# function dash_layout()
#     blocks = DashWrapper.Block[]
#     push!(
#         blocks,
#         DashWrapper.Block(
#             (0, 0),
#             (1, 5),
#             "JuDGE Scenario Tree",
#             div_name = "tree-div",
#         ),
#     )
#     push!(
#         blocks,
#         DashWrapper.Block(
#             (0, 5),
#             (1, 1),
#             "Settings",
#             div_name = "settings-div",
#         ),
#     )
#     return blocks
# end
#
# function export_tree(
#     some_tree::AbstractTree;
#     scale_edges = nothing,
#     scale_nodes::Float64 = 0.0,
#     max_size::Float64 = 50.0,
#     truncate::Int = -1,
#     rel_angle::Bool = false,
#     style::Symbol = :standard,
#     box_size::Int = 800,
#     skip_root::Bool = false,
#     data::Dict = Dict(),
# )
#     maxdata = Dict{Symbol,Any}()
#     mindata = Dict{Symbol,Any}()
#     scale_range = Dict{Symbol,Dict{Symbol,Any}}()
#
#     function data_scale(node)
#         #        if data_from_original
#         #            node = getID(node)
#         #        end
#         if skip_root && node.parent == nothing
#             return
#         end
#         for sym in keys(data[node])
#             if typeof(data[node][sym]) == Float64
#                 if !haskey(scale_range, sym)
#                     scale_range[sym] = Dict{Symbol,Any}()
#                     scale_range[sym][:min] = data[node][sym]
#                     scale_range[sym][:max] = data[node][sym]
#                 else
#                     if data[node][sym] < scale_range[sym][:min]
#                         scale_range[sym][:min] = data[node][sym]
#                     elseif data[node][sym] > scale_range[sym][:max]
#                         scale_range[sym][:max] = data[node][sym]
#                     end
#                 end
#             elseif typeof(data[node][sym]) <: Dict
#                 for key in keys(data[node][sym])
#                     if !haskey(scale_range, sym)
#                         scale_range[sym] = Dict{Symbol,Any}()
#                         scale_range[sym][:max] = data[node][sym][key]
#                         scale_range[sym][:min] = data[node][sym][key]
#                     else
#                         if data[node][sym][key] < scale_range[sym][:min]
#                             scale_range[sym][:min] = data[node][sym][key]
#                         elseif data[node][sym][key] > scale_range[sym][:max]
#                             scale_range[sym][:max] = data[node][sym][key]
#                         end
#                     end
#                 end
#             elseif typeof(data[node][sym]) <: Array
#                 for val in data[node][sym]
#                     if !haskey(scale_range, sym)
#                         scale_range[sym] = Dict{Symbol,Any}()
#                         scale_range[sym][:max] = val
#                         scale_range[sym][:min] = val
#                     else
#                         if val < scale_range[sym][:min]
#                             scale_range[sym][:min] = val
#                         elseif val > scale_range[sym][:max]
#                             scale_range[sym][:max] = val
#                         end
#                     end
#                 end
#             end
#         end
#     end
#
#     function arc_json(node, parent)
#         return Dict("from" => get_id[parent] + 1, "to" => get_id[node] + 1)
#     end
#
#     get_id = Dict{AbstractTree,Int}()
#
#     angles = Dict{AbstractTree,Float64}()
#     position = Dict{AbstractTree,Tuple{Float64,Float64}}()
#     position2 = Dict{AbstractTree,Tuple{Float64,Float64}}()
#
#     function setpositions(
#         node::AbstractTree,
#         rel_angle::Bool,
#         l::Float64,
#         scale::Float64,
#         odd::Bool,
#     )
#         if typeof(node) == Leaf
#             return
#         end
#         a = -2 * pi / (length(node.children))
#         current = 0.0
#         if rel_angle
#             current = angles[node] + pi + a / 2#0.0
#         elseif (length(node.children) % 2) == 1
#             current +=
#                 pi / 2 +
#                 (length(node.children) - 1) * pi / length(node.children)
#         elseif odd
#             current += 3 * pi / 2 - pi / length(node.children)
#         else
#             current += 3 * pi / 2 - 2 * pi / length(node.children)
#         end
#
#         for child in node.children
#             angles[child] = current
#             position[child] = (
#                 position[node][1] + l * cos(current),
#                 position[node][2] - l * sin(current),
#             )
#             if length(node.children) == 2
#                 if odd
#                     setpositions(child, rel_angle, l, scale, !odd)
#                 else
#                     setpositions(child, rel_angle, l * scale^2, scale, !odd)
#                 end
#             else
#                 setpositions(child, rel_angle, l * scale, scale, !odd)
#             end
#             current += a
#         end
#     end
#
#     function setpositions2(node::AbstractTree, leaf_sep::Float64)
#         function locate(node::AbstractTree, vert::Float64, horz::Float64)
#             if typeof(node) == Leaf ||
#                (truncate != -1 && JuDGE.depth(node) >= truncate)
#                 vert += leaf_sep
#                 position2[node] = (horz, vert)
#                 return (vert, vert)
#             else
#                 verts = Float64[]
#                 for child in node.children
#                     pos, vert = locate(child, vert, horz + parch_sep)
#                     push!(verts, pos)
#                 end
#                 position2[node] = (horz, sum(verts) / length(verts))
#                 return (position2[node][2], vert)
#             end
#         end
#         num_leaf = length(JuDGE.get_leafnodes(node, truncate = truncate))
#         max_depth = JuDGE.depth(collect(node, truncate = truncate)[end])
#         parch_sep = 0.8 * leaf_sep * num_leaf / max_depth
#         return locate(node, 0.0, 0.0)
#     end
#
#     scale_factors = [
#         [1.0],
#         [1.0, 0.87, 0.83, 0.78, 0.74, 0.71, 0.695, 0.685],
#         [1.0, 0.65, 0.52, 0.48],
#         [1.0, 0.45, 0.42, 0.42, 0.41, 0.4],
#         [1.0, 0.44, 0.37, 0.36],
#         [1.0, 0.42, 0.34, 0.33],
#         [1.0, 0.35, 0.3],
#         [1.0, 0.3, 0.26],
#         [1.0, 0.27, 0.23],
#         [1.0, 0.24, 0.22],
#     ]
#
#     if scale_edges == nothing
#         if typeof(some_tree) == Leaf
#             scale_edges = 1.0
#         else
#             dg = length(some_tree.children)
#             if dg <= 10
#                 dp = min(
#                     JuDGE.depth(collect(some_tree, truncate = truncate)[end]),
#                     length(scale_factors[dg]),
#                 )
#                 scale_edges = scale_factors[dg][dp]
#             else
#                 dp = JuDGE.depth(collect(some_tree, truncate = truncate)[end])
#                 scale_edges = 0.22 * 0.98^(dg^dp - 11)
#             end
#         end
#     end
#
#     angles[some_tree] = rel_angle ? 0.0 : -pi / 2
#     position[some_tree] = (0.0, 0.0)
#
#     setpositions(some_tree, rel_angle, 700.0, scale_edges, true)
#     setpositions2(some_tree, 70.0)
#     index = 0
#     if skip_root
#         index -= 1
#     end
#
#     arcs = Dict[]
#     nodes = Dict[]
#
#     for node in collect(some_tree, truncate = truncate)
#         get_id[node] = index
#         parent = node.parent
#         if parent != nothing && (!skip_root || parent.parent != nothing)
#             push!(arcs, arc_json(node, parent))
#         end
#         index += 1
#     end
#
#     for node in collect(some_tree, truncate = truncate)
#         if !skip_root || node != some_tree
#             data_scale(node)
#             temp = Dict{Symbol,Any}()
#             temp[:id] = get_id[node] + 1
#             temp[:label] = node.name
#             temp[:level] = JuDGE.depth(node)
#             temp[:leaf] =
#                 typeof(node) == Leaf ||
#                 (truncate != -1 && JuDGE.depth(node) >= truncate)
#             temp[:posX] = position[node][1]
#             temp[:posY] = position[node][2]
#             temp[:posX2] = position2[node][1]
#             temp[:posY2] = position2[node][2]
#             if haskey(data, node)
#                 temp[:data] = data[node]
#             end
#             push!(nodes, temp)
#         end
#     end
#
#     if scale_nodes == 0.0
#         scale_nodes = min(
#             1.0,
#             exp(
#                 log(
#                     (
#                         400 *
#                         scale_edges^JuDGE.depth(
#                             collect(some_tree, truncate = truncate)[end],
#                         )
#                     ) / max_size,
#                 ) /
#                 JuDGE.depth(collect(some_tree, truncate = truncate)[end]),
#             ),
#         )
#     end
#     min_size = 25
#
#     return Dict(
#         "nodes" => nodes,
#         "arcs" => arcs,
#         "node_scale" => scale_nodes,
#         "min_size" => min_size,
#         "max_size" => max_size,
#         "scale" => scale_range,
#     )
# end

function run_file(filename)
    if Sys.iswindows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    elseif Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
    else
        error("Unable to show plot. Try opening the file $(filename) manually.")
    end
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
