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
                        if length(key.I) == 1
                            strkey = string(key.I[1])
                        else
                            strkey = string(key.I)
                            strkey = replace(strkey, ")" => "")
                            strkey = replace(strkey, "(" => "")
                            strkey = replace(strkey, "\"" => "")
                            strkey = replace(strkey, ", " => ",")
                        end
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
                    strkey = ss[2][1:end-1]
                    strkey = replace(strkey, ")" => "")
                    strkey = replace(strkey, "(" => "")
                    strkey = replace(strkey, ", " => ",")
                    solution[node][Symbol(prefix * ss[1])][strkey] =
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
                # @constraint(deteq.problem,deteq.problem.ext[:master_vars][node][name][key]==
                #     JuMP.value(jmodel.master_problem.ext[:expansions][node][name][index])
                # )
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
            # @constraint(deteq.problem,deteq.problem.ext[:vars][node][variable]==JuMP.value(temp[string(variable)]))
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
