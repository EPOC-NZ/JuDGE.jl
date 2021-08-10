###
# asserts for constructor of JuDGEModel
###

# Definition of an important set of constraints
function NormalConstraints()
    return [
        (GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64}),
        (GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64}),
        (GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64}),
        (VariableRef, MOI.LessThan{Float64}),
        (VariableRef, MOI.GreaterThan{Float64}),
        (VariableRef, MOI.EqualTo{Float64}),
    ]
end

function check_specification_is_legal(sub_problems)
    message = "JuDGE Specification Error: "
    if !same_expansions_at_each_node(sub_problems)
        message *= "Every subproblem must have the same expansion variables"
        return error(message)
    end
    if !same_shutdowns_at_each_node(sub_problems)
        message *= "Every subproblem must have the same shutdown variables"
        return error(message)
    end
    if !check_objectives(sub_problems)
        message *= "The subproblem objective must be minimizing and either be of type VariableRef or GenericAffExpr"
        return error(message)
    end

    check_costs(sub_problems)
    return check_constraints(sub_problems) #must be at the end of this subroutine!!
end

#### Expansion variables at every node have to be the "same"
function same_expansions_at_each_node(subproblems)
    expansion_structures =
        collect(map(get_structure_of_expansions, values(subproblems)))
    return all(map(x -> x == expansion_structures[1], expansion_structures))
end

function same_shutdowns_at_each_node(subproblems)
    shutdown_structures =
        collect(map(get_structure_of_shutdowns, values(subproblems)))
    return all(map(x -> x == shutdown_structures[1], shutdown_structures))
end

function check_objectives(subproblems)
    for (n, sp) in subproblems
        if !(objective_function_type(sp) <: GenericAffExpr) &&
           objective_function_type(sp) != VariableRef
            return false
        elseif objective_sense(sp) != MOI.MIN_SENSE
            return false
        end
    end
    return true
end

function check_costs(subproblems)
    for (n, sp) in subproblems
        check_sp_costs(sp)
    end
end

function check_sp_costs(model)
    if haskey(model.ext, :capitalcosts) &&
       model.ext[:capitalcosts][:constant] != 0
        #if typeof(model.ext[:capitalcosts])!=AffExpr
        #    error("@capitalcosts must be provided a linear expression (AffExpr)")
        #elseif model.ext[:capitalcosts][:constant]!=0
        error("@capitalcosts should not contain any constant terms")
        # else
        #     av=all_variables(model)
        #     for v in av
        #         found=false
        #         for e in model.ext[:expansions]
        #             if v==e[2] || v in e[2]
        #                 found=true
        #                 break
        #             end
        #         end
        #         if !found && v in keys(model.ext[:capitalcosts].terms)
        #             error("@capitalcosts should only contain expansion variables"*string(v))
        #         end
        #     end
        # end
    elseif haskey(model.ext, :ongoingcosts) &&
           model.ext[:ongoingcosts][:constant] != 0
        #if typeof(model.ext[:ongoingcosts])!=AffExpr
        #    error("@ongoingcosts must be provided a linear expression (AffExpr)")
        #elseif model.ext[:ongoingcosts].constant!=0
        error("@ongoingcosts should not contain any constant terms")
        # else
        #     av=all_variables(model)
        #     for v in av
        #         found=false
        #         for e in model.ext[:expansions]
        #             if v==e[2] || v in e[2]
        #                 found=true
        #                 break
        #             end
        #         end
        #         if !found && v in keys(model.ext[:ongoingcosts].terms)
        #             error("@ongoingcosts should only contain expansion variables")
        #         end
        #     end
        # end
    end
    return nothing
end

function get_structure_of_expansions(subproblem)
    expkeys = get_expansion_keys(subproblem)

    function thinwrapper(var::AbstractArray)
        collect(eachindex(var))
    end
    function thinwrapper(var)
        ()
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in expkeys)
end

function get_structure_of_shutdowns(subproblem)
    shutkeys = get_shutdown_keys(subproblem)

    function thinwrapper(var::AbstractArray)
        collect(eachindex(var))
    end
    function thinwrapper(var)
        ()
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in shutkeys)
end

function get_expansion_keys(model)
    filter(keys(model.obj_dict)) do key
        for (exp, var) in model.ext[:expansions]
            if var === model.obj_dict[key] &&
               model.ext[:options][exp][1] == :expansion
                return true
            end
        end
        return false
    end
end

function get_shutdown_keys(model)
    filter(keys(model.obj_dict)) do key
        for (exp, var) in model.ext[:expansions]
            if var === model.obj_dict[key] &&
               model.ext[:options][exp][1] == :shutdown
                return true
            end
        end
        return false
    end
end

function check_constraints(subproblems)
    for (n, sp) in subproblems
        check_sp_constraints(sp)
    end
end

function check_sp_constraints(model)
    loct = list_of_constraint_types(model)
    expansion_keys = get_expansion_keys(model)
    shutdown_keys = get_shutdown_keys(model)

    all_expansion_variables = Vector{VariableRef}()
    all_shutdown_variables = Vector{VariableRef}()

    for exp_key in expansion_keys
        var = model[exp_key]
        if model.ext[:options][exp_key][4] != :Con
            if typeof(var) == VariableRef
                push!(all_expansion_variables, var)
            elseif typeof(var) <: AbstractArray
                for i in eachindex(var)
                    push!(all_expansion_variables, var[i])
                end
            end
        end
    end

    for shut_key in shutdown_keys
        var = model[shut_key]
        if typeof(var) == VariableRef
            push!(all_shutdown_variables, var)
        elseif typeof(var) <: AbstractArray
            for i in eachindex(var)
                push!(all_shutdown_variables, var[i])
            end
        end
    end

    warnings = [
        "Positive coefficient for shutdown variable on LHS of >= constraint",
        "Positive coefficient for expansion variable on LHS of <= constraint",
        "Expansion or shutdown variable in == constraint",
        "Negative coefficient for expansion variable on LHS of >= constraint",
        "Negative coefficient for shutdown variable on LHS of <= constraint",
    ]

    status = [false, false, false, false, false]

    for ct in loct
        for con in all_constraints(model, ct[1], ct[2])
            con_obj = JuMP.constraint_object(con)
            if typeof(con_obj.func) == VariableRef
                # if typeof(con_obj.set)==MOI.GreaterThan{Float64}
                #     if con_obj.func in all_shutdown_variables
                #         status[1]=true
                #     end
                # elseif typeof(con_obj.set)==MOI.LessThan{Float64}
                #     if con_obj.func in all_expansion_variables
                #        status[2]=true
                #     end
                if typeof(con_obj.set) == MOI.EqualTo{Float64}
                    if con_obj.func in all_expansion_variables ||
                       con_obj.func in all_shutdown_variables
                        status[3] = true
                    end
                end
            elseif typeof(con_obj.func) <: GenericAffExpr
                if typeof(con_obj.set) == MOI.GreaterThan{Float64}
                    for (v, c) in con_obj.func.terms
                        if c > 0.0 && v in all_shutdown_variables
                            status[1] = true
                        elseif c < 0.0 && v in all_expansion_variables
                            status[4] = true
                        end
                    end
                elseif typeof(con_obj.set) == MOI.LessThan{Float64}
                    for (v, c) in con_obj.func.terms
                        if c < 0.0 && v in all_shutdown_variables
                            status[5] = true
                        elseif c > 0.0 && v in all_expansion_variables
                            status[2] = true
                        end
                    end
                elseif typeof(con_obj.set) == MOI.EqualTo{Float64}
                    for (v, c) in con_obj.func.terms
                        if v in all_shutdown_variables ||
                           v in all_expansion_variables
                            status[3] = true
                        end
                    end
                end
            elseif typeof(con_obj.func) <: GenericQuadExpr
                LHS = QuadExpr(
                    AffExpr(0.0),
                    OrderedDict{UnorderedPair{VariableRef},Float64}(),
                )
                for (v, c) in con_obj.func.terms
                    if v.a in all_shutdown_variables ||
                       v.b in all_shutdown_variables ||
                       v.a in all_expansion_variables ||
                       v.b in all_expansion_variables
                        error(
                            "JuDGE Specification Error: Expansion or shutdown variable in quadratic term",
                        )
                    end
                end

                if typeof(con_obj.set) == MOI.GreaterThan{Float64}
                    for (v, c) in con_obj.func.aff.terms
                        if c > 0.0 && v in all_shutdown_variables
                            status[1] = true
                        elseif c < 0.0 && v in all_expansion_variables
                            status[4] = true
                        end
                    end
                elseif typeof(con_obj.set) == MOI.LessThan{Float64}
                    for (v, c) in con_obj.func.aff.terms
                        if c < 0.0 && v in all_shutdown_variables
                            status[5] = true
                        elseif c > 0.0 && v in all_expansion_variables
                            status[2] = true
                        end
                    end
                elseif typeof(con_obj.set) == MOI.EqualTo{Float64}
                    for (v, c) in con_obj.func.aff.terms
                        if v in all_shutdown_variables ||
                           v in all_expansion_variables
                            status[3] = true
                        end
                    end
                end
            elseif typeof(con_obj.func) ==
                   Vector{GenericAffExpr{Float64,VariableRef}}
                for aff in con_obj.func
                    for (v, c) in aff.terms
                        if v in all_shutdown_variables ||
                           v in all_expansion_variables
                            error(
                                "JuDGE Specification Error: Expansion or shutdown variable in indicator constraint",
                            )
                        end
                    end
                end
            else
                @warn(
                    "JuDGE Specification Error: Unable to verify constraint of type " *
                    string(typeof(con_obj.func))
                )
            end
        end
    end

    for i in 1:5
        if status[i]
            @warn(warnings[i])
        end
    end

    return nothing
end
