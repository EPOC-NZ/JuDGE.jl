# Builds a deterministic equivalent by joining the subproblems into
# a single MIP, then adding the constraints that there can only be
# a single investment along any path in the tree.
using DataStructures

struct DetEqModel
    problem::JuMP.Model
    tree::AbstractTree
    probabilities::Dict{AbstractTree,Float64}
    risk::Union{Risk,Vector{Risk}}
end

"""
	DetEqModel(tree::AbstractTree,
               probabilities,
               sub_problem_builder::Function,
               solver
               discount_factor=1.0,
               risk=RiskNeutral,
               sideconstraints=nothing,
               check=true)

Define a deterministic equivalent model for the stochastic capacity expansion
problem.

### Required arguments
`tree` is a reference to a scenario tree

`probabilities` is either a function, which returns a dictionary of the probabilities
of all nodes in a tree, or simply the dictionary itself

`sub_problem_builder` is a function mapping a node to a JuMP model for each subproblems

`solver` is a reference to the optimizer used for this problem (with appropriate settings)

### Optional arguments
`discount_factor` is a number between 0 and 1 defining a constant discount factor along each arc
in the scenario tree

`risk` is a tuple with the two CVaR parameters: (λ, α)

`sideconstraints` is a function which specifies side constraints in the master problem, see
[Tutorial 9: Side-constraints](@ref) for further details.

`check` is a boolean, which can be set to `false` to disable the validation of the JuDGE model.

### Examples
	deteq = DetEqModel(tree, ConditionallyUniformProbabilities, sub_problems,
                                    Gurobi.Optimizer)
	judge = DetEqModel(tree, probabilities, sub_problems, CPLEX.Optimizer,
                                    discount_factor=0.9, risk=(0.5,0.1)))
"""
function DetEqModel(
    tree::AbstractTree,
    probabilities,
    sub_problem_builder::Function,
    solver;
    discount_factor = 1.0,
    risk::Union{Risk,Vector{Risk}} = RiskNeutral(),
    sideconstraints = nothing,
    check = true,
    perfect_foresight = false,
)
    println("")
    println(
        "Establishing deterministic equivalent model for tree: " * string(tree),
    )
    if typeof(probabilities) <: Function
        probabilities = probabilities(tree)
    end
    if typeof(probabilities) != Dict{AbstractTree,Float64}
        error(
            "\'probabilities\' needs to be a dictionary mapping AbstractTree to Float64\nor a function that generates such a dictionary",
        )
    end

    nodes = collect(tree)

    sub_problems = Dict(i => sub_problem_builder(i) for i in nodes)

    if check
        print("Checking sub-problem format...")
        check_specification_is_legal(sub_problems)
        println("Passed")
    else
        println("Skipping checks of sub-problem format")
    end
    if typeof(sub_problem_builder) <: Function
        JuDGE.scale_objectives(tree, sub_problems, discount_factor)
    end

    if !perfect_foresight
        print("Building deterministic equivalent problem...")
        problem = build_deteq(
            sub_problems,
            tree,
            probabilities,
            solver,
            discount_factor,
            risk,
            sideconstraints,
        )
        println("Complete")
        return DetEqModel(problem, tree, probabilities, risk)
    else
        scenarios = Dict{AbstractTree,DetEqModel}()
        pr = Dict{AbstractTree,Float64}()
        @info "Building deterministic equivalent problems"
        scen_trees = get_scenarios(tree)
        for t in scen_trees
            leaf = getID(get_leafnodes(t)[1])
            sps = Dict(i => sub_problems[getID(i)] for i in collect(t))
            probs = Dict(i => 1.0 for i in collect(t))
            problem = build_deteq(
                sps,
                t,
                probs,
                solver,
                discount_factor,
                risk,
                sideconstraints,
            )
            scenarios[leaf] = DetEqModel(problem, t, probs, risk)
            pr[leaf] = probabilities[leaf]
        end
        @info "Deterministic equivalent problems built"
        return scenarios, pr
    end
end

function build_deteq(
    sub_problems::T where {T<:Dict},
    tree::T where {T<:AbstractTree},
    probabilities::Dict{AbstractTree,Float64},
    solver,
    discount_factor::Float64,
    risk::Any,
    sideconstraints,
)
    model = JuMP.Model(solver)

    @objective(model, Min, 0)

    model.ext[:vars] = Dict()
    model.ext[:master_vars] = Dict{AbstractTree,Dict{Symbol,Any}}()
    model.ext[:master_names] = Dict{AbstractTree,Dict{Symbol,Any}}()

    scen_con = Dict{Leaf,ConstraintRef}()
    scen_var = Dict{Leaf,VariableRef}()

    leafs = get_leafnodes(tree)
    for leaf in leafs
        scen_var[leaf] = @variable(model)
        JuMP.set_name(scen_var[leaf], string("scenario_obj#", leaf.name))
        scen_con[leaf] = @constraint(model, 0 == scen_var[leaf])
        set_objective_coefficient(model, scen_var[leaf], probabilities[leaf])
    end

    risk_objectives = AffExpr[]
    remain = 1.0
    if typeof(risk) == Risk
        if (risk.α == 1.0 || risk.λ == 0.0) && risk.bound == nothing
            risk = []
        else
            risk = [risk]
        end
    end
    for i in 1:length(risk)
        risk_objective = AffExpr(0.0)
        eta = @variable(model)
        for leaf in leafs
            v = @variable(model)
            w = @variable(model)
            set_lower_bound(v, 0.0)
            set_lower_bound(w, 0.0)
            if risk[i].offset == nothing
                @constraint(model, v >= eta - scen_var[leaf])
                @constraint(model, w >= scen_var[leaf] - eta)
                add_to_expression!(
                    risk_objective,
                    scen_var[leaf] * probabilities[leaf],
                )
            else
                @constraint(
                    model,
                    v >= eta - scen_var[leaf] + risk[i].offset[leaf]
                )
                @constraint(
                    model,
                    w >= scen_var[leaf] - risk[i].offset[leaf] - eta
                )
                add_to_expression!(
                    risk_objective,
                    (scen_var[leaf] - risk[i].offset[leaf]) *
                    probabilities[leaf],
                )
            end
            add_to_expression!(
                risk_objective,
                probabilities[leaf] * (v + w / risk[i].α * (1 - risk[i].α)),
            )
        end
        remain -= risk[i].λ
        push!(risk_objectives, risk_objective)
    end

    for (node, sp) in sub_problems
        leafnodes = get_leafnodes(node)
        model.ext[:vars][node] = Dict()
        for variable in all_variables(sp)
            model.ext[:vars][node][variable] =
                JuDGE.copy_variable!(model, variable)

            var_name =
                isempty(JuMP.name(variable)) ?
                string("_[", index(variable).value, "]") : JuMP.name(variable)

            JuMP.set_name(
                model.ext[:vars][node][variable],
                string(var_name, "#", node.name),
            )

            if variable == sp.ext[:objective]
                for leaf in leafnodes
                    set_normalized_coefficient(
                        scen_con[leaf],
                        model.ext[:vars][node][variable],
                        1.0,
                    )
                end
            end
        end
        copy_values = false
        loct = list_of_constraint_types(sp)

        for ct in loct
            for con in all_constraints(sp, ct[1], ct[2])
                con_obj = JuMP.constraint_object(con)
                if typeof(con_obj.func) == VariableRef
                    LHS = AffExpr(0.0)
                    add_to_expression!(
                        LHS,
                        1,
                        model.ext[:vars][node][JuMP.constraint_object(
                            con,
                        ).func],
                    )
                elseif typeof(con_obj.func) <: GenericAffExpr
                    LHS = AffExpr(0.0)
                    for (v, c) in con_obj.func.terms
                        add_to_expression!(LHS, c, model.ext[:vars][node][v])
                    end
                elseif typeof(con_obj.func) <: GenericQuadExpr
                    LHS = QuadExpr(
                        AffExpr(0.0),
                        OrderedDict{UnorderedPair{VariableRef},Float64}(),
                    )
                    for (v, c) in con_obj.func.terms
                        LHS.terms[UnorderedPair{VariableRef}(
                            model.ext[:vars][node][v.a],
                            model.ext[:vars][node][v.b],
                        )] = c
                    end
                    for (v, c) in con_obj.func.aff.terms
                        add_to_expression!(
                            LHS.aff,
                            c,
                            model.ext[:vars][node][v],
                        )
                    end
                elseif typeof(con_obj.func) ==
                       Vector{GenericAffExpr{Float64,VariableRef}}
                    group = Vector{AffExpr}()
                    for aff in con_obj.func
                        LHS = AffExpr(0.0)
                        for (v, c) in aff.terms
                            add_to_expression!(
                                LHS,
                                c,
                                model.ext[:vars][node][v],
                            )
                        end
                        push!(group, LHS)
                    end
                else
                    error(
                        "Unsupported constraint type found: " *
                        string(typeof(con_obj.func)),
                    )
                end
                set = con_obj.set
                if typeof(set) == MOI.GreaterThan{Float64}
                    @constraint(model, LHS >= set.lower)
                elseif typeof(set) == MOI.LessThan{Float64}
                    @constraint(model, LHS <= set.upper)
                elseif typeof(set) == MOI.EqualTo{Float64}
                    @constraint(model, LHS == set.value)
                elseif typeof(set) == MOI.SecondOrderCone
                    @constraint(model, group in SecondOrderCone())
                elseif typeof(set) ==
                       MOI.Indicator{MOI.ACTIVATE_ON_ZERO,MOI.EqualTo{Float64}}
                    @constraint(
                        model,
                        !collect(keys(group[1].terms))[1] =>
                            {group[2] == set.set.value}
                    )
                elseif typeof(set) ==
                       MOI.Indicator{MOI.ACTIVATE_ON_ZERO,MOI.LessThan{Float64}}
                    @constraint(
                        model,
                        !collect(keys(group[1].terms))[1] =>
                            {group[2] <= set.set.value}
                    )
                elseif typeof(set) == MOI.Indicator{
                    MOI.ACTIVATE_ON_ZERO,
                    MOI.GreaterThan{Float64},
                }
                    @constraint(
                        model,
                        !collect(keys(group[1].terms))[1] =>
                            {group[2] >= set.set.value}
                    )
                elseif typeof(set) ==
                       MOI.Indicator{MOI.ACTIVATE_ON_ONE,MOI.EqualTo{Float64}}
                    @constraint(
                        model,
                        collect(keys(group[1].terms))[1] =>
                            {group[2] == set.set.value}
                    )
                elseif typeof(set) ==
                       MOI.Indicator{MOI.ACTIVATE_ON_ONE,MOI.LessThan{Float64}}
                    @constraint(
                        model,
                        collect(keys(group[1].terms))[1] =>
                            {group[2] <= set.set.value}
                    )
                elseif typeof(set) == MOI.Indicator{
                    MOI.ACTIVATE_ON_ONE,
                    MOI.GreaterThan{Float64},
                }
                    @constraint(
                        model,
                        collect(keys(group[1].terms))[1] =>
                            {group[2] >= set.set.value}
                    )
                elseif typeof(set) != MOI.ZeroOne && typeof(set) != MOI.Integer
                    error(
                        "Unsupported constraint type found: " *
                        string(typeof(set)),
                    )
                else
                    continue
                end
            end
        end
    end

    for (node, sp) in sub_problems
        model.ext[:master_vars][node] = Dict{Symbol,Any}()
        model.ext[:master_names][node] = Dict{Symbol,Any}()
        for (name, exps) in sp.ext[:expansions]
            if isa(exps, VariableRef)
                variable = sp.ext[:expansions][name]
                model.ext[:master_vars][node][name] =
                    JuDGE.copy_variable!(model, variable)

                var_name =
                    isempty(JuMP.name(variable)) ?
                    string("_[", index(variable).value, "]") :
                    JuMP.name(variable)

                JuMP.set_name(
                    model.ext[:master_vars][node][name],
                    string(var_name, "#", node.name),
                )

                model.ext[:master_names][node][name] = string(name)
            elseif typeof(exps) <: AbstractArray
                variables = sp.ext[:expansions][name]
                model.ext[:master_vars][node][name] = Dict()
                model.ext[:master_names][node][name] = Dict()
                for index in keys(exps)
                    key = densekey_to_tuple(index)
                    model.ext[:master_vars][node][name][key] =
                        JuDGE.copy_variable!(model, variables[index])

                    var_name =
                        isempty(JuMP.name(variables[index])) ?
                        string("_[", index(variables[index]).value, "]") :
                        JuMP.name(variables[index])

                    JuMP.set_name(
                        model.ext[:master_vars][node][name][key],
                        string(var_name, "#", node.name),
                    )
                    if typeof(exps) <: JuMP.Containers.SparseAxisArray ||
                       typeof(exps) <: JuMP.Containers.DenseAxisArray
                        model.ext[:master_names][node][name][key] =
                            string(name) * "[" * string(index.I) * "]"
                    else
                        model.ext[:master_names][node][name][key] =
                            string(name) * "[" * string(index) * "]"
                    end
                end
            end
            if sp.ext[:options][name][5] != nothing
                if typeof(exps) <: AbstractArray
                    for index in keys(exps)
                        key = densekey_to_tuple(index)
                        set_lower_bound(
                            model.ext[:master_vars][node][name][key],
                            sp.ext[:options][name][5],
                        )
                    end
                else
                    set_lower_bound(
                        model.ext[:master_vars][node][name],
                        sp.ext[:options][name][5],
                    )
                end
            end
            if sp.ext[:options][name][6] != nothing
                if typeof(exps) <: AbstractArray
                    for index in keys(exps)
                        key = densekey_to_tuple(index)
                        set_upper_bound(
                            model.ext[:master_vars][node][name][key],
                            sp.ext[:options][name][6],
                        )
                    end
                else
                    set_upper_bound(
                        model.ext[:master_vars][node][name],
                        sp.ext[:options][name][6],
                    )
                end
            end
        end
    end

    for leaf in get_leafnodes(tree)
        nodes = history(leaf)
        for n in eachindex(nodes)
            node = nodes[n]
            sp = sub_problems[node]
            df = discount_factor^depth(node)
            for (name, exps) in sp.ext[:expansions]
                interval =
                    max(
                        1,
                        n - sp.ext[:options][name][3] -
                        sp.ext[:options][name][2] + 1,
                    ):n-sp.ext[:options][name][2]
                disc = Dict{Int,Float64}()
                for i in interval
                    disc[i] = df / discount_factor^(i - 1)
                end
                if isa(exps, VariableRef)
                    variable = sp.ext[:expansions][name]
                    cost_coef = df * coef(sp.ext[:capitalcosts], variable)
                    for j in interval
                        cost_coef +=
                            disc[j] * coef(sp.ext[:ongoingcosts], variable)
                    end
                    set_normalized_coefficient(
                        scen_con[leaf],
                        model.ext[:master_vars][node][name],
                        cost_coef,
                    )
                elseif typeof(exps) <: AbstractArray
                    variables = sp.ext[:expansions][name]
                    for index in keys(exps)
                        key = densekey_to_tuple(index)
                        cost_coef =
                            df * coef(sp.ext[:capitalcosts], variables[index])
                        for j in interval
                            cost_coef +=
                                disc[j] *
                                coef(sp.ext[:ongoingcosts], variables[index])
                        end
                        set_normalized_coefficient(
                            scen_con[leaf],
                            model.ext[:master_vars][node][name][key],
                            cost_coef,
                        )
                    end
                end
            end
        end
    end

    for (node, sp) in sub_problems
        past = history(node)
        leafnodes = get_leafnodes(node)
        for (name, exps) in sp.ext[:expansions]
            interval =
                sp.ext[:options][name][2]+1:min(
                    sp.ext[:options][name][2] + sp.ext[:options][name][3],
                    length(past),
                )
            if isa(exps, VariableRef)
                if sp.ext[:options][name][1] == :cumulative
                    if length(interval) == 0
                        expr = 0
                    else
                        expr = sum(
                            model.ext[:master_vars][past[index]][name] for
                            index in interval
                        )
                    end
                elseif sp.ext[:options][name][1] == :state
                    if node.parent == nothing
                        expr =
                            model.ext[:master_vars][node][name] -
                            sp.ext[:options][name][7]
                    else
                        expr =
                            model.ext[:master_vars][node][name] -
                            model.ext[:master_vars][node.parent][name]
                    end
                end
                if sp.ext[:options][name][8][1] != Inf
                    v1 = @variable(model)
                    set_lower_bound(v1, 0)
                    expr -= v1
                    for leaf in leafnodes
                        set_normalized_coefficient(
                            scen_con[leaf],
                            v1,
                            sp.ext[:options][name][8][1],
                        )
                    end
                end
                if sp.ext[:options][name][8][2] != Inf
                    v2 = @variable(model)
                    set_lower_bound(v2, 0)
                    expr += v2
                    for leaf in leafnodes
                        set_normalized_coefficient(
                            scen_con[leaf],
                            v2,
                            sp.ext[:options][name][8][2],
                        )
                    end
                end
                @constraint(model, model.ext[:vars][node][exps] == expr)
                # if typeof(node)==Leaf && sp.ext[:options][name][1]==:shutdown
                #     @constraint(model,sum(model.ext[:master_vars][n][name] for n in history_function(node))<=1)
                # end
            elseif typeof(exps) <: AbstractArray
                for i in keys(exps)
                    key = densekey_to_tuple(i)
                    if sp.ext[:options][name][1] == :cumulative
                        if length(interval) == 0
                            expr = 0
                        else
                            expr = sum(
                                model.ext[:master_vars][past[index]][name][key]
                                for index in interval
                            )
                        end
                    elseif sp.ext[:options][name][1] == :state
                        if node.parent == nothing
                            expr =
                                model.ext[:master_vars][node][name][key] -
                                sp.ext[:options][name][7]
                        else
                            expr =
                                model.ext[:master_vars][node][name][key] -
                                model.ext[:master_vars][node.parent][name][key]
                        end
                    end
                    if sp.ext[:options][name][8][1] != Inf
                        v1 = @variable(model)
                        set_lower_bound(v1, 0)
                        expr -= v1
                        for leaf in leafnodes
                            set_normalized_coefficient(
                                scen_con[leaf],
                                v1,
                                sp.ext[:options][name][8][1],
                            )
                        end
                    end
                    if sp.ext[:options][name][8][2] != Inf
                        v2 = @variable(model)
                        set_lower_bound(v2, 0)
                        expr += v2
                        for leaf in leafnodes
                            set_normalized_coefficient(
                                scen_con[leaf],
                                v2,
                                sp.ext[:options][name][8][2],
                            )
                        end
                    end
                    @constraint(model, model.ext[:vars][node][exps[i]] == expr)

                    # if typeof(node)==Leaf && sp.ext[:options][name][1]==:shutdown
                    #     @constraint(model,sum(model.ext[:master_vars][n][name][i] for n in history_function(node))<=1)
                    # end
                end
            end
        end
    end

    if typeof(sideconstraints) <: Function
        map(Main.eval, unpack_expansions(model.ext[:master_vars])) #bring expansion variables into global scope
        sideconstraints(model, tree)
        map(Main.eval, clear_expansions(model.ext[:master_vars]))
    end

    objective_fn = objective_function(model) * remain

    for i in 1:length(risk_objectives)
        objective_fn += risk_objectives[i] * risk[i].λ
        if risk[i].bound != nothing
            if risk[i].penalty != nothing
                surplus = @variable(model)
                set_lower_bound(surplus, 0)
                @constraint(
                    model,
                    risk_objectives[i] <= risk[i].bound + surplus
                )
                objective_fn += risk[i].penalty * surplus
            else
                @constraint(model, risk_objectives[i] <= risk[i].bound)
            end
        end
    end

    set_objective_function(model, objective_fn)

    model.ext[:scenario_obj] = scen_var

    return model
end

"""
	solve(deteq::DetEqModel)

Solve a determinisitc equivalent model.

### Required Arguments
`deteq` is the determinisitc equivalent model that we wish to solve.

### Example
    JuDGE.solve(deteq)
"""
function solve(deteq::DetEqModel)
    print("Solving deterministic equivalent formulation...")
    optimize!(deteq.problem)
    if termination_status(deteq.problem) == MOI.OPTIMAL
        println("Solved.")
    else
        println("Not solved: " * string(termination_status(deteq.problem)))
    end
end

function get_objval(deteq::DetEqModel; risk = deteq.risk)
    scenario_objs = Dict{Leaf,Float64}()

    for (leaf, var) in deteq.problem.ext[:scenario_obj]
        scenario_objs[leaf] = JuMP.value(var)
    end

    return compute_objval(scenario_objs, deteq.probabilities, risk)
end
