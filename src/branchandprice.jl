struct BranchConstraint
    node::Union{JuMP.Model,AbstractTree}
    expression::Union{VariableRef,AffExpr}
    relation::Symbol
    rhs::Float64

    function BranchConstraint(
        expression::Union{VariableRef,AffExpr},
        relation::Symbol,
        rhs::Float64,
        node::Union{JuMP.Model,AbstractTree},
    )
        if relation ∉ [:le, :eq, :ge]
            error("relation in BranchConstraint must be ':le', ':ge' or ':eq'")
        end
        return new(node, expression, relation, rhs)
    end
end

struct Branch
    filter::Union{Nothing,Function}
    constraints::Vector{BranchConstraint}

    function Branch(constraint::BranchConstraint)
        return new(nothing, [constraint])
    end
    function Branch(constraints::Vector{BranchConstraint})
        return new(nothing, constraints)
    end
    function Branch(constraint::BranchConstraint, filter::Function)
        return new(filter, [constraint])
    end
    function Branch(constraints::Vector{BranchConstraint}, filter::Function)
        return new(filter, constraints)
    end
end

function copy_model(
    jmodel::JuDGEModel,
    branch::Union{Nothing,Branch},
    warm_start::Bool,
)
    newmodel = JuDGEModel(
        jmodel.tree,
        jmodel.probabilities,
        jmodel.master_problem,
        jmodel.sub_problems,
        jmodel.master_solver,
        Bounds(jmodel.bounds.UB, jmodel.bounds.LB),
        jmodel.discount_factor,
        jmodel.risk,
        jmodel.sideconstraints,
        jmodel.ext,
    )

    if typeof(branch) != Nothing
        newmodel.bounds.UB = Inf
        push!(newmodel.ext[:branches], branch)
    end

    return newmodel
end

"""
	variable_branch(jmodel::JuDGEModel, inttol::Float64)

This is an in-built function that is called during branch-and-price to perform a branch.
Users can define their own functions that follow this format to create new branching strategies.

### Required Arguments
`jmodel` is the JuDGE model

`inttol` is the maximum permitted deviation from binary/integer for a value to still be considered
binary/integer feasible.
"""
function variable_branch(jmodel::JuDGEModel, inttol::Float64)
    #optimal_value = jmodel.ext[:best_integer_solution]
    function get_fractionality(var::VariableRef)
        val = JuMP.value(var)
        return min(val - floor(val), ceil(val) - val), val
        #return max(val - optimal_value[var], optimal_value[var] - val), val
    end

    function return_branch(var::VariableRef, val::Float64)
        ex = @expression(master, var)
        branch1 = BranchConstraint(ex, :ge, ceil(val), master)
        branch2 = BranchConstraint(ex, :le, floor(val), master)
        return [Branch(branch1), Branch(branch2)]
    end

    master = jmodel.master_problem
    expansions = master.ext[:expansions]
    nodes = collect(jmodel.tree, order = :breadth)

    for node in nodes
        cutoff = inttol
        best_val = 0
        v = nothing
        for x in keys(expansions[node])
            if master.ext[:options][x][4] != :Con
                var = expansions[node][x]
                if typeof(var) <: AbstractArray
                    for key in get_keys(var)
                        fractionality, val = get_fractionality(var[key])
                        if fractionality > cutoff
                            cutoff = fractionality
                            best_val = val
                            v = var[key]
                        end
                    end
                else
                    fractionality, val = get_fractionality(var)
                    if fractionality > cutoff
                        cutoff = fractionality
                        best_val = val
                        v = var
                    end
                end
            end
        end

        if v !== nothing
            return return_branch(v, best_val)
        end
    end

    slacks = master.ext[:cover_slacks]

    for node in nodes
        cutoff = inttol
        best_val = 0
        v = nothing
        for x in keys(expansions[node])
            if master.ext[:options][x][4] != :Con
                var = expansions[node][x]
                slack = slacks[node][x]
                if typeof(var) <: AbstractArray
                    for key in get_keys(slack)
                        for s in keys(slack[key])
                            fractionality, val =
                                get_fractionality(slack[key][s])
                            if fractionality > cutoff
                                cutoff = fractionality
                                best_val = val
                                v = slack[key][s]
                            end
                        end
                    end
                else
                    for s in keys(slack)
                        fractionality, val = get_fractionality(slack[s])
                        if fractionality > cutoff
                            cutoff = fractionality
                            best_val = val
                            v = slack[s]
                        end
                    end
                end
            end
        end

        if v !== nothing
            return return_branch(v, best_val)
        end
    end

    subproblems = jmodel.sub_problems

    for node in nodes
        cutoff = inttol
        index = 0
        best_val = 0
        if subproblems[node].ext[:form] == :mixed
            for i in eachindex(subproblems[node].ext[:discrete_branch])
                fractionality, val =
                    get_fractionality(master.ext[:discrete_var][node][i])
                if fractionality > cutoff
                    cutoff = fractionality
                    best_val = val
                    index = i
                end
            end

            if index != 0
                function func1(col::Column)
                    return (
                        col.node == node && col.solution[index] >= best_val
                    ) ? :ban : nothing
                end
                function func2(col::Column)
                    return (
                        col.node == node && col.solution[index] < best_val
                    ) ? :ban : nothing
                end

                ex = @expression(
                    subproblems[node],
                    subproblems[node].ext[:discrete_branch][index]
                )
                branch1 = BranchConstraint(ex, :le, ceil(best_val) - 1, node)
                branch2 = BranchConstraint(ex, :ge, ceil(best_val), node)

                return [Branch(branch1, func1), Branch(branch2, func2)]
            end
        end
    end

    return
end

# function variable_branch(jmodel::JuDGEModel, inttol::Float64)
#     master = jmodel.master_problem
#     subproblems = jmodel.sub_problems
#     tree = jmodel.tree
#     expansions = master.ext[:expansions]
#     slacks = master.ext[:cover_slacks]

#     for node in collect(tree, order = :breadth)
#         if subproblems[node].ext[:form] == :mixed
#             cutoff = inttol
#             index = 0
#             best_val = 0
#             for i in eachindex(subproblems[node].ext[:discrete_branch])
#                 val = JuMP.value(master.ext[:discrete_var][node][i])
#                 fractionality = min(val - floor(val), ceil(val) - val)
#                 if fractionality > cutoff
#                     cutoff = fractionality
#                     best_val = val
#                     index = i
#                 end
#             end

#             if index != 0
#                 function func1(col::Column)
#                     return (
#                         col.node == node && col.solution[index] >= best_val
#                     ) ? :ban : nothing
#                 end
#                 function func2(col::Column)
#                     return (
#                         col.node == node && col.solution[index] < best_val
#                     ) ? :ban : nothing
#                 end

#                 ex = @expression(
#                     subproblems[node],
#                     subproblems[node].ext[:discrete_branch][index]
#                 )
#                 branch1 = BranchConstraint(ex, :le, ceil(best_val) - 1, node)
#                 branch2 = BranchConstraint(ex, :ge, ceil(best_val), node)

#                 return [Branch(branch1, func1), Branch(branch2, func2)]
#             end
#         end
#         for x in keys(expansions[node])
#             if master.ext[:options][x][4] != :Con
#                 var = expansions[node][x]
#                 slack = slacks[node][x]

#                 if typeof(var) <: AbstractArray
#                     for key in get_keys(var)
#                         val = JuMP.value(var[key])
#                         fractionality = min(val - floor(val), ceil(val) - val)
#                         if fractionality > inttol
#                             ex = @expression(master, var[key])
#                             branch1 =
#                                 BranchConstraint(ex, :ge, ceil(val), master)
#                             branch2 =
#                                 BranchConstraint(ex, :le, floor(val), master)
#                             return [Branch(branch1), Branch(branch2)]
#                         end
#                     end

#                     for key in get_keys(var)
#                         for v in keys(slack[key])
#                             val = JuMP.value(slack[key][v])
#                             fractionality = min(val - floor(val), ceil(val) - val)
#                             if fractionality > inttol
#                                 ex = @expression(master, slack[key][v])
#                                 branch1 =
#                                     BranchConstraint(ex, :ge, ceil(val), master)
#                                 branch2 =
#                                     BranchConstraint(ex, :le, floor(val), master)
#                                 return [Branch(branch1), Branch(branch2)]
#                             end                            
#                         end
#                     end                    
#                 else
#                     val = JuMP.value(var)
#                     fractionality = min(val - floor(val), ceil(val) - val)
#                     if fractionality > inttol
#                         ex = @expression(master, var)
#                         branch1 = BranchConstraint(ex, :ge, ceil(val), master)
#                         branch2 = BranchConstraint(ex, :le, floor(val), master)
#                         return [Branch(branch1), Branch(branch2)]
#                     end

#                     for v in keys(slack)
#                         val = JuMP.value(slack[v])
#                         fractionality = min(val - floor(val), ceil(val) - val)
#                         if fractionality > inttol
#                             ex = @expression(master, slack[v])
#                             branch1 =
#                                 BranchConstraint(ex, :ge, ceil(val), master)
#                             branch2 =
#                                 BranchConstraint(ex, :le, floor(val), master)
#                             return [Branch(branch1), Branch(branch2)]
#                         end                            
#                     end                    
#                 end
#             end
#         end
#     end
#     return
# end

function perform_branch(
    jmodel::JuDGEModel,
    branches::Vector{Branch},
    warm_starts::Bool,
)
    newmodels = Vector{JuDGEModel}()

    jmodel.ext[:state] =
        jmodel.ext[:state] == :incumbent ? :incumbent_branched : :branched
    for i in eachindex(branches)
        push!(newmodels, copy_model(jmodel, branches[i], warm_starts))
        newmodels[end].ext[:path] = copy(jmodel.ext[:path])
        push!(newmodels[end].ext[:path], i)
        newmodels[end].ext[:state] = :pending
    end
    return newmodels
end

function set_banned_variables!(jmodel::JuDGEModel; remove_all = false)
    if :branches in keys(jmodel.ext)
        for node in collect(jmodel.tree)
            for col in jmodel.master_problem.ext[:columns][node]
                if haskey(node.ext, :col_max)
                    set_upper_bound(col.var, node.ext[:col_max])
                else
                    set_upper_bound(col.var, 1.0)
                end
                set_lower_bound(col.var, 0.0)
                if !remove_all
                    for branch in jmodel.ext[:branches]
                        if branch.filter !== nothing
                            if branch.filter(col) == :ban
                                set_upper_bound(col.var, 0.0)
                                break
                            end
                        end
                    end
                end
            end
        end
    end
end

function remove_branch_constraints!(jmodel::JuDGEModel)
    for con in jmodel.master_problem.ext[:branch_cons]
        delete(jmodel.master_problem, con)
    end
    jmodel.master_problem.ext[:branch_cons] = ConstraintRef[]
    return nothing
end

function add_branch_constraints!(jmodel::JuDGEModel)
    b_con = Tuple{JuMP.Model,ConstraintRef}[]

    for branch in jmodel.ext[:branches]
        for constraint in branch.constraints
            if constraint.node != jmodel.master_problem
                con = nothing
                sp = jmodel.sub_problems[constraint.node]
                if constraint.relation == :eq
                    con =
                        @constraint(sp, constraint.expression == constraint.rhs)
                elseif constraint.relation == :le
                    con =
                        @constraint(sp, constraint.expression <= constraint.rhs)
                elseif constraint.relation == :ge
                    con =
                        @constraint(sp, constraint.expression >= constraint.rhs)
                end
                push!(b_con, (sp, con))
            else
                con = nothing
                if constraint.relation == :eq
                    con = @constraint(
                        jmodel.master_problem,
                        constraint.expression == constraint.rhs
                    )
                elseif constraint.relation == :le
                    con = @constraint(
                        jmodel.master_problem,
                        constraint.expression <= constraint.rhs
                    )
                elseif constraint.relation == :ge
                    con = @constraint(
                        jmodel.master_problem,
                        constraint.expression >= constraint.rhs
                    )
                end
                push!(jmodel.master_problem.ext[:branch_cons], con)
            end
        end
    end
    return b_con
end

"""
	branch_and_price(models::Union{JuDGEModel,Vector{JuDGEModel}};
		branch_method::Function=JuDGE.variable_branch,search::Symbol=:lowestLB,
		termination::Termination=Termination(),
		max_no_int::Int=typemax(Int),
		blocks::Union{Nothing,Vector{Vector{AbstractTree}}}=nothing,
		warm_starts::Bool=false,
		optimizer_attributes::Union{Nothing,Function}=nothing,
		mp_callback::Union{Nothing,Function}=nothing,
		bp_callback::Union{Nothing,Function}=nothing,
		heuristic::Union{Nothing,Function}=nothing,
		verbose::Int=2)

Solve a JuDGEModel `judge` without branch and price.

### Required Arguments
`judge` is the JuDGE model that we wish to solve.

### Optional Arguments
`branch_method` is a function specifies the way that constraints are added to create new nodes in the
branch_and_price tree.

`search` specifies the order in which nodes are solved in the (branch-and-price) tree. Options are:
`:lowestLB`, `:depth_first_dive`, `:depth_first`, `:breadth_first`.

`termination` is a `Termination` object containing all the stopping conditions.

`max_no_int` is the maximum number of iterations yielding a fractional solution before a MIP solve is
performed on the master. By default, the MIP solves will not occur until the relaxed bound gap is less
than the `relgap` / `absgap` stopping conditions. To override this, set `max_no_int` to the negative
of the number desired value.

`blocks` specifies the groups of nodes to solve in each iteration (these groups can be generated using
`JuDGE.get_groups()`, or created manually), after all nodes have been solved, a full pricing iteration
is used to compute an updated lower bound. See `advanced.jl` for more details.

`warm_starts` boolean specifing whether to use warm starts for subproblems and binary solves of master
problem.

`optimizer_attributes` can be set to a specific function that dynamically changes optimizer attributes
for the subproblems; this should only be used by people who have examined the `advanced.jl` example.

`mp_callback` is a user-defined function that specifies termination conditions for MIP solves of the
master problem. See examples/advanced.jl.

`bp_callback` is a user-defined function that allows you to modify the `Termination` conditions for
`JuDGE.solve` and the `search` policy during the branch-and-price process.

`heuristic` is a user-defined function that typically would perform an improvement heuristic when a
new incumbent is found.

`verbose` if 0, most output from `JuDGE.solve` will be suppressed, if 1, the subproblem solve process
will be suppressed. Default is 2.

### Examples
	JuDGE.branch_and_price(jmodel, termination=Termination(abstol=10^-6))
	JuDGE.branch_and_price(jmodel, search=:depth_first_dive, verbose=0)
"""
function branch_and_price(
    models::Union{JuDGEModel,Vector{JuDGEModel}};
    branch_method::Function = JuDGE.variable_branch,
    search::Symbol = :lowestLB,
    termination::Termination = Termination(),
    max_no_int::Int = typemax(Int),
    warm_starts::Bool = false,
    blocks::Union{Nothing,Vector{Vector{AbstractTree}}} = nothing,
    verbose::Int = 2,
    optimizer_attributes::Union{Nothing,Function} = nothing,
    mp_callback::Union{Nothing,Function} = nothing,
    bp_callback::Union{Nothing,Function} = nothing,
    heuristic::Union{Nothing,Function} = nothing,
)
    initial_time = time()

    if typeof(models) == JuDGEModel
        models = [models]
    end

    models[1].master_problem.ext[:branches] = Vector{Any}()
    models[1].master_problem.ext[:log] = Vector{Vector{ConvergenceState}}()
    models[1].ext[:path] = [1]
    models[1].ext[:index] = 1
    models[1].ext[:processing] = nothing
    models[1].ext[:state] = :pending

    if termination.allow_frac == :binary_solve
        termination.allow_frac = :binary_solve_return_relaxation
    end
    UB = Inf
    LB = Inf
    otherLB = Inf
    i = 1
    bestRef = models[1]
    best = copy_model(models[1], nothing, warm_starts)
    while true
        model = models[i]

        N = length(models)
        while model.bounds.LB > UB
            if verbose > 0
                print("\n")
            end
            println(
                "Model " *
                string(model.ext[:index]) *
                " dominated. UB: " *
                string(UB) *
                ", LB:" *
                string(LB),
            )
            model.ext[:state] = :dominated
            if i == N
                break
            end
            i += 1
            model = models[i]
        end

        bestLB = 0.0
        LB = Inf
        for j in i:N
            if models[j].bounds.LB < LB
                LB = models[j].bounds.LB
                bestLB = j
            end
        end
        if otherLB < LB
            LB = otherLB
        end

        if search == :lowestLB
            model = models[bestLB]
            deleteat!(models, bestLB)
            insert!(models, i, model)
        end

        if verbose > 0
            print("\n")
        end

        if verbose > 0 || i == 1
            println(
                "Branch & Price  |   Upper Bound   Lower Bound  |  Absolute Diff   Relative Diff  |  Fractional  |      Time",
            )
        end

        displayBP(UB, LB, time() - initial_time, model.ext[:index], N)
        # println(
        #     "Model " *
        #     string(i) *
        #     " of " *
        #     string(N) *
        #     ". UB: " *
        #     string(UB) *
        #     ", LB:" *
        #     string(LB) *
        #     ", Time: " *
        #     string(Int(floor((time() - initial_time) * 1000 + 0.5)) / 1000) *
        #     "s",
        # )

        if verbose > 0
            print("\n")
        end

        model.ext[:state] = :processing
        models[1].ext[:processing] = i
        flag = solve(
            model,
            termination = termination,
            warm_starts = warm_starts,
            prune = UB,
            optimizer_attributes = optimizer_attributes,
            mp_callback = mp_callback,
            max_no_int = max_no_int,
            blocks = blocks,
            heuristic = heuristic,
            verbose = verbose,
        )
        models[1].ext[:processing] = nothing

        push!(model.master_problem.ext[:log], copy(model.log))

        if bp_callback !== nothing
            (termination, search) =
                bp_callback(termination, search, model.master_problem.ext[:log])
        end

        if model.bounds.UB < UB
            if bestRef.bounds.LB > model.bounds.UB
                bestRef.ext[:state] = :dominated
            elseif bestRef.ext[:state] == :incumbent_branched
                bestRef.ext[:state] = :branched
            else
                bestRef.ext[:state] = :dominated
            end
            model.ext[:state] = :incumbent
            UB = model.bounds.UB
            bestRef = model
            best = copy_model(model, nothing, warm_starts)
        end

        if flag == :user_interrupt ||
           (haskey(models[1].ext, :stop) && models[1].ext[:stop] == :all)
            break
        elseif flag == :dominated
            model.ext[:state] = :dominated
        end

        bestLB = 0
        LB = otherLB
        for j in i:N
            if models[j].bounds.LB < LB
                LB = models[j].bounds.LB
                bestLB = j
            end
        end

        status = termination_status(model.master_problem)
        if model.bounds.LB <= UB &&
           (
               LB + termination.abstol < UB &&
               UB - LB > termination.reltol * abs(LB)
           ) &&
           (
               status != MOI.INFEASIBLE_OR_UNBOUNDED &&
               status != MOI.INFEASIBLE &&
               status != MOI.DUAL_INFEASIBLE &&
               status != MOI.NUMERICAL_ERROR &&
               flag != :sp_infeasible
           )
            if verbose > 0
                println("\nAttempting to branch.")
            end
            branches = branch_method(model, termination.inttol)
            if branches !== nothing && length(branches) > 0
                if verbose > 0
                    println(
                        "Adding " *
                        string(length(branches)) *
                        " new nodes to B&P tree.",
                    )
                end

                newmodels = perform_branch(model, branches, warm_starts)
                for ii in eachindex(newmodels)
                    newmodels[ii].ext[:index] = length(models) + ii
                    newmodels[ii].ext[:edge] =
                        [model.ext[:index], newmodels[ii].ext[:index]]
                end
                i += 1
                if search == :depth_first_dive
                    for j in 1:length(newmodels)
                        insert!(models, i, newmodels[length(newmodels)+1-j])
                    end
                elseif search == :breadth_first || search == :lowestLB
                    append!(models, newmodels)
                elseif search ∈ [:depth_first, :depth_first_resurface]
                    insert!(models, i, newmodels[1])
                    deleteat!(newmodels, 1)
                    append!(models, newmodels)
                end
            else
                if fractionalcount(model, termination.inttol) > 0
                    error(
                        "Solution is fractional, but branching method returned no branches.",
                    )
                end

                if i == length(models)
                    break
                end
                if model.bounds.LB <= UB
                    model.ext[:state] = :integer
                end
                i += 1
                # if flag == :sp_infeasible # cannot enter this?
                #     i += 1
                # end
                #     break
                # else
                # if model.bounds.LB < otherLB
                #     otherLB = model.bounds.LB
                # end

            end
        else
            if (
                status != MOI.INFEASIBLE_OR_UNBOUNDED &&
                status != MOI.INFEASIBLE &&
                status != MOI.DUAL_INFEASIBLE &&
                status != MOI.NUMERICAL_ERROR &&
                flag != :sp_infeasible
            ) && model.bounds.LB < otherLB
                otherLB = model.bounds.LB
            end

            i += 1

            if UB - LB <= termination.abstol ||
               (UB - LB <= termination.reltol * abs(UB) && UB != Inf)
                bestRef.ext[:state] = :converged
                break
            elseif i - 1 == length(models)
                break
            end
        end
    end
    bestLB = 0
    LB = otherLB
    for j in i:length(models)
        if models[j].bounds.LB < LB
            LB = models[j].bounds.LB
            bestLB = j
        end
        if models[j].bounds.LB > UB
            models[j].ext[:state] = :dominated
        end
    end

    if termination.allow_frac ∉ [:no_binary_solve, :first_fractional]
        if verbose == 2
            print("Performing final MIP solve")
        end
        remove_branch_constraints!(best)
        set_banned_variables!(best; remove_all = true)
        solve_binary(
            best,
            termination.abstol,
            termination.reltol,
            warm_starts,
            nothing,
        )

        if verbose == 2
            overprint("")
        end
    else
        optimize!(best.master_problem)
    end
    UB = best.bounds.UB
    models[end].ext[:finished] = true
    println(
        "\nObjective value of best integer-feasible solution: " * string(UB),
    )
    println("Objective value of lower bound: " * string(min(LB, UB)))
    println(
        "Solve time: " *
        string(Int(floor((time() - initial_time) * 1000 + 0.5)) / 1000) *
        "s",
    )
    return best
end
