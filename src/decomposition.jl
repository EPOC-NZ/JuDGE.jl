mutable struct Bounds
    UB::Float64
    LB::Float64
end

struct Column
    node::AbstractTree
    coeffs::Dict{Symbol,Any}
    obj::Float64
    var::VariableRef
    solution::Vector{Float64}
end

struct JuDGEModel
    tree::AbstractTree
    master_problem::JuMP.Model
    sub_problems::Dict{AbstractTree,JuMP.Model}
    bounds::Bounds
    discount_factor::Float64
    master_solver::Any
    probabilities::Dict{AbstractTree,Float64}
    risk::Union{Risk,Vector{Risk}}
    sideconstraints::Union{Nothing,Function}
    log::Vector{ConvergenceState}
    ext::Dict{Symbol,Any}
end

# pretty printing
function Base.show(io::IO, ::MIME"text/plain", judge::JuDGEModel)
    print(io, "JuDGE Model with:\n")
    println(io, "  Tree: ", judge.tree)
    print(io, "  Expansion variables: ")
    keys = judge.sub_problems[judge.tree].ext[:expansions]
    for (i, ii) in keys
        print(io, "$(i) ")
    end
end

function Base.show(io::IO, judge::JuDGEModel)
    print(io, "JuDGE Model with:\n")
    println(io, "  Tree: ", judge.tree)
    print(io, "  Expansion variables: ")
    keys = judge.sub_problems[judge.tree].ext[:expansions]
    for (i, ii) in keys
        print(io, "$(i) ")
    end
end

function JuDGEModel(
    tree::T where {T<:AbstractTree},
    probabilities::Dict{AbstractTree,Float64},
    master_problem::JuMP.Model,
    sub_problems::Dict{AbstractTree,JuMP.Model},
    solver,
    bounds::Bounds,
    discount_factor::Float64,
    risk::Union{Risk,Vector{Risk}},
    sideconstraints::Union{Function,Nothing},
    ext::Dict{Symbol,Any},
)
    new_ext = Dict{Symbol,Any}()
    new_ext[:branches] = copy(ext[:branches])
    new_ext[:optimizer_settings] = deepcopy(ext[:optimizer_settings])
    return JuDGEModel(
        tree,
        master_problem,
        sub_problems,
        bounds,
        discount_factor,
        solver,
        probabilities,
        risk,
        sideconstraints,
        Vector{ConvergenceState}(),
        new_ext,
    )
end

"""
	JuDGEModel(tree::AbstractTree,
               probabilities,
               sub_problem_builder::Function,
               solver;
               discount_factor=1.0,
               risk=RiskNeutral(),
               sideconstraints=nothing,
               check=true,
               perfect_foresight=false)

Define a JuDGE model.

### Required arguments
`tree` is a reference to a scenario tree

`probabilities` is either a function, which returns a dictionary of the probabilities
of all nodes in a tree, or simply the dictionary itself

`sub_problem_builder` is a function mapping a node to a JuMP model for each subproblems

`solver` is a reference to the optimizer used for the master problem (with appropriate settings);
 this can also be a tuple containing two optimizers (one for solving the relaxation, and one for
 solving the binary model)

### Optional arguments
`discount_factor` is a number between 0 and 1 defining a constant discount factor along each arc
in the scenario tree

`risk` can be either a `Risk` object, or a vector of such objects.

`sideconstraints` is a function which specifies side constraints in the master problem, see
[Tutorial 9: Side-constraints](@ref) for further details

`check` is a boolean, which can be set to `false` to disable the validation of the JuDGE model.

`perfect_foresight` is a boolean; this is an experimental feature, which creates an array of
JuDGE models, one for each leaf node. This will enable users to easily compute the EVPI for
the stochastic program. Also can be used for regret-based risk implementations. See the example
EVPI_and_VSS.jl.

### Examples
	judge = JuDGEModel(tree, ConditionallyUniformProbabilities, sub_problems,
                                    Gurobi.Optimizer)
	judge = JuDGEModel(tree, probabilities, sub_problems, CPLEX.Optimizer,
                                    discount_factor=0.9, risk=Risk(0.5,0.1)))
"""
function JuDGEModel(
    tree::T where {T<:AbstractTree},
    probabilities,
    sub_problem_builder::Function,
    solver;
    discount_factor::Float64 = 1.0,
    risk::Union{Risk,Vector{Risk}} = RiskNeutral(),
    sideconstraints::Union{Function,Nothing} = nothing,
    check::Bool = true,
    perfect_foresight::Bool = false,
)
    println("")
    if !perfect_foresight
        @info "Establishing JuDGE model for tree: " * string(tree)
    else
        @info "Establishing perfect foresight models for tree: " * string(tree)
    end

    if typeof(probabilities) <: Function
        probabilities = probabilities(tree)
    end
    if typeof(probabilities) != Dict{AbstractTree,Float64}
        error(
            "\'probabilities\' needs to be a dictionary mapping AbstractTree to Float64\nor a function that generates such a dictionary",
        )
    end

    nodes = collect(tree)

    sub_problems = Dict{AbstractTree,JuMP.Model}()
    @info("Building JuMP Model for node ")
    print("\e[F")
    for i in nodes
        overprint("$(i.name)", hpos = 38)
        sub_problems[i] = sub_problem_builder(getID(i))
    end
    overprint("s...Complete\n", hpos = 28)

    if check
        @info "Checking sub-problem format..."
        print("\e[F")
        check_specification_is_legal(sub_problems)
        overprint("Complete\n", hpos = 39)
    else
        @info "Skipping checks of sub-problem format"
    end

    scale_objectives(tree, sub_problems, discount_factor)

    for sp in values(sub_problems)
        sp.ext[:all_vars] = copy(sp.obj_dict)
    end

    if !perfect_foresight
        @info "Building master problem..."
        print("\e[F")
        master_problem = build_master(
            sub_problems,
            tree,
            probabilities,
            solver,
            discount_factor,
            risk,
            sideconstraints,
        )
        overprint("Complete\n", hpos = 35)
        ext = Dict{Symbol,Any}()
        ext[:branches] = Branch[]
        ext[:optimizer_settings] = Dict{Symbol,Any}()

        return JuDGEModel(
            tree,
            master_problem,
            sub_problems,
            Bounds(Inf, -Inf),
            discount_factor,
            solver,
            probabilities,
            risk,
            sideconstraints,
            Vector{ConvergenceState}(),
            ext,
        )
    else
        scenarios = Dict{AbstractTree,JuDGEModel}()
        pr = Dict{AbstractTree,Float64}()
        @info "Building master problem for scenario given by node "
        print("\e[F")
        scen_trees = get_scenarios(tree)
        for t in scen_trees
            leaf = getID(get_leafnodes(t)[1])
            overprint("$(leaf.name)", hpos = 60)
            sps = Dict(i => sub_problems[getID(i)] for i in collect(t))
            probs = Dict(i => 1.0 for i in collect(t))
            master_problem = build_master(
                sps,
                t,
                probs,
                solver,
                discount_factor,
                risk,
                sideconstraints,
            )
            ext = Dict{Symbol,Any}()
            ext[:branches] = Branch[]
            ext[:optimizer_settings] = Dict{Symbol,Any}()
            scenarios[leaf] = JuDGEModel(
                t,
                master_problem,
                sps,
                Bounds(Inf, -Inf),
                discount_factor,
                solver,
                probs,
                risk,
                sideconstraints,
                Vector{ConvergenceState}(),
                ext,
            )
            pr[leaf] = probabilities[leaf]
        end
        overprint("s...Complete\n", hpos = 32)
        return scenarios, pr
    end
end

function add_variable_as_column(master, column)
    for constr in master.ext[:convexcombination][column.node]
        set_normalized_coefficient(constr, column.var, 1.0)
    end

    for (var, val) in column.coeffs
        if typeof(val) == Float64
            set_normalized_coefficient(
                master.ext[:coverconstraint][column.node][var],
                column.var,
                val,
            )
        else
            for (i, val2) in val
                set_normalized_coefficient(
                    master.ext[:coverconstraint][column.node][var][i],
                    column.var,
                    val2,
                )
            end
        end
    end

    for node in collect(column.node)
        if typeof(node) == Leaf
            set_normalized_coefficient(
                master.ext[:scenprofit_con][node],
                column.var,
                column.obj,
            )
        end
    end
    return nothing
end

# constraint that's used to recreate a discrete subproblem variable in the master problem
function add_mixed_cover(master, sp, column)
    if !(column.node in keys(master.ext[:discrete_con]))
        master.ext[:discrete_var][column.node] = Dict{Int,VariableRef}()
        master.ext[:discrete_con][column.node] = Dict{Int,ConstraintRef}()
        for i in 1:length(sp.ext[:discrete])
            master.ext[:discrete_var][column.node][i] = @variable(master)
            master.ext[:discrete_con][column.node][i] = @constraint(
                master,
                0.0 == master.ext[:discrete_var][column.node][i]
            )
        end
    end

    for i in 1:length(sp.ext[:discrete])
        set_normalized_coefficient(
            master.ext[:discrete_con][column.node][i],
            column.var,
            column.solution[i],
        )
    end
    return nothing
end

function add_column(
    master::JuMP.Model,
    sub_problem::JuMP.Model,
    node::AbstractTree;
    branches = nothing,
)
    coeffs = Dict{Symbol,Any}()

    sol = [0.0]

    if sub_problem.ext[:form] == :mixed
        sol = JuMP.value.(sub_problem.ext[:discrete])
    end

    for (name, variable) in sub_problem.ext[:expansions]
        if typeof(variable) <: AbstractArray
            coeffs[name] = Dict{Any,Float64}()
            vals = JuMP.value.(variable)
            for i in get_keys(variable)
                if sub_problem.ext[:options][name][4] != :Con
                    rval = round(vals[i])
                    if rval != 0.0
                        coeffs[name][i] = rval
                    end
                else
                    coeffs[name][i] = vals[i]
                end
            end
        else
            if sub_problem.ext[:options][name][4] != :Con
                rval = round(JuMP.value(variable))
                if rval != 0.0
                    coeffs[name] = rval
                end
            else
                coeffs[name] = JuMP.value(variable)
            end
        end
    end

    UB = haskey(node.ext, :col_max) ? node.ext[:col_max] : 1.0

    column = Column(
        node,
        coeffs,
        JuMP.value(sub_problem.ext[:objective]),
        JuMP.add_variable(
            master,
            JuMP.build_variable(error, UnitIntervalInformation(UB = UB)),
        ),
        sol,
    )

    if branches !== nothing
        for b in branches
            if b.filter !== nothing && b.filter(column) == :ban
                set_upper_bound(column.var, 0.0)
                println("New column banned on current branch")
                break
            end
        end
    end

    add_variable_as_column(master, column)
    push!(master.ext[:columns][node], column)
    if sub_problem.ext[:form] == :mixed
        add_mixed_cover(master, sub_problem, column)
    end

    JuMP.set_name(
        column.var,
        string("col$(length(master.ext[:columns][node]))", "#", node.name),
    )
    return column
end

"""
	solve(judge::JuDGEModel;
	      termination::Termination=Termination(),
	      max_no_int::Int=typemax(Int),
	      blocks::Union{Nothing,Vector{Vector{AbstractTree}}}=nothing,
	      warm_starts::Bool=false,
	      optimizer_attributes::Union{Nothing,Function}=nothing,
	      mp_callback::Union{Nothing,Function}=nothing,
	      prune::Float64=Inf,
	      heuristic::Union{Nothing,Function}=nothing,
	      verbose::Int=2)

Solve a JuDGEModel `judge` without branch-and-price.

### Required Arguments
`judge` is the JuDGE model that we wish to solve.

### Optional Arguments

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

`prune` is used to stop the algorithm before convergence, if a known upper bound for the problem is
specified.

`heuristic` is a user-defined function that typically would perform an improvement heuristic when a
new incumbent is found.

`skip_nodes` is the maximum number of nodes to skip in any iteration. Nodes are added to the `skip_list`
if they return a 0 reduced cost. Once the `length(skip_list)==skip_nodes` the list is emptied, and a
complete iteration is run in order to update the lower bound.

`verbose` if 0, all output from solve will be suppressed, if 1, the subproblem solve process will be
suppressed. Default is 2.

### Examples
    JuDGE.solve(jmodel, termination=Termination(rlx_abstol=10^-6))
	JuDGE.solve(jmodel, termination=Termination(rlx_abstol=10^-6), max_no_int=-5)
"""
function solve(
    judge::JuDGEModel;
    termination::Termination = Termination(),
    max_no_int::Int = typemax(Int),
    blocks::Union{Nothing,Vector{Vector{AbstractTree}}} = nothing,
    warm_starts::Bool = false,
    optimizer_attributes::Union{Nothing,Function} = nothing,
    mp_callback::Union{Nothing,Function} = nothing,
    prune::Float64 = Inf,
    heuristic::Union{Nothing,Function} = nothing,
    skip_nodes::Int = 0,
    verbose::Int = 2,
)
    current = InitialConvergenceState()
    empty!(judge.log)
    push!(judge.log, current)
    if verbose > 0 && prune == Inf
        @info "Solving JuDGE model for tree: " * string(judge.tree)
        if verbose == 2
            display(termination)
        end
        println("")
        println(
            "Relaxed ObjVal  |   Upper Bound   Lower Bound  |  Absolute Diff   Relative Diff  |  Fractional  |      Time     Iter",
        )
    end
    # set up times for use in convergence
    initial_time = time()
    obj = Inf
    heuristicobj = Inf
    nodes = collect(judge.tree)
    exit_flag = nothing
    no_int_count = 0
    if optimizer_attributes !== nothing
        optimizer_attributes(judge, false, true)
    end

    skip_list = []
    skip_counter = skip_nodes

    max_char = length(nodes[end].name) + length(string(length(nodes)))
    function get_whitespace(name::String, number::Int)
        blank = "  "
        spaces = max_char - length(name) - length(string(number))
        for i in 1:spaces
            blank *= " "
        end
        return blank
    end

    function create_column(node::AbstractTree)
        column =
            add_column(judge.master_problem, judge.sub_problems[node], node)
        if warm_starts
            set_start_value(column.var, 0.0)
        end
        return nothing
    end

    if blocks === nothing
        blocks = [nodes]
    end

    if judge.master_problem.ext[:mip]
        remove_binary(judge)
    end
    set_banned_variables!(judge)

    remove_branch_constraints!(judge)
    b_con = add_branch_constraints!(judge)

    optimize!(judge.master_problem)
    status = termination_status(judge.master_problem)

    block = length(blocks) == 1 ? 0 : 1

    if status == MOI.NUMERICAL_ERROR
        println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
        @goto terminate
    elseif status ∈
           [MOI.INFEASIBLE_OR_UNBOUNDED, MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE]
        block = 0
    end

    while true
        if block <= 0
            nodes2 = copy(nodes)
        else
            nodes2 = copy(blocks[block])
        end
        if skip_counter == 0 || length(skip_list) == length(nodes)
            empty!(skip_list)
            skip_counter = skip_nodes
        else
            skip_counter -= 1
            for n in skip_list
                index = findfirst(x -> x == n, nodes2)
                if index !== nothing
                    deleteat!(nodes2, index)
                end
            end
        end
        # perform the main iterations
        for i in eachindex(nodes2)
            node = nodes2[i]
            sp = judge.sub_problems[node]
            updateduals(judge.master_problem, sp, node, status, current.iter)
            if verbose == 2
                overprint(
                    "Solving subproblem for node " *
                    node.name *
                    get_whitespace(node.name, i) *
                    string(i) *
                    "/" *
                    string(length(nodes2)),
                )
            end

            optimize!(sp)

            if termination_status(sp) ∉
               [MOI.OPTIMAL, MOI.INTERRUPTED, MOI.TIME_LIMIT]
                if verbose == 2
                    println("")
                end
                @warn(
                    "Solve for subproblem " *
                    node.name *
                    " exited with status " *
                    string(termination_status(sp))
                )
                exit_flag = :sp_infeasible
                @goto terminate
            end
            if warm_starts
                vars = all_variables(sp)
                set_start_value.(vars, JuMP.value.(vars))
            end
        end
        if verbose == 2
            overprint("")
        end

        frac = 0

        if status ∉
           [MOI.INFEASIBLE_OR_UNBOUNDED, MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE]
            if status != MOI.OPTIMAL
                @warn(
                    "Master problem did not solve to optimality: " *
                    string(status)
                )
            elseif block == 0 && length(skip_list) == 0
                getlowerbound(judge)
            end
        elseif judge.bounds.LB > -Inf
            println("\nMaster problem is infeasible or unbounded")
            judge.bounds.LB = Inf
            judge.bounds.UB = Inf
            exit_flag = :master_infeasible
            break
        elseif current.iter > 2
            println("\nMaster problem is infeasible or unbounded")
            exit_flag = :master_infeasible
            break
        end

        num_var = num_variables(judge.master_problem)
        if status ∈
           [MOI.INFEASIBLE_OR_UNBOUNDED, MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE]
            for node in nodes2
                create_column(node)
            end
        else
            for node in nodes2
                if objective_value(judge.sub_problems[node]) < -10^-10
                    create_column(node)
                elseif skip_nodes != 0
                    push!(skip_list, node)
                end
            end
        end
        if verbose == 2
            overprint("Solving master problem")
            optimize!(judge.master_problem)
            overprint("")
        else
            optimize!(judge.master_problem)
        end
        status = termination_status(judge.master_problem)
        if status == MOI.NUMERICAL_ERROR
            println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
            exit_flag = :master_numerical_error
            break
        elseif status ∉ [
            MOI.INFEASIBLE_OR_UNBOUNDED,
            MOI.INFEASIBLE,
            MOI.DUAL_INFEASIBLE,
        ]
            frac = fractionalcount(judge, termination.inttol)
            obj = objective_value(judge.master_problem)
            if frac == 0
                if heuristic !== nothing && heuristicobj > obj
                    if heuristic(judge) < 0.0 &&
                       termination.allow_frac ∉
                       [:first_fractional, :no_binary_solve]
                        solve_master_binary(
                            judge,
                            initial_time,
                            termination,
                            warm_starts,
                            nothing,
                            verbose,
                        )
                    else
                        optimize!(judge.master_problem)
                    end
                    obj = objective_value(judge.master_problem)
                    heuristicobj = obj
                end
                judge.bounds.UB = obj

                no_int_count = 0
            else
                no_int_count += 1
            end
        else
            block = -1
        end
        current = ConvergenceState(
            obj,
            judge.bounds.UB,
            judge.bounds.LB,
            time() - initial_time,
            current.iter + 1,
            frac,
        )
        if verbose > 0
            display(current)
        end
        push!(judge.log, current)
        if prune < judge.bounds.LB
            if verbose > 0
                println("\nDominated by incumbent.")
            end
            break
        elseif has_converged(termination, current)
            if frac > 0
                solve_master_binary(
                    judge,
                    initial_time,
                    termination,
                    warm_starts,
                    nothing,
                    verbose,
                )
            end
            if heuristic !== nothing
                if heuristic(judge) < 0.0 &&
                   termination.allow_frac ∉
                   [:first_fractional, :no_binary_solve]
                    solve_master_binary(
                        judge,
                        initial_time,
                        termination,
                        warm_starts,
                        nothing,
                        verbose,
                    )
                    heuristicobj = judge.bounds.UB
                else
                    optimize!(judge.master_problem)
                end
            end
            if verbose > 0
                println("\nConvergence criteria met.")
            end
            break
        elseif termination.allow_frac == :first_fractional && frac > 0
            if verbose > 0
                println("\nFractional solution found.")
            end
            break
        elseif (
            (
                max_no_int > 0 &&
                no_int_count >= max_no_int &&
                current.num_frac > 0
            ) && (
                current.rlx_abs < termination.abstol ||
                current.rlx_rel < termination.reltol
            )
        ) || (
            max_no_int < 0 &&
            no_int_count >= -max_no_int &&
            current.num_frac > 0
        )
            current = solve_master_binary(
                judge,
                initial_time,
                termination,
                warm_starts,
                mp_callback,
                verbose,
            )
            if heuristic !== nothing && heuristicobj > judge.bounds.UB
                if heuristic(judge) < 0.0
                    current = solve_master_binary(
                        judge,
                        initial_time,
                        termination,
                        warm_starts,
                        mp_callback,
                        verbose,
                    )
                else
                    optimize!(judge.master_problem)
                end
                heuristicobj = judge.bounds.UB
            end

            if has_converged(termination, current)
                if verbose > 0
                    println("\nConvergence criteria met.")
                end
                break
            elseif termination.allow_frac == :binary_solve
                remove_binary(judge)
                optimize!(judge.master_problem)
            end
            no_int_count = 0
        end

        if optimizer_attributes === nothing
            if block == 0 && num_var == num_variables(judge.master_problem)
                if length(nodes2) == length(nodes)
                    solve_master_binary(
                        judge,
                        initial_time,
                        termination,
                        warm_starts,
                        nothing,
                        verbose,
                    )
                    if verbose > 0
                        println("\nStalled: exiting.")
                    end
                    break
                else
                    skip_counter = 0
                end
            end
        elseif optimizer_attributes(
            judge,
            num_var == num_variables(judge.master_problem),
            false,
        )
            if length(nodes2) == length(nodes)
                if verbose > 0
                    println("\nStalled: exiting.")
                end
                break
            else
                skip_counter = 0
            end
        end
        if length(blocks) == 1
            block = 0
        else
            block = (block + 1) % (length(blocks) + 1)
        end
    end
    @label terminate
    for (sp, con) in b_con
        delete(sp, con)
    end
    return exit_flag
end

function getlowerbound(judge::JuDGEModel)
    lb = objective_value(judge.master_problem)

    for (node, sp) in judge.sub_problems
        if :sum_max in keys(node.ext)
            lb += objective_bound(sp) * node.ext[:sum_max]
        else
            lb += objective_bound(sp)
        end
    end
    if lb > judge.bounds.LB
        judge.bounds.LB = lb
    end
end

function fractionalcount(jmodel::JuDGEModel, inttol::Float64)
    count = 0

    for (node, sp) in jmodel.sub_problems
        if sp.ext[:form] == :binarycolumns
            for col in jmodel.master_problem.ext[:columns][node]
                val = JuMP.value(col.var)
                count += min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
            end
        end
    end

    for node in collect(jmodel.tree)
        for x in keys(jmodel.master_problem.ext[:expansions][node])
            if jmodel.sub_problems[jmodel.tree].ext[:options][x][4] != :Con
                var = jmodel.master_problem.ext[:expansions][node][x]
                slacks = jmodel.master_problem.ext[:cover_slacks][node][x]
                if typeof(var) <: AbstractArray
                    for key in get_keys(var)
                        val = JuMP.value(var[key])
                        count +=
                            min(val - floor(val), ceil(val) - val) > inttol ?
                            1 : 0
                        for v in keys(slacks[key])
                            val = JuMP.value(slacks[key][v])
                            count +=
                                min(val - floor(val), ceil(val) - val) >
                                inttol ? 1 : 0
                        end
                    end
                else
                    val = JuMP.value(var)
                    count +=
                        min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
                    for v in keys(slacks)
                        val = JuMP.value(slacks[v])
                        count +=
                            min(val - floor(val), ceil(val) - val) > inttol ?
                            1 : 0
                    end
                end
            end
        end
        if jmodel.sub_problems[node].ext[:form] == :mixed
            for (i, var) in jmodel.master_problem.ext[:discrete_var][node]
                val = JuMP.value(var)
                count += min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
            end
        end
    end

    return count
end

function updateduals(master, sub_problem, node, status, iter)
    if status != MOI.INFEASIBLE_OR_UNBOUNDED &&
       status != MOI.INFEASIBLE &&
       status != MOI.DUAL_INFEASIBLE
        for (name, var) in sub_problem.ext[:expansions]
            if typeof(var) <: AbstractArray
                for i in get_keys(var)
                    set_objective_coefficient(
                        sub_problem,
                        var[i],
                        -dual(master.ext[:coverconstraint][node][name][i]),
                    )
                end
            else
                set_objective_coefficient(
                    sub_problem,
                    var,
                    -dual(master.ext[:coverconstraint][node][name]),
                )
            end
        end
        nodes = collect(node)
        total = 0.0
        for n in nodes
            if typeof(n) == Leaf
                total -= dual(master.ext[:scenprofit_con][n])
            end
        end
        set_objective_coefficient(
            sub_problem,
            sub_problem.ext[:objective],
            total,
        )

        cc_sum = 0.0
        for constr in master.ext[:convexcombination][node]
            set = typeof(constraint_object(constr).set)
            if set <: MOI.LessThan
                cc_sum += dual(constr)
            elseif set <: MOI.GreaterThan
                cc_sum -= dual(constr)
            elseif set <: MOI.EqualTo
                cc_sum += dual(constr)
            end
        end

        # for col in master.ext[:columns][node]
        #     cc_sum += dual(LowerBoundRef(col.var))# - dual(UpperBoundRef(col.var))
        # end

        set_objective_function(
            sub_problem,
            objective_function(sub_problem) -
            objective_function(sub_problem).constant - cc_sum,
        )
    else
        if iter % 2 == 0
            oc = -10^14
        else
            oc = 10^14
        end
        for (name, var) in sub_problem.ext[:expansions]
            if sub_problem.ext[:options][name][8][1] == Inf
                flip = -1
            else
                flip = 1
            end
            if typeof(var) <: AbstractArray
                for i in get_keys(var)
                    set_objective_coefficient(sub_problem, var[i], flip * oc)
                end
            else
                set_objective_coefficient(sub_problem, var, flip * oc)
            end
        end
        set_objective_coefficient(sub_problem, sub_problem.ext[:objective], 1.0)
    end
end

function get_objval(jmodel::JuDGEModel; risk = jmodel.risk)
    scenario_objs = Dict{Leaf,Float64}()

    for (leaf, var) in jmodel.master_problem.ext[:scenprofit_var]
        scenario_objs[leaf] = JuMP.value(var)
    end

    return compute_objval(scenario_objs, jmodel.probabilities, risk)
end

"""
	resolve_subproblems(judge::JuDGEModel)

Once a JuDGE model has converged, it is necessary to re-solve the subproblems to find the optimal decisions within each node.

### Required Arguments
`jmodel` is the JuDGE model that we wish to solve.

### Optional Arguments
`force_match` if true will force all expansion and shutdown variables in the subproblem to exactly match the master problem.

### Examples
    resolve_subproblems(judge)
"""
function resolve_subproblems(jmodel::JuDGEModel; force_match::Bool = false)
    fix_expansions(jmodel, force_match)
    return resolve_fixed(jmodel)
end

function fix_expansions(jmodel::JuDGEModel, force_match::Bool)
    if termination_status(jmodel.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel.master_problem) != MOI.LOCALLY_SOLVED &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED
        error("You need to first solve the decomposed model.")
    end

    for node in collect(jmodel.tree)
        sp = jmodel.sub_problems[node]
        set_objective_function(
            sp,
            objective_function(sp) - objective_function(sp).constant,
        )
        set_objective_coefficient(sp, sp.ext[:objective], 1.0)
        for (name, var) in jmodel.master_problem.ext[:expansions][node]
            var2 = sp.ext[:expansions][name]
            slacks = jmodel.master_problem.ext[:cover_slacks][node][name]
            if typeof(var) <: AbstractArray
                for i in get_keys(var)
                    value = 0.0
                    if sp.ext[:options][name][1] == :state
                        prev = node.parent
                        if prev === nothing
                            var3 =
                                jmodel.master_problem.ext[:expansions][node][name][i]
                            value = JuMP.value(var3) - sp.ext[:options][name][7]
                        else
                            var3 =
                                jmodel.master_problem.ext[:expansions][node][name][i]
                            var4 =
                                jmodel.master_problem.ext[:expansions][prev][name][i]
                            value = JuMP.value(var3) - JuMP.value(var4)
                        end
                    else
                        con_obj = constraint_object(
                            jmodel.master_problem.ext[:coverconstraint][node][name][i],
                        )
                        for prev in history(node)
                            var3 =
                                jmodel.master_problem.ext[:expansions][prev][name][i]
                            if var3 in keys(con_obj.func.terms)
                                value +=
                                    JuMP.value(var3) * -con_obj.func.terms[var3]
                            end
                        end
                    end

                    if haskey(slacks[i], 1)
                        value -= JuMP.value(slacks[i][1])
                    end
                    if haskey(slacks[i], 2)
                        value += JuMP.value(slacks[i][2])
                    end

                    if sp.ext[:options][name][4] ∈ [:Bin, :Int]
                        # if sp.ext[:options][name][8][2] == 0
                        #     value=floor(value)
                        # elseif sp.ext[:options][name][8][1] == 0
                        #     value=ceil(value)
                        # else
                        value = round(value)
                        # end
                    end

                    JuMP.fix(var2[i], value, force = true)
                    # if sp.ext[:options][name][8][2] == 0
                    #     if force_match
                    #         JuMP.fix(var2[i], value, force = true)
                    #     else
                    #         JuMP.set_lower_bound(var2[i], value)
                    #     end
                    # elseif sp.ext[:options][name][8][1] == 0
                    #     if force_match
                    #         JuMP.fix(var2[i], value, force = true)
                    #     else
                    #         JuMP.set_upper_bound(var2[i], value)
                    #     end
                    # elseif sp.ext[:options][name][8][1] == Inf && sp.ext[:options][name][8][2] == Inf
                    #     JuMP.fix(var2[i], value, force = true)
                    # end
                    set_objective_coefficient(sp, var2[i], 0.0)
                end
            elseif isa(var, VariableRef)
                value = 0.0
                if sp.ext[:options][name][1] == :state
                    prev = node.parent
                    if prev === nothing
                        var3 =
                            jmodel.master_problem.ext[:expansions][node][name]
                        value = JuMP.value(var3) - sp.ext[:options][name][7]
                    else
                        var3 =
                            jmodel.master_problem.ext[:expansions][node][name]
                        var4 =
                            jmodel.master_problem.ext[:expansions][prev][name]
                        value = JuMP.value(var3) - JuMP.value(var4)
                    end
                else
                    con_obj = constraint_object(
                        jmodel.master_problem.ext[:coverconstraint][node][name],
                    )
                    for prev in history(node)
                        var3 =
                            jmodel.master_problem.ext[:expansions][prev][name]
                        if var3 in keys(con_obj.func.terms)
                            value +=
                                JuMP.value(var3) * -con_obj.func.terms[var3]
                        end
                    end
                end

                if haskey(slacks, 1)
                    value -= JuMP.value(slacks[1])
                end
                if haskey(slacks, 2)
                    value += JuMP.value(slacks[2])
                end

                if sp.ext[:options][name][4] ∈ [:Bin, :Int]
                    value = round(value)
                end

                JuMP.fix(var2, value, force = true)

                # if sp.ext[:options][name][8][2] == 0
                #     if force_match
                #         JuMP.fix(var2, value, force = true)
                #     else
                #         JuMP.set_lower_bound(var2, value)
                #     end
                # elseif sp.ext[:options][name][8][1] == 0
                #     if force_match
                #         JuMP.fix(var2, value, force = true)
                #     else
                #         JuMP.set_upper_bound(var2, value)
                #     end
                # elseif sp.ext[:options][name][8][1] == Inf && sp.ext[:options][name][8][2] == Inf
                #     JuMP.fix(var2, value, force = true)
                # end
                set_objective_coefficient(sp, var2, 0.0)
            end
        end
    end
    return nothing
end

function unfix_expansions(jmodel::JuDGEModel)
    if termination_status(jmodel.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel.master_problem) != MOI.LOCALLY_SOLVED &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED
        error("You need to first solve the decomposed model.")
    end

    for node in collect(jmodel.tree)
        sp = jmodel.sub_problems[node]
        for (name, var) in jmodel.master_problem.ext[:expansions][node]
            var2 = sp.ext[:expansions][name]
            slacks = jmodel.master_problem.ext[:cover_slacks][node][name]
            if typeof(var) <: AbstractArray
                for i in get_keys(var)
                    JuMP.unfix(var2[i])
                end
            elseif isa(var, VariableRef)
                JuMP.unfix(var2)
            end
        end
    end
    return nothing
end

"""
set_policy!(
    jmodel::JuDGEModel,
    jmodel2::JuDGEModel,
    mapping::Union{Symbol,Dict{AbstractTree,AbstractTree}})

Fixes the policy of a JuDGEModel object based on another JuDGEModel object.

### Required Arguments
`jmodel` is the JuDGE model for which we wish to set the policy.
`jmodel2` is the JuDGE model from which we wish to copy the policy.
`mapping` is can either be set to the symbol `:by_depth` or `:by_nodeID` or be an explicit
dictionary mapping the nodes in `jmodel.tree` the nodes in `jmodel2.tree`. For any node in `jmodel.tree`
that is not mapped, no policy is set.
"""
function set_policy!(
    jmodel::JuDGEModel,
    jmodel2::JuDGEModel,
    mapping::Union{Symbol,Dict{AbstractTree,AbstractTree}},
)
    function set_var(
        var,
        var2,
        hist::Vector{AbstractTree},
        mapping,
        name,
        node::AbstractTree,
        options,
        i,
        rounded::Bool,
    )
        if var2 !== nothing
            if rounded
                if i == 0
                    val2 = round(JuMP.value(var2))
                else
                    val2 = round(JuMP.value(var2[i]))
                end
            else
                if i == 0
                    val2 = JuMP.value(var2)
                else
                    val2 = JuMP.value(var2[i])
                end
            end
            if i == 0
                JuMP.fix(var, val2, force = true)
            else
                JuMP.fix(var[i], val2, force = true)
            end
        end
        val = nothing
        if options[1] == :cumulative
            missing_node = false
            for n in hist
                if haskey(mapping, n)
                    if val === nothing
                        val = 0.0
                    end

                    v =
                        jmodel2.master_problem.ext[:expansions][mapping[n]][name]
                    if rounded
                        if i == 0
                            val += round(JuMP.value(v))
                        else
                            val += round(JuMP.value(v[i]))
                        end
                    else
                        if i == 0
                            val += JuMP.value(v)
                        else
                            val += JuMP.value(v[i])
                        end
                    end
                else
                    missing_node = true
                    break
                end
            end
            if missing_node
                val = nothing
            end
        elseif options[1] == :state
            if node.parent !== nothing
                if haskey(mapping, node) && haskey(mapping, node.parent)
                    v =
                        i == 0 ?
                        jmodel2.master_problem.ext[:expansions][mapping[node]][name] :
                        jmodel2.master_problem.ext[:expansions][mapping[node]][name][i]
                    v2 =
                        i == 0 ?
                        jmodel2.master_problem.ext[:expansions][mapping[node.parent]][name] :
                        jmodel2.master_problem.ext[:expansions][mapping[node.parent]][name][i]
                    if rounded
                        val = round(JuMP.value(v) - JuMP.value(v2))
                    else
                        val = JuMP.value(v) - JuMP.value(v2)
                    end
                end
            else
                if haskey(mapping, node)
                    v =
                        i == 0 ?
                        jmodel2.master_problem.ext[:expansions][mapping[node]][name] :
                        jmodel2.master_problem.ext[:expansions][mapping[node]][name][i]
                    if rounded
                        val = round(JuMP.value(v) - options[7])
                    else
                        val = JuMP.value(v) - options[7]
                    end
                end
            end
        end

        slacks =
            i == 0 ? jmodel.master_problem.ext[:cover_slacks][node][name] :
            jmodel.master_problem.ext[:cover_slacks][node][name][i]

        sp_var =
            i == 0 ? jmodel.sub_problems[node].ext[:expansions][name] :
            jmodel.sub_problems[node].ext[:expansions][name][i]

        bc = nothing

        if val !== nothing
            if haskey(slacks, 1) || haskey(slacks, 2)
                if !haskey(slacks, 2)
                    bc = BranchConstraint(sp_var, :le, val, node)
                elseif !haskey(slacks, 1)
                    bc = BranchConstraint(sp_var, :ge, val, node)
                else
                    bc = nothing
                end
            else
                bc = BranchConstraint(sp_var, :eq, val, node)
            end
        end
        return bc
    end

    if termination_status(jmodel2.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel2.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel2.master_problem) != MOI.LOCALLY_SOLVED &&
       termination_status(jmodel2.master_problem) != MOI.INTERRUPTED
        error("You need to first solve the model that sets the policy.")
    end

    if typeof(mapping) == Symbol
        mode = mapping
        mapping = Dict{AbstractTree,AbstractTree}()
        if mode == :by_depth
            if length(get_leafnodes(jmodel2.tree)) == 1
                for node in collect(jmodel.tree)
                    dpth = depth(node)
                    for node2 in collect(jmodel2.tree)
                        if depth(node2) == dpth
                            mapping[node] = node2
                            break
                        end
                    end
                end
            else
                error("Cannot map nodes by depth since this is ambiguous.")
            end
        elseif mode == :by_nodeID
            for node in collect(jmodel.tree)
                for node2 in collect(jmodel2.tree)
                    if getID(node) == getID(node2)
                        mapping[node] = node2
                        break
                    end
                end
            end
        else
            error("Invalid mapping mode. Use ':by_depth' or ':by_nodeID'.")
        end
    end
    bcs = BranchConstraint[]
    for node in collect(jmodel.tree)
        node2 = nothing
        if haskey(mapping, node)
            node2 = mapping[node]
            #continue
        end

        for (name, var) in jmodel.master_problem.ext[:expansions][node]
            var2 =
                node2 === nothing ? nothing :
                jmodel2.master_problem.ext[:expansions][node2][name]
            options = jmodel.sub_problems[jmodel.tree].ext[:options][name]
            i_min = max(1, depth(node) - options[3] - options[2] + 2)
            i_max = depth(node) + 1 - options[2]
            interval = i_min:i_max

            pre = reverse(history(node))
            hist = AbstractTree[]
            for i in interval
                push!(hist, pre[i])
            end

            if options[4] == :Con
                if typeof(var) <: AbstractArray
                    for i in get_keys(var)
                        bc = set_var(
                            var,
                            var2,
                            hist,
                            mapping,
                            name,
                            node,
                            options,
                            i,
                            false,
                        )
                        if bc !== nothing
                            push!(bcs, bc)
                        end
                    end
                elseif isa(var, VariableRef)
                    bc = set_var(
                        var,
                        var2,
                        hist,
                        mapping,
                        name,
                        node,
                        options,
                        0,
                        false,
                    )
                    if bc !== nothing
                        push!(bcs, bc)
                    end
                end
            else
                if typeof(var) <: AbstractArray
                    for i in get_keys(var)
                        bc = set_var(
                            var,
                            var2,
                            hist,
                            mapping,
                            name,
                            node,
                            options,
                            i,
                            true,
                        )
                        if bc !== nothing
                            push!(bcs, bc)
                        end
                    end
                elseif isa(var, VariableRef)
                    bc = set_var(
                        var,
                        var2,
                        hist,
                        mapping,
                        name,
                        node,
                        options,
                        0,
                        true,
                    )
                    if bc !== nothing
                        push!(bcs, bc)
                    end
                end
            end
        end
    end
    jmodel.ext[:branches] = [Branch(bcs)]
    return
end

function resolve_fixed(jmodel::JuDGEModel)
    for n in collect(jmodel.tree)
        JuMP.optimize!(jmodel.sub_problems[n])
    end

    scenario_objs = Dict{Leaf,Float64}()
    for (leaf, con) in jmodel.master_problem.ext[:scenprofit_con]
        obj = 0.0
        for node in history(leaf)
            obj += objective_value(jmodel.sub_problems[node])
            for key in keys(jmodel.master_problem.ext[:expansions][node])
                var = jmodel.master_problem.ext[:expansions][node][key]
                slacks = jmodel.master_problem.ext[:cover_slacks][node][key]
                if isa(var, VariableRef)
                    obj += JuMP.value(var) * normalized_coefficient(con, var)
                    if haskey(slacks, 1) && JuMP.value(slacks[1]) > 10e-8
                        obj +=
                            JuMP.value(slacks[1]) *
                            normalized_coefficient(con, slacks[1])
                    end
                    if haskey(slacks, 2) && JuMP.value(slacks[2]) > 10e-8
                        obj +=
                            JuMP.value(slacks[2]) *
                            normalized_coefficient(con, slacks[2])
                    end
                elseif typeof(var) <: AbstractArray
                    for v in get_keys(var)
                        obj +=
                            JuMP.value(var[v]) *
                            normalized_coefficient(con, var[v])
                        if haskey(slacks[v], 1)
                            obj +=
                                JuMP.value(slacks[v][1]) *
                                normalized_coefficient(con, slacks[v][1])
                        end
                        if haskey(slacks[v], 2)
                            obj +=
                                JuMP.value(slacks[v][2]) *
                                normalized_coefficient(con, slacks[v][2])
                        end
                    end
                end
            end
        end
        scenario_objs[leaf] = obj
    end

    return compute_objval(scenario_objs, jmodel.probabilities, jmodel.risk)
end
