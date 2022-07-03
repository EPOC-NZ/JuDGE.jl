using Random
using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function createRandomArray(seed::Int, n::Int, a::Int, b::Int)
    Random.seed!(seed)

    array = []
    for i in 1:n
        push!(array, rand(a))
    end

    for i in 1:2*n
        push!(array, rand(b))
    end

    return array
end

function setup_model(seed::Int, numitems::Int, array)

    # how many investments?
    numinvest = 2

    # size of tree?
    degree = 3
    height = 3

    mytree = narytree(height, degree)
    probabilities = ConditionallyUniformProbabilities(mytree)
    exptree = narytree(height, 1)
    nodes = collect(mytree)
    totalnodes = length(nodes)

    if array == nothing
        array = createRandomArray(seed, totalnodes, numinvest, numitems)
    end

    # create a dictionary to store all the nodal data
    data = Dict{Symbol,Dict{AbstractTree,Any}}()

    data[:investcost] = Dict{AbstractTree,Any}()
    for i in 1:totalnodes
        data[:investcost][nodes[i]] =
            (popfirst!(array) * 2 + 2 * [2.0, 3.5]) *
            (1 - ((i - 1) / (totalnodes * 1.2))) / 2
    end

    investvol = [20, 25]
    initialcap = 80

    data[:itemvolume] = Dict{AbstractTree,Any}()
    for i in 1:totalnodes
        data[:itemvolume][nodes[i]] =
            (
                ((popfirst!(array) .- 0.5) * 2) * 2 +
                collect(range(4, 22, length = numitems))
            ) * 20 / numitems
    end

    data[:itemcost] = Dict{AbstractTree,Any}()
    for i in 1:totalnodes
        data[:itemcost][nodes[i]] =
            ((popfirst!(array) .- 0.5) * 2) * 0.5 +
            collect(range(0.5, 1, length = numitems))
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, 0 <= bag[1:numinvest] <= 100, Int)
        @capitalcosts(
            model,
            sum(data[:investcost][node][i] * bag[i] for i in 1:numinvest)
        )
        @variable(model, y[1:numitems], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * data[:itemvolume][node][i] for i in 1:numitems) <=
            initialcap + sum(bag[i] * investvol[i] for i in 1:numinvest)
        )
        @objective(
            model,
            Min,
            sum(-data[:itemcost][node][i] * y[i] for i in 1:numitems)
        )
        return model
    end

    return (mytree, probabilities, sub_problems)
end

function policy_comparison(
    tree::AbstractTree,
    probabilities::Dict{AbstractTree,Float64},
    sub_problems::T where {T<:Function},
    visualise::Bool,
    reltol::Float64,
)
    height = JuDGE.depth(collect(tree)[end])

    # formulates a separate JuDGEModel for each path to a leaf node
    @info("Setting up wait-and-see problems")
    WS, WS_pr = JuDGEModel(
        tree,
        probabilities,
        sub_problems,
        JuDGE_MP_Solver,
        perfect_foresight = true,
    )

    # solve all the wait and see problems
    @info("Solving wait-and-see problems")
    for (leaf, ws) in WS
        JuDGE.solve(ws, verbose = 0)
    end

    # formulate a standard JuDGEModel
    @info("Setting up here-and-now problem")
    HN = JuDGEModel(tree, probabilities, sub_problems, JuDGE_MP_Solver)

    # solve the standard JuDGEMdoel using branch-and-price
    @info("Solving here-and-now problem")
    JuDGE.branch_and_price(
        HN,
        termination = Termination(inttol = 10^-7, reltol = reltol),
        verbose = 0,
    )

    # create a cloned tree that has additional branches representing the conditional expectation of the future data
    tree2 = JuDGE.append_expected_branches(tree)

    # add data to the dictionary for the new nodes in the tree
    JuDGE.define_conditional_means!(tree2, probabilities, sub_problems.data)

    vis = Dict{AbstractTree,Dict{Symbol,Any}}()
    JuDGE.add_to_dictionary!(vis, sub_problems.data[:investcost], :investcost)
    JuDGE.add_to_dictionary!(vis, sub_problems.data[:itemvolume], :itemvolume)
    JuDGE.add_to_dictionary!(vis, sub_problems.data[:itemcost], :itemcost)
    if visualise
        JuDGE.visualize_tree(
            tree2,
            vis,
            filename = "augmented",
            rel_angle = true,
        )
    else
        JuDGE.remove_from_dictionary!(vis, :itemcost)
        JuDGE.get_active_columns(HN; inttol = 10^-7)
    end
    # code for EEV and rolling horizon
    RH = nothing
    RH_old = nothing
    EEV = nothing
    for iter in 0:height
        # based on the iteration, probabilities within the tree are updated, and nodes with probability 0 are removed
        tree3, pr = JuDGE.refine_tree(tree2, probabilities, iter)

        # visualise the trees that are used to find and simulate the EV and RH policies
        if visualise
            vis = Dict{AbstractTree,Dict{Symbol,Any}}()
            JuDGE.add_to_dictionary!(
                vis,
                sub_problems.data[:investcost],
                :investcost,
            )
            JuDGE.add_to_dictionary!(
                vis,
                sub_problems.data[:itemvolume],
                :itemvolume,
            )
            JuDGE.add_to_dictionary!(
                vis,
                sub_problems.data[:itemcost],
                :itemcost,
            )
            #JuDGE.add_to_dictionary!(vis, pr, :probability)
            JuDGE.visualize_tree(
                tree3,
                vis,
                filename = "iter" * string(iter),
                rel_angle = true,
            )
        end
        # in iteration 1, we can find the EEV by simulating the solution from iteration 0, in the real tree
        if iter == 1
            EEV = JuDGEModel(tree, probabilities, sub_problems, JuDGE_MP_Solver)
            JuDGE.set_policy!(EEV, RH_old, :by_depth)
            JuDGE.branch_and_price(
                EEV,
                termination = Termination(inttol = 10^-7, reltol = reltol),
                verbose = 0,
            )
        end

        RH = JuDGEModel(tree3, pr, sub_problems, JuDGE_MP_Solver)

        if RH_old != nothing
            JuDGE.set_policy!(RH, RH_old, :by_nodeID)
        end

        JuDGE.branch_and_price(
            RH,
            termination = Termination(inttol = 10^-7, reltol = reltol),
            verbose = 0,
        )
        RH_old = RH
    end

    @info("HN obj: " * string(JuDGE.get_objval(HN)))
    @info(
        "WS obj: " *
        string(sum(JuDGE.get_objval(ws) * WS_pr[node] for (node, ws) in WS))
    )
    @info("EEV obj: " * string(JuDGE.get_objval(EEV)))
    @info("RH obj: " * string(JuDGE.get_objval(RH)))
    @info(
        "EVPI: " * string(
            JuDGE.get_objval(HN) -
            sum(JuDGE.get_objval(ws) * WS_pr[node] for (node, ws) in WS),
        )
    )
    @info(
        "VSS (EV policy): " *
        string(JuDGE.get_objval(EEV) - JuDGE.get_objval(HN))
    )
    @info(
        "VSS (RH policy): " *
        string(JuDGE.get_objval(RH) - JuDGE.get_objval(HN))
    )
    @info("HN expansions:")
    JuDGE.print_expansions(HN, inttol = 10e-6)

    @info("EEV expansions:")
    JuDGE.print_expansions(EEV, inttol = 10e-6)

    @info("RH expansions:")
    JuDGE.print_expansions(RH, inttol = 10e-6)

    return [JuDGE.get_objval(HN), JuDGE.get_objval(EEV), JuDGE.get_objval(RH)]
end

function evpi_vss(
    seed::Int,
    numitems::Int,
    visualise::Bool,
    reltol::Float64;
    array = nothing,
)
    array = deepcopy(array)
    (tree, probabilities, sub_problems) = setup_model(seed, numitems, array)
    return policy_comparison(
        tree,
        probabilities,
        sub_problems,
        visualise,
        reltol,
    )
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    evpi_vss(1, 20, true, 0.001)
end
