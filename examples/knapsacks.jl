using Random
using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

# this test is a test based off of a presentation by Prof. Andy Philpott
function knapsack_fixed()
    mytree = narytree(2, 2)
    function invest_cost(node)
        if node == get_node(mytree, [1])
            180.0
        elseif node == get_node(mytree, [1, 1])
            50.0
        elseif node == get_node(mytree, [1, 2])
            60.0
        elseif node == get_node(mytree, [1, 1, 1])
            40.0
        elseif node == get_node(mytree, [1, 1, 2])
            60.0
        elseif node == get_node(mytree, [1, 2, 1])
            10.0
        elseif node == get_node(mytree, [1, 2, 2])
            10.0
        end
    end

    function item_volume(node)
        if node == get_node(mytree, [1])
            [6, 2, 1, 1, 1]
        elseif node == get_node(mytree, [1, 1])
            [8, 2, 2, 2, 1]
        elseif node == get_node(mytree, [1, 2])
            [8, 1, 1, 1, 3]
        elseif node == get_node(mytree, [1, 1, 1])
            [4, 4, 3, 1, 2]
        elseif node == get_node(mytree, [1, 1, 2])
            [1, 3, 1, 1, 2]
        elseif node == get_node(mytree, [1, 2, 1])
            [7, 3, 1, 1, 1]
        elseif node == get_node(mytree, [1, 2, 2])
            [2, 5, 2, 1, 2]
        end
    end

    function item_reward(node)
        if node == get_node(mytree, [1])
            [60, 20, 10, 15, 10]
        elseif node == get_node(mytree, [1, 1])
            [8, 10, 20, 20, 10]
        elseif node == get_node(mytree, [1, 2])
            [8, 10, 15, 10, 30]
        elseif node == get_node(mytree, [1, 1, 1])
            [40, 40, 35, 10, 20]
        elseif node == get_node(mytree, [1, 1, 2])
            [15, 35, 15, 15, 20]
        elseif node == get_node(mytree, [1, 2, 1])
            [70, 30, 15, 15, 10]
        elseif node == get_node(mytree, [1, 2, 2])
            [25, 50, 25, 15, 20]
        end
    end

    ### with judge
    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, bag, Bin)
        @capitalcosts(model, bag * invest_cost(node))
        @variable(model, y[1:5], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag
        )
        @objective(model, Min, sum(-item_reward(node)[i] * y[i] for i in 1:5))
        return model
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
    )
    JuDGE.solve(judy, verbose = 1)

    println("Objective: " * string(JuDGE.get_objval(judy)))
    JuDGE.print_expansions(judy, onlynonzero = false)

    println("Re-solved Objective: " * string(resolve_subproblems(judy)))

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )

    return objective_value(judy.master_problem)
end

# this test is a test based off of the test in the old JuDGE
function knapsack_random()
    Random.seed!(1)
    # how many investments?
    numinvest = 2

    # number of items to pick from in the knapsack?
    numitems = 20

    # size of tree?
    degree = 3
    height = 3

    totalnodes = Int((degree^(height + 1) - 1) / (degree - 1))

    investcost = zeros(totalnodes, numinvest)
    for i in 1:totalnodes
        investcost[i, :] =
            (rand(numinvest) * 2 + 2 * [2.0, 3.5]) *
            (1 - ((i - 1) / (totalnodes * 1.2)))
    end

    investvol = [40, 50]
    initialcap = 80

    itemvolume = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemvolume[i, :] =
            ((rand(numitems) .- 0.5) * 2) * 2 +
            collect(range(4, 22, length = numitems))
    end

    itemcost = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemcost[i, :] =
            ((rand(numitems) .- 0.5) * 2) * 0.5 +
            collect(range(0.5, 1, length = numitems))
    end

    mytree = narytree(height, degree)

    nodes = collect(mytree)
    function data(node, input)
        return input[findall(x -> x == node, nodes)[1], :]
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, bag[1:numinvest], Bin)
        @capitalcosts(
            model,
            sum(data(node, investcost)[i] * bag[i] for i in 1:numinvest)
        )
        @variable(model, y[1:numitems], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * data(node, itemvolume)[i] for i in 1:numitems) <=
            initialcap + sum(bag[i] * investvol[i] for i in 1:numinvest)
        )
        @objective(
            model,
            Min,
            sum(-data(node, itemcost)[i] * y[i] for i in 1:numitems)
        )
        return model
    end

    function format_output(s::Symbol, values)
        if s == :bag
            return sum(values[i] * investvol[i] for i in 1:numinvest)
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
    )
    JuDGE.solve(judy)

    println("Objective: " * string(JuDGE.get_objval(judy)))
    JuDGE.print_expansions(judy, format = format_output)

    println("Re-solved Objective: " * string(resolve_subproblems(judy)))

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )

    return objective_value(judy.master_problem)
end

function knapsack_branch_and_price()
    Random.seed!(50)
    # how many investments?
    numinvest = 5

    # number of items to pick from in the knapsack?
    numitems = 10

    # size of tree?
    degree = 3
    height = 4

    totalnodes = Int((degree^(height + 1) - 1) / (degree - 1))

    investcost = zeros(totalnodes, numinvest)
    for i in 1:totalnodes
        investcost[i, :] =
            ([1, 1.8, 3.5, 6.8, 13.5]) * (1 - ((i - 1) / (totalnodes * 1.2)))
    end

    investvol = [1, 2, 4, 8, 16]
    initialcap = 0

    itemvolume = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemvolume[i, :] =
            ((rand(numitems)) * 2) + collect(range(4, 22, length = numitems))
    end

    itemcost = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemcost[i, :] = ((rand(numitems) .- 0.5) * 2) * 2# + collect(range(0.5,1,length = numitems))
    end

    mytree = narytree(height, degree)

    nodes = collect(mytree)
    function data(node, input)
        return input[findall(x -> x == node, nodes)[1], :]
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, bag[1:numinvest], Bin)
        @capitalcosts(
            model,
            sum(data(node, investcost)[i] * bag[i] for i in 1:numinvest)
        )
        @variable(model, y[1:numitems], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * data(node, itemvolume)[i] for i in 1:numitems) <=
            initialcap + sum(bag[i] * investvol[i] for i in 1:numinvest)
        )
        @objective(
            model,
            Min,
            sum(-data(node, itemcost)[i] * y[i] for i in 1:numitems)
        )
        return model
    end

    function format_output(s::Symbol, values)
        if s == :bag
            return sum(values[i] * investvol[i] for i in 1:numinvest)
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
    )

    best = JuDGE.branch_and_price(
        judy,
        termination = Termination(rlx_abstol = 10^-7, inttol = 10^-6),
    )

    println("Objective: " * string(JuDGE.get_objval(best)))
    JuDGE.print_expansions(best, format = format_output)

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )

    return objective_value(best.master_problem)
end

function knapsack_risk_averse()
    Random.seed!(50)
    # how many investments?
    numinvest = 5

    # number of items to pick from in the knapsack?
    numitems = 10

    # size of tree?
    degree = 3
    height = 4

    totalnodes = Int((degree^(height + 1) - 1) / (degree - 1))

    investcost = zeros(totalnodes, numinvest)
    for i in 1:totalnodes
        investcost[i, :] =
            ([1, 1.8, 3.5, 6.8, 13.5]) * (1 - ((i - 1) / (totalnodes * 1.2)))
    end

    investvol = [1, 2, 4, 8, 16]
    initialcap = 0

    itemvolume = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemvolume[i, :] =
            ((rand(numitems)) * 2) + collect(range(4, 22, length = numitems))
    end

    itemcost = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemcost[i, :] = ((rand(numitems) .- 0.5) * 2) * 2# + collect(range(0.5,1,length = numitems))
    end

    mytree = narytree(height, degree)

    nodes = collect(mytree)
    function data(node, input)
        return input[findall(x -> x == node, nodes)[1], :]
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, bag[1:numinvest], Bin)
        @capitalcosts(
            model,
            sum(data(node, investcost)[i] * bag[i] for i in 1:numinvest)
        )
        @variable(model, y[1:numitems], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * data(node, itemvolume)[i] for i in 1:numitems) <=
            initialcap + sum(bag[i] * investvol[i] for i in 1:numinvest)
        )
        @objective(
            model,
            Min,
            sum(-data(node, itemcost)[i] * y[i] for i in 1:numitems)
        )
        return model
    end

    function format_output(s::Symbol, values)
        if s == :bag
            return sum(values[i] * investvol[i] for i in 1:numinvest)
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        risk = Risk(0.5, 0.05),
    )

    best = JuDGE.branch_and_price(
        judy,
        termination = Termination(rlx_abstol = 10^-6, inttol = 10^-6),
    )
    println("Objective: " * string(JuDGE.get_objval(best)))
    JuDGE.print_expansions(best, format = format_output)

    println("Re-solved Objective: " * string(resolve_subproblems(best)))

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        risk = Risk(0.5, 0.05),
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )
    JuDGE.print_expansions(deteq, format = format_output)

    return objective_value(best.master_problem)
end

function knapsack_delayed_investment(; CVaR = RiskNeutral())
    Random.seed!(100)
    # how many investments?
    numinvest = 2

    # number of items to pick from in the knapsack?
    numitems = 25

    # size of tree?
    degree = 3
    height = 3

    totalnodes = Int((degree^(height + 1) - 1) / (degree - 1))

    investcost = zeros(totalnodes, numinvest)
    for i in 1:totalnodes
        investcost[i, :] =
            (rand(numinvest) * 2 + 2 * [2.0, 3.5]) * (1 - ((i) / (totalnodes)))
    end

    # investvol = [40,45,50,70]
    investvol = [40, 50]
    initialcap = 80

    itemvolume = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemvolume[i, :] =
            ((rand(numitems) .- 0.5) * 2) * 2 +
            collect(range(4, 22, length = numitems))
    end

    itemcost = zeros(totalnodes, numitems)
    for i in 1:totalnodes
        itemcost[i, :] =
            ((rand(numitems) .- 0.5) * 2) * 0.5 +
            collect(range(0.5, 1, length = numitems))
    end

    mytree = narytree(height, degree)

    nodes = collect(mytree)
    function data(node, input)
        return input[findall(x -> x == node, nodes)[1], :]
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, bag[1:numinvest], Bin, lag = 1)
        @capitalcosts(
            model,
            sum(data(node, investcost)[i] * bag[i] for i in 1:numinvest)
        )
        @variable(model, y[1:numitems], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * data(node, itemvolume)[i] for i in 1:numitems) <=
            initialcap + sum(bag[i] * investvol[i] for i in 1:numinvest)
        )
        @objective(
            model,
            Min,
            sum(-data(node, itemcost)[i] * y[i] for i in 1:numitems)
        )
        return model
    end

    function format_output(s::Symbol, values)
        if s == :bag_bought
            return sum(values[i] * investvol[i] for i in 1:numinvest)
        elseif s == :bag_arrived
            return 0.0
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        risk = CVaR,
    )

    judy = JuDGE.branch_and_price(
        judy,
        termination = Termination(rlx_abstol = 10^-6, inttol = 10^-6),
    )
    println("Objective: " * string(JuDGE.get_objval(judy, risk = CVaR)))
    JuDGE.print_expansions(judy, format = format_output)

    println("Re-solved Objective: " * string(resolve_subproblems(judy)))

    solution = JuDGE.solution_to_dictionary(judy)
    JuDGE.visualize_tree(mytree, solution, filename = "knapsack")

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        risk = CVaR,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq, risk = CVaR)),
    )
    return objective_value(judy.master_problem)
end

function knapsack_shutdown()
    mytree = narytree(2, 2)

    divest_revenue = Dict{AbstractTree,Float64}()
    divest_revenue[get_node(mytree, [1])] = 60.0
    divest_revenue[get_node(mytree, [1, 1])] = 58.0
    divest_revenue[get_node(mytree, [1, 2])] = 35.0
    divest_revenue[get_node(mytree, [1, 1, 1])] = 40.0
    divest_revenue[get_node(mytree, [1, 1, 2])] = 35.0
    divest_revenue[get_node(mytree, [1, 2, 1])] = 25.0
    divest_revenue[get_node(mytree, [1, 2, 2])] = 15.0

    item_volume = Dict{AbstractTree,Array{Float64,1}}()
    item_volume[get_node(mytree, [1])] = [6, 2, 1, 1, 1]
    item_volume[get_node(mytree, [1, 1])] = [8, 2, 2, 2, 1]
    item_volume[get_node(mytree, [1, 2])] = [8, 1, 1, 1, 3]
    item_volume[get_node(mytree, [1, 1, 1])] = [4, 4, 3, 1, 2]
    item_volume[get_node(mytree, [1, 1, 2])] = [1, 3, 1, 1, 2]
    item_volume[get_node(mytree, [1, 2, 1])] = [7, 3, 1, 1, 1]
    item_volume[get_node(mytree, [1, 2, 2])] = [2, 5, 2, 1, 2]

    item_reward = Dict{AbstractTree,Array{Float64,1}}()
    item_reward[get_node(mytree, [1])] = [60, 20, 10, 15, 10]
    item_reward[get_node(mytree, [1, 1])] = [8, 10, 20, 20, 10]
    item_reward[get_node(mytree, [1, 2])] = [8, 10, 15, 10, 30]
    item_reward[get_node(mytree, [1, 1, 1])] = [40, 40, 35, 10, 20]
    item_reward[get_node(mytree, [1, 1, 2])] = [15, 35, 15, 15, 20]
    item_reward[get_node(mytree, [1, 2, 1])] = [70, 30, 15, 15, 10]
    item_reward[get_node(mytree, [1, 2, 2])] = [25, 50, 25, 15, 20]

    ### with judge
    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @shutdown(model, bag, Bin)
        @capitalcosts(model, -bag * divest_revenue[node])
        @variable(model, y[1:5], Bin)
        @constraint(
            model,
            BagExtension,
            sum(y[i] * item_volume[node][i] for i in 1:5) <= 4 - 2 * bag
        )
        @objective(model, Min, sum(-item_reward[node][i] * y[i] for i in 1:5))
        return model
    end

    function format_output(s::Symbol, value)
        if s == :bag
            return -2.0 * value
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
    )
    JuDGE.solve(judy, verbose = 1)

    println("Objective: " * string(JuDGE.get_objval(judy)))
    JuDGE.print_expansions(judy, format = format_output)

    println("Re-solved Objective: " * string(resolve_subproblems(judy)))

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )

    return objective_value(judy.master_problem)
end

function knapsack_budget()
    mytree = narytree(2, 2)

    invest_cost =
        Dict(zip(collect(mytree, order = :breadth), [15, 8, 8, 4, 4, 4, 4]))

    item_volume = Dict(
        zip(
            collect(mytree, order = :breadth),
            [
                [4, 3, 3, 1, 2],
                [5, 3, 4, 2, 1],
                [5, 4, 2, 7, 2],
                [5, 4, 1, 8, 2],
                [3, 1, 5, 6, 3],
                [2, 5, 8, 4, 6],
                [7, 5, 4, 2, 3],
            ],
        ),
    )

    item_reward = Dict(
        zip(
            collect(mytree, order = :breadth),
            [
                [32, 9, 9, 4, 8],
                [30, 12, 40, 10, 9],
                [20, 28, 12, 42, 12],
                [40, 28, 9, 24, 10],
                [15, 7, 20, 48, 12],
                [10, 30, 54, 32, 30],
                [32, 25, 24, 14, 24],
            ],
        ),
    )

    num_items = 5
    num_invest = 6
    initial_volume = 6
    invest_volume = [2, 2, 2, 3, 3, 3]

    ### with judge
    function sub_problems(node)
        sp = Model(JuDGE_SP_Solver)
        @expansion(sp, invest[1:num_invest], Bin)
        @capitalcosts(
            sp,
            sum(invest[i] * invest_volume[i] for i in 1:num_invest) *
            invest_cost[node]
        )
        @variable(sp, y[1:num_items], Bin)
        @constraint(
            sp,
            BagExtension,
            sum(y[i] * item_volume[node][i] for i in 1:num_items) <=
            initial_volume +
            sum(invest_volume[i] * invest[i] for i in 1:num_invest)
        )
        @objective(
            sp,
            Min,
            sum(-item_reward[node][i] * y[i] for i in 1:num_items)
        )
        return sp
    end

    function format_output(s::Symbol, value)
        if s == :invest
            return sum(invest_volume[i] * value[i] for i in 1:num_invest)
        end
        return nothing
    end

    function budget(model, tree)
        for node in collect(tree)
            @constraint(
                model,
                sum(
                    invest_cost[node] * invest_volume[i] * invest[node][i] for
                    i in 1:num_invest
                ) <= 40
            )
        end
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        sideconstraints = budget,
    )
    #JuDGE.solve(judy)
    judy = JuDGE.branch_and_price(judy, verbose = 0)
    println("Objective: " * string(JuDGE.get_objval(judy)))
    JuDGE.print_expansions(judy, format = format_output)

    println("Re-solved Objective: " * string(resolve_subproblems(judy)))

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        sideconstraints = budget,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(JuDGE.get_objval(deteq)),
    )

    return objective_value(judy.master_problem)
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    knapsack_fixed()
    knapsack_random()
    knapsack_branch_and_price()
    knapsack_risk_averse()
    knapsack_delayed_investment()
    knapsack_delayed_investment(CVaR = Risk(0.95, 0.05))
    knapsack_shutdown()
    knapsack_budget()
end
