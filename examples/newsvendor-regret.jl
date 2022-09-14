using JuMP, JuDGE, Plots

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function newsvendor_regret(;
    depth = 1,
    cost = 5.0,
    price = 8.0,
    demands = [50.0, 60.0, 70.0],
    CVaR = RiskNeutral(),
    mytree = nothing,
)
    if mytree == nothing
        mytree = narytree(depth, length(demands))
    end

    demand = Dict{AbstractTree,Float64}()
    nodes = collect(mytree, order = :breadth)
    for i in eachindex(nodes)
        if i == 1
            demand[nodes[i]] = 0
        else
            demand[nodes[i]] = demands[(i-2)%length(demands)+1]
        end
    end

    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)
        @expansion(model, 0 <= papers_ordered <= 1000, lag = 1, duration = 1)
        @capitalcosts(model, papers_ordered * cost)
        @variable(model, sales >= 0)#, Int)
        @constraint(model, maxsale1, sales <= papers_ordered)
        @constraint(model, maxsale2, sales <= demand[node])
        @objective(model, Min, -price * sales)
        return model
    end

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        risk = CVaR,
        check = false,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(objective_value(deteq.problem)),
    )

    return deteq
end

tree = narytree(5, 3)

deteq1 = newsvendor_regret(
    depth = 5,
    cost = 5.0,
    price = 8.0,
    demands = [10, 20, 30],
    CVaR = RiskNeutral(),
    mytree = tree,
)

deteq2 = newsvendor_regret(
    depth = 5,
    cost = 5.0,
    price = 8.0,
    demands = [10, 20, 30],
    CVaR = Risk(0.6, 0.1),
    mytree = tree,
)

deteq3 = newsvendor_regret(
    depth = 5,
    cost = 5.0,
    price = 8.0,
    demands = [10, 20, 30],
    CVaR = [
        Risk(0.3, 0.1),
        Risk(0.7, 0.5, offset = JuDGE.get_scen_objs(deteq1)),
    ],
    mytree = tree,
)

#plotly()
pyplot()

plot(JuDGE.scenarios_CDF(deteq1, tol = 1e-8))
plot!(JuDGE.scenarios_CDF(deteq2, tol = 1e-8))
plot!(JuDGE.scenarios_CDF(deteq3, tol = 1e-8))

riskprobs = JuDGE.get_risk_probs(deteq3)
scenprof = JuDGE.get_scen_objs(deteq3)
regret = JuDGE.get_regret(deteq3, deteq1)

plot(
    [(scenprof[leaf], riskprobs[leaf]) for leaf in keys(scenprof)],
    seriestype = :scatter,
)

plot(
    [(scenprof[leaf], regret[leaf]) for leaf in keys(scenprof)],
    marker_z = [riskprobs[leaf] for leaf in keys(scenprof)],
    color = :jet,
    alpha = 0.3,
    seriestype = :scatter,
)

plot(
    [
        (scenprof[leaf], regret[leaf], riskprobs[leaf]) for
        leaf in keys(scenprof)
    ],
    marker_z = [riskprobs[leaf] for leaf in keys(scenprof)],
    color = :jet,
    alpha = 0.3,
    seriestype = :scatter,
)

solution = JuDGE.solution_to_dictionary(deteq3)
JuDGE.add_to_dictionary!(solution, riskprobs, :marginal_pr)
JuDGE.add_to_dictionary!(
    solution,
    JuDGE.get_risk_probs(deteq3, :conditional),
    :conditional_pr,
)
JuDGE.visualize_tree(tree, solution, filename = "newsvendor")
