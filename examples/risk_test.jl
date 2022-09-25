using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function transport_risk_test(
    tree::AbstractTree,
    supply_nodes::Vector{Symbol},
    demand_nodes::Vector{Symbol},
    demand::Dict{AbstractTree,Dict{Symbol,T}} where {T<:Real},
    model_type::Symbol,
)
    cost = Dict(
        zip(
            [(supply_nodes[i], demand_nodes[j]) for i in 1:3, j in 1:3],
            [3 4 5; 4 3 4; 5 4 3],
        ),
    )

    function sub_problem(n::AbstractTree)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @expansion(sp, 0 <= small_upgrade[supply_nodes] <= 2, lag = 0, ub = 1)
        @expansion(sp, 0 <= large_upgrade[supply_nodes] <= 1, lag = 1)

        @capitalcosts(sp, 9 * sum(small_upgrade) + 20 * sum(large_upgrade))

        @variable(sp, flow[supply_nodes, demand_nodes] >= 0)

        @constraint(
            sp,
            supply_limit[i in supply_nodes],
            sum(flow[i, j] for j in demand_nodes) <=
            10 + 10 * small_upgrade[i] + 30 * large_upgrade[i]
        )
        @constraint(
            sp,
            demand_required[j in demand_nodes],
            sum(flow[i, j] for i in supply_nodes) == demand[n][j]
        )

        @objective(
            sp,
            Min,
            sum(
                cost[i, j] * flow[i, j] for i in supply_nodes, j in demand_nodes
            )
        )

        return sp
    end

    if model_type == :deteq
        model = DetEqModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_DE_Solver,
            risk = Risk(0.19, 1 / 3),
        )
    elseif model_type == :decomp
        model = JuDGEModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_MP_Solver,
            risk = Risk(0.19, 1 / 3),
        )
    else
        error("Invalid type of model.")
    end
    JuDGE.solve(model)
    JuDGE.print_expansions(model, inttol = 1e-7)
    JuDGE.get_scen_objs(model)
    pr = JuDGE.get_risk_probs(model)
    obj1 = JuDGE.get_objval(model)

    if model_type == :deteq
        model2 = DetEqModel(tree, pr, sub_problem, JuDGE_DE_Solver)
    else
        model2 = JuDGEModel(tree, pr, sub_problem, JuDGE_MP_Solver)
    end

    JuDGE.solve(model2)
    JuDGE.print_expansions(model2, inttol = 1e-7)
    JuDGE.get_scen_objs(model2)
    obj2 = JuDGE.get_objval(model2)

    return obj1, obj2
end

function transport_risk_test2(
    tree::AbstractTree,
    supply_nodes::Vector{Symbol},
    demand_nodes::Vector{Symbol},
    demand::Dict{AbstractTree,Dict{Symbol,T}} where {T<:Real},
    model_type::Symbol,
)
    cost = Dict(
        zip(
            [(supply_nodes[i], demand_nodes[j]) for i in 1:3, j in 1:3],
            [3 4 5; 4 3 4; 5 4 3],
        ),
    )

    function sub_problem(n::AbstractTree)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @expansion(
            sp,
            0 <= small_upgrade[supply_nodes] <= 2,
            Int,
            lag = 0,
            ub = 1
        )
        @expansion(sp, large_upgrade[supply_nodes], Bin, lag = 1)

        @capitalcosts(sp, 9 * sum(small_upgrade) + 20 * sum(large_upgrade))

        @variable(sp, flow[supply_nodes, demand_nodes] >= 0)

        @constraint(
            sp,
            supply_limit[i in supply_nodes],
            sum(flow[i, j] for j in demand_nodes) <=
            10 + 10 * small_upgrade[i] + 30 * large_upgrade[i]
        )
        @constraint(
            sp,
            demand_required[j in demand_nodes],
            sum(flow[i, j] for i in supply_nodes) == demand[n][j]
        )

        @objective(
            sp,
            Min,
            sum(
                cost[i, j] * flow[i, j] for i in supply_nodes, j in demand_nodes
            )
        )

        return sp
    end

    if model_type == :deteq
        model = DetEqModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_DE_Solver,
        )
        model2 = DetEqModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_DE_Solver,
            risk = Risk(0.9, 1 / 3),
        )
        JuDGE.solve(model)
        JuDGE.solve(model2)
    elseif model_type == :decomp
        model = JuDGEModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_MP_Solver,
        )
        model2 = JuDGEModel(
            tree,
            ConditionallyUniformProbabilities,
            sub_problem,
            JuDGE_MP_Solver,
            risk = Risk(0.9, 1 / 3),
        )
        model = JuDGE.branch_and_price(model)
        model2 = JuDGE.branch_and_price(model2)
    else
        error("Invalid type of model.")
    end

    return maximum(values(JuDGE.get_regret(model2, model))),
    JuDGE.get_risk_probs(model2)[tree[1, 2]],
    sum(sum.(JuDGE.scenarios_CDF(model2)))
end
