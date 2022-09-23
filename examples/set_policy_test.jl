using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function test_set_policy(
    D1::Vector{Int},
    D2::Vector{Int},
    lag::Int,
    duration::Int,
)
    function subproblem(node)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @expansion(sp, Z, Bin, lag = lag, duration = duration)

        @variable(sp, x >= 0)
        @variable(sp, y >= 0)

        @constraint(sp, x <= Z)
        @constraint(sp, x + y == D[node])

        #@capitalcosts(sp,Z*10)
        @ongoingcosts(sp, Z * 40 * (JuDGE.depth(node) + 1.2))
        @objective(sp, Min, 1000.0y)

        return sp
    end

    tree = narytree(3, 1)
    D = Dict(zip(collect(tree), D1))
    model = JuDGEModel(
        tree,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model)
    resolve_subproblems(model)

    tree2 = narytree(3, 2)
    D = Dict(zip(collect(tree2), D2))
    model2 = JuDGEModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model2)
    obj2 = JuDGE.get_objval(model2)

    model3 = JuDGEModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )
    JuDGE.set_policy!(model3, model, :by_depth)
    JuDGE.solve(model3)
    obj3 = JuDGE.get_objval(model3)

    return obj3 - obj2
end

function test_set_policy_state(D1::Vector{Int}, D2::Vector{Int})
    function subproblem(node)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @state(sp, -1 <= Z, initial = 2, lb = 0.0)

        @variable(sp, x >= 0)
        @variable(sp, y >= 0)

        @constraint(sp, -x >= Z)
        @constraint(sp, x + y == D[node])

        @objective(sp, Min, x + 1000.0y)

        return sp
    end

    tree = narytree(3, 1)
    D = Dict(zip(collect(tree), D1))
    model = JuDGEModel(
        tree,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )

    JuDGE.solve(model)
    resolve_subproblems(model)

    tree2 = narytree(3, 2)
    D = Dict(zip(collect(tree2), D2))
    model2 = JuDGEModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model2)
    obj2 = JuDGE.get_objval(model2)

    model3 = JuDGEModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_MP_Solver,
        discount_factor = 0.5,
    )
    JuDGE.set_policy!(model3, model, :by_depth)
    JuDGE.solve(model3)
    obj3 = JuDGE.get_objval(model3)

    return obj3 - obj2
end

function test_set_policy2(
    D1::Vector{Int},
    D2::Vector{Int},
    lag::Int,
    duration::Int,
)
    function subproblem(node)
        sp = JuMP.Model()

        @expansion(sp, Z, Bin, lag = lag, duration = duration)

        @variable(sp, x >= 0)
        @variable(sp, y >= 0)

        @constraint(sp, x <= Z)
        @constraint(sp, x + y == D[node])

        #@capitalcosts(sp,Z*10)
        @ongoingcosts(sp, Z * 40 * (JuDGE.depth(node) + 1.2))
        @objective(sp, Min, 1000.0y)

        return sp
    end

    tree = narytree(3, 1)
    D = Dict(zip(collect(tree), D1))
    model = DetEqModel(
        tree,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model)

    tree2 = narytree(3, 2)
    D = Dict(zip(collect(tree2), D2))
    model2 = DetEqModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model2)
    obj2 = JuDGE.get_objval(model2)

    model3 = DetEqModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )
    JuDGE.set_policy!(model3, model, :by_depth)
    JuDGE.solve(model3)
    obj3 = JuDGE.get_objval(model3)

    return obj3 - obj2
end

function test_set_policy_state2(D1::Vector{Int}, D2::Vector{Int})
    function subproblem(node)
        sp = JuMP.Model()

        @state(sp, -1 <= Z, initial = 2, lb = 0.0)

        @variable(sp, x >= 0)
        @variable(sp, y >= 0)

        @constraint(sp, -x >= Z)
        @constraint(sp, x + y == D[node])

        @objective(sp, Min, x + 1000.0y)

        return sp
    end

    tree = narytree(3, 1)
    D = Dict(zip(collect(tree), D1))
    model = DetEqModel(
        tree,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )

    JuDGE.solve(model)

    tree2 = narytree(3, 2)
    D = Dict(zip(collect(tree2), D2))
    model2 = DetEqModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )
    JuDGE.solve(model2)
    obj2 = JuDGE.get_objval(model2)

    model3 = DetEqModel(
        tree2,
        ConditionallyUniformProbabilities,
        subproblem,
        JuDGE_DE_Solver,
        discount_factor = 0.5,
    )
    JuDGE.set_policy!(model3, model, :by_depth)
    JuDGE.solve(model3)
    obj3 = JuDGE.get_objval(model3)

    return obj3 - obj2
end
