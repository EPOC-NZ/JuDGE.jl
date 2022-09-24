using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function test_variable_types(
    tree::AbstractTree,
    demand::Dict{AbstractTree,<:Real},
)
    function sub_problems(n::AbstractTree)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @expansion(sp, e >= 0, Bin)
        @expansion(sp, a[1:4] >= 0, Bin)
        @expansion(sp, 0 <= b[2:5] <= 10, Int)
        @expansion(sp, 0 <= c[[:one, :two]] <= 3)
        @expansion(sp, 0 <= d[i in 1:4, 1:i] <= 2)

        @variable(sp, x[1:4] >= 0)
        @variable(sp, y[1:4] >= 0)
        @variable(sp, z[1:4] >= 0)

        @capitalcosts(
            sp,
            sum(a) +
            3 * sum(b) +
            sum(c) +
            sum(d[i, j] * j for i in 1:4 for j in 1:i) +
            3 * e
        )

        @constraint(sp, meet_demand[i in 1:4], x[i] + y[i] + z[i] == demand[n])

        @constraint(
            sp,
            available_x[i in 1:4],
            x[i] <= a[i] + b[i+1] + c[:one] + e
        )
        @constraint(
            sp,
            available_y[i in 1:4],
            y[i] <= sum(d[i, j] for j in 1:i) + c[:two] + e
        )

        @objective(sp, Min, sum(y) + 10 * sum(z))

        return sp
    end

    model = JuDGEModel(
        tree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        discount_factor = 0.9,
    )

    JuDGE.solve(model)
    resolve_subproblems(model)

    model2 = DetEqModel(
        tree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        discount_factor = 0.9,
    )
    JuDGE.set_starting_solution!(model2, model)

    JuDGE.solve(model2)

    JuDGE.write_solution_to_file(
        model,
        joinpath(@__DIR__, "variable_types_decomp.csv"),
    )
    JuDGE.write_solution_to_file(
        model2,
        joinpath(@__DIR__, "variable_types_deteq.csv"),
    )

    return JuDGE.get_objval(model), JuDGE.get_objval(model2)
end
