using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function deterministic(;
    formulation = :decomp,
    var_type = :expansion,
    discount = 1.0,
)
    depth = 2
    mytree = narytree(depth, 1)

    demand = Dict(zip(collect(mytree), [2, 1, 3]))

    function sub_problems(n)
        sp = JuMP.Model(JuDGE_SP_Solver)

        if var_type == :expansion
            @expansion(sp, big, Bin)
            #@expansion(sp,0<=small<=2,Int)
            @expansion(sp, small[1:2], Bin)
        elseif var_type == :enforced
            @enforced(sp, big, Bin)
            #@enforced(sp,0<=small<=2,Int)
            @enforced(sp, small[1:2], Bin)
        end

        @capitalcosts(sp, 16 * big + 10 * sum(small))

        @variable(sp, supply >= 0)
        @variable(sp, shortage >= 0)

        @constraint(sp, supply <= 1 + 2 * big + sum(small))
        @constraint(sp, demand[n] <= supply + shortage)

        @objective(sp, Min, 2 * big + sum(small) + 5 * supply + 20 * shortage)

        return sp
    end

    model = nothing
    if formulation == :decomp
        model = JuDGEModel(
            mytree,
            UniformLeafProbabilities,
            sub_problems,
            JuDGE_MP_Solver,
            discount_factor = discount,
        )
        model = JuDGE.branch_and_price(model, verbose = 0)
        println("\nRe-solved Objective: " * string(resolve_subproblems(model)))
    elseif formulation == :deteq
        model = DetEqModel(
            mytree,
            UniformLeafProbabilities,
            sub_problems,
            JuDGE_DE_Solver,
            discount_factor = discount,
        )
        JuDGE.solve(model)
        println("\nObjective: " * string(JuDGE.get_objval(model)))
    else
        @error("Invalid formulation type")
    end
    JuDGE.print_expansions(model, onlynonzero = true, inttol = 10^-5)

    return JuDGE.get_objval(model)
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    deterministic(formulation = :decomp, var_type = :expansion, discount = 0.85)
    deterministic(formulation = :decomp, var_type = :enforced, discount = 0.85)
end
