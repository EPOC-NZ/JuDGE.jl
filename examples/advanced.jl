using Random, Test
using Gurobi, JuMP, JuDGE

include("solvers/setup_gurobi.jl")

function knapsack_advanced(
    seed::Int,
    tree_size::Tuple{Int,Int},
    numitems::Int,
    custom_solve_settings::Symbol,
)
    Random.seed!(seed)
    # how many investments?
    numinvest = 6

    # size of tree?
    height, degree = tree_size

    totalnodes = Int((degree^(height + 1) - 1) / (degree - 1))

    investcost = zeros(totalnodes, numinvest)
    for i in 1:totalnodes
        investcost[i, :] =
            ([1, 1.8, 3.5, 6.8, 13.5, 25]) *
            (1 - ((i - 1) / (totalnodes * 1.2)))
    end

    investvol = [1, 2, 4, 8, 16, 32]
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
        print("Formulating subproblem for node ")
        println(node.name)
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
        set_optimizer_attribute(model, "MIPGap", 0.0005)
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

    if custom_solve_settings == :callbacks
        function oa(judge::JuDGEModel, stalled::Bool, initialize::Bool)
            function slowed()
                if judge.log[end].iter > 3
                    if judge.log[end].rlx_abs / judge.log[end-3].rlx_abs > 0.8
                        return true
                    end
                end
                return false
            end

            function relgap_factor(probability)
                if judge.ext[:optimizer_settings][:level] == 4
                    0.0001
                else
                    10^-4 / probability *
                    20^(4 - judge.ext[:optimizer_settings][:level])
                end
            end

            if initialize
                if !haskey(judge.ext[:optimizer_settings], :level)
                    judge.ext[:optimizer_settings][:level] = 1
                end

                println(
                    "\u1b[1F\u1b[130GInitializing optimizer settings at level " *
                    string(judge.ext[:optimizer_settings][:level]),
                )

                for (node, sp) in judge.sub_problems
                    set_optimizer_attribute(
                        sp,
                        "MIPGap",
                        relgap_factor(judge.probabilities[node]),
                    )
                end
                return false
            elseif judge.ext[:optimizer_settings][:level] == 4
                if stalled
                    return true
                else
                    return false
                end
            elseif stalled || slowed()
                judge.ext[:optimizer_settings][:level] += 1

                for (node, sp) in judge.sub_problems
                    set_optimizer_attribute(
                        sp,
                        "MIPGap",
                        relgap_factor(judge.probabilities[node]),
                    )
                end

                if judge.ext[:optimizer_settings][:level] == 4
                    println("\u1b[1F\u1b[130GFinal optimizer settings")
                else
                    println(
                        "\u1b[1F\u1b[130GIncremented optimizer settings to level " *
                        string(judge.ext[:optimizer_settings][:level]),
                    )
                end
            end

            return false
        end

        function mp_cb(judge::JuDGEModel, abs_tol::Float64, rel_tol::Float64)
            function earlytermination(cb_data, cb_where)
                if cb_where == GRB_CB_MIP
                    objbst = Ref{Cint}()
                    GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, objbst)
                    objbnd = Ref{Cint}()
                    GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBND, objbnd)
                    runtime = Ref{Cint}()
                    GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
                    JuDGE.printright(
                        "Incumbent: " *
                        string(objbst[]) *
                        " Bound: " *
                        string(objbnd[]) *
                        " Runtime: " *
                        string(runtime[]),
                    )
                    if (
                           objbnd[] > judge.bounds.LB + abs_tol &&
                           objbnd[] >
                           judge.bounds.LB + rel_tol * abs(judge.bounds.LB)
                       ) ||
                       objbst[] < judge.bounds.LB + abs_tol ||
                       objbst[] <
                       judge.bounds.LB + rel_tol * abs(judge.bounds.LB)
                        GRBterminate(
                            backend(judge.master_problem).optimizer.model.inner,
                        )
                    end
                end
                return
            end
            return MOI.set(
                judge.master_problem,
                Gurobi.CallbackFunction(),
                earlytermination,
            )
        end

        judy = JuDGE.branch_and_price(
            judy,
            termination = Termination(
                rlx_reltol = 10^-6,
                reltol = 5 * 10^-3,
                inttol = 10^-6,
            ),
            optimizer_attributes = oa,
            mp_callback = nothing,
            max_no_int = 5,
        )
    elseif custom_solve_settings == :partialpricing
        blocks = JuDGE.get_groups(mytree, combine = height - 1)
        judy = JuDGE.branch_and_price(
            judy,
            termination = Termination(
                rlx_reltol = 10^-6,
                reltol = 5 * 10^-3,
                inttol = 10^-6,
            ),
            max_no_int = 5,
            blocks = blocks,
        )
    else
        judy = JuDGE.branch_and_price(
            judy,
            termination = Termination(
                rlx_reltol = 10^-6,
                reltol = 5 * 10^-3,
                inttol = 10^-6,
            ),
            max_no_int = 5,
        )
    end

    println("Objective: " * string(objective_value(judy.master_problem)))
    JuDGE.print_expansions(judy, format = format_output, inttol = 10^-6)

    return objective_value(judy.master_problem)
end

# Solve without callbacks for medium tree and 100 items
@time @test knapsack_advanced(50, (4, 3), 100, :standard) ≈ -14.05 atol = 2e-2

# Solve with callbacks for medium tree and 100 items
@time @test knapsack_advanced(50, (4, 3), 100, :callbacks) ≈ -14.05 atol = 2e-2

# Solve with parial pricing for medium tree and 100 items
@time @test knapsack_advanced(50, (4, 3), 100, :partialpricing) ≈ -14.05 atol =
    2e-2
