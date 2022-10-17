using Random
using JuMP
using JuDGE
using Test

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function cutting_stock(;
    seed::Int = 200,
    L = 10,
    sizes = [2, 3, 5],
    test = false,
)
    mytree = narytree(0, 0)

    n = length(sizes)

    if test
        demand = [26, 31, 109]
    else
        Random.seed!(seed)
        demand = rand(10:200, n)
    end

    num_sizes = length(sizes)
    max_repeat = L ÷ minimum(sizes)
    max_patterns = 2

    mytree.ext[:sum_max] = sum(ceil(demand[i] / (L ÷ sizes[i])) for i in 1:n)#max_patterns
    mytree.ext[:sum_min] = 1
    mytree.ext[:col_max] =
        maximum(ceil(demand[i] / (L ÷ sizes[i])) for i in 1:n)

    function initialise_columns(judge::JuDGEModel)
        patterns = []
        for i in 1:n
            z = zeros(n)
            z[i] = L ÷ sizes[i]
            push!(patterns, z)
        end

        for p in patterns
            sp = sub_problems(judge.tree, pattern = p, mipgap = 0.0)
            sp.ext[:form] = :binary
            sp.ext[:objective] = @variable(sp, obj)
            sp.ext[:objective_con] = @constraint(
                sp,
                sp.ext[:objective] - objective_function(sp) == 0
            )
            sp.ext[:discrete] = VariableRef[]
            sp.ext[:discrete_branch] = VariableRef[]
            optimize!(sp)
            col = JuDGE.add_column(judge.master_problem, sp, judge.tree)
        end
    end

    function sub_problems(node; pattern = nothing, mipgap = 0.00)
        sp = JuMP.Model(JuDGE_SP_Solver)
        max_branch = L * 4
        @shutdown(
            sp,
            0 <= made[1:num_sizes] <= L ÷ minimum(sizes),
            Int,
            lb = 0.0,
            ub = 1000
        )

        @variable(sp, dif[1:num_sizes, 1:max_branch] >= 0)
        @variable(sp, val[1:num_sizes, 1:max_branch])
        @variable(sp, switch[1:num_sizes, 1:max_branch], Bin)

        @constraint(
            sp,
            dif1[i in 1:num_sizes, j in 1:max_branch],
            dif[i, j] <= switch[i, j]
        )
        @constraint(
            sp,
            dif2[i in 1:num_sizes, j in 1:max_branch],
            dif[i, j] <=
            made[i] - val[i, j] + (L ÷ minimum(sizes)) * (1 - switch[i, j])
        )

        @constraint(
            sp,
            block[j in 1:max_branch],
            sum(dif[i, j] for i in 1:num_sizes) >= 1
        )

        @constraint(sp, sum(sizes[i] * made[i] for i in 1:num_sizes) <= L)

        if typeof(pattern) != Nothing
            @constraint(sp, set[i = 1:num_sizes], made[i] == pattern[i])
        end

        @objective(sp, Min, 1.0)
        return sp
    end

    function meet_demand(model, tree)
        @constraint(
            model,
            meet_demand[i in 1:num_sizes],
            made[tree][i] >= demand[i]
        )
    end

    function pattern_branch(jmodel::JuDGEModel, inttol::Float64)
        master = jmodel.master_problem
        subproblems = jmodel.sub_problems
        tree = jmodel.tree

        index = nothing
        maxvalue = inttol
        for col in master.ext[:columns][tree]
            val = value(col.var) + 0.5
            v = val > ceil(val) - inttol ? 0 : val - floor(val)
            if v > maxvalue
                maxvalue = v
                index = col
            end
        end

        if index !== nothing
            temp = index.coeffs[:made]

            branch1 = JuDGE.BranchConstraint(
                index.var,
                :ge,
                ceil(value(index.var)),
                master,
            )

            branch2 = JuDGE.BranchConstraint[]
            push!(
                branch2,
                JuDGE.BranchConstraint(
                    index.var,
                    :le,
                    floor(value(index.var)),
                    master,
                ),
            )
            j = length(jmodel.ext[:branches]) + 1
            for i in 1:num_sizes
                if haskey(temp, i)
                    push!(
                        branch2,
                        JuDGE.BranchConstraint(
                            subproblems[tree][:val][i, j],
                            :eq,
                            temp[i],
                            tree,
                        ),
                    )
                else
                    push!(
                        branch2,
                        JuDGE.BranchConstraint(
                            subproblems[tree][:val][i, j],
                            :eq,
                            0.0,
                            tree,
                        ),
                    )
                end
            end
            f = function f(col::JuDGE.Column)
                at_least_one_less = false

                for i in 1:num_sizes
                    if !haskey(col.coeffs[:made], i) && !haskey(temp, i)
                        continue
                    end
                    if !haskey(temp, i) || (
                        haskey(col.coeffs[:made], i) &&
                        col.coeffs[:made][i] > temp[i]
                    )
                        return
                    elseif !haskey(col.coeffs, i) || (
                        haskey(temp, i) && col.coeffs[:made][i] < temp[i]
                    )
                        at_least_one_less = true
                    end
                end
                if at_least_one_less
                    return :ban
                else
                    return
                end
            end
            println(f(index))
            return [JuDGE.Branch(branch1), JuDGE.Branch(branch2, f)]
        end
    end

    function print_patterns(jmodel::JuDGEModel)
        active = JuDGE.get_active_columns(judy, inttol = 10^-6)
        println("")
        print("Demands ")
        for i in 1:n
            print(" ")
            print(sizes[i])
            print(": ")
            print(demand[i])
            if i != n
                print(", ")
            end
        end
        println("\n")
        for (col, count) in active[mytree]
            print("[")
            for (index, number) in col.coeffs[:made]
                for n in 1:number
                    print(" ")
                    print(sizes[index])
                    print(" ")
                end
            end
            print("]: ")
            println(count)
        end
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        sideconstraints = meet_demand,
        check = false,
    )
    initialise_columns(judy)

    judy = JuDGE.branch_and_price(
        judy,
        termination = Termination(abstol = 0.99),
        branch_method = pattern_branch,
        verbose = 1,
    )
    print_patterns(judy)

    return JuDGE.get_objval(judy)
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    cutting_stock() # 70
    cutting_stock(
        seed = 12,
        L = 100,
        sizes = [2, 3, 5, 7, 11, 13, 17, 19, 24, 37, 47, 62],
    ) # 238
end
