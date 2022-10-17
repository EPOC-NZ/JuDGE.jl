using Random
using DelimitedFiles
using JuMP
using JuDGE
using Test

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function VRP(;
    seed::Int = 648,
    num_packages = 20,
    num_drivers = 4,
    max_deliveries = 6,
    use_heuristic = true,
    use_branch_callback = false,
    reltol = 0.001,
    rlx_reltol = 0.0001,
)
    Random.seed!(seed)
    package_values = rand(num_packages + 1)
    positions = rand(num_packages + 1, 2)
    package_values[1] = 0.0
    positions[1, 1] = 0.5
    positions[1, 2] = 0.5

    dist = Dict{Any,Float64}()

    for i in 1:num_packages+1
        for j in 1:num_packages+1
            dist[i, j] =
                (
                    (positions[i, 1] - positions[j, 1])^2 +
                    (positions[i, 2] - positions[j, 2])^2
                )^0.5
        end
    end
    dist

    mytree = narytree(0, 0)

    mytree.ext[:sum_max] = num_drivers
    mytree.ext[:sum_min] = 1

    num_packages += 1
    max_deliveries += 1

    function initialise_columns(judge::JuDGEModel)
        sequences = []
        package = 2
        for i in 1:num_drivers
            sequence = [1]
            for j in 1:(num_packages-1)÷num_drivers
                push!(sequence, package)
            end
            push!(sequences, sequence)
            package += 1
        end

        for s in sequences
            sp = sub_problems(judge.tree, deliveries = s, mipgap = 0.0)
            sp.ext[:form] = :binarycolumns
            sp.ext[:objective] = @variable(sp, obj)
            sp.ext[:objective_con] = @constraint(
                sp,
                sp.ext[:objective] - objective_function(sp) == 0
            )
            optimize!(sp)
            col = JuDGE.add_column(
                judge.master_problem,
                sp,
                judge.tree,
                branches = judge.ext[:branches],
            )
        end
    end

    function sub_problems(node; deliveries = nothing, mipgap = 0.00)
        sp = JuMP.Model(JuDGE_SP_Solver)
        @shutdown(sp, delivered[1:num_packages], Bin)

        @enforced(sp, 1 <= numdrivers <= 1, Int, lb = 1, ub = num_drivers)

        @capitalcosts(
            sp,
            -sum(
                package_values[i] * delivered[i] * 10 for i in 2:num_packages
            ) + 2 * numdrivers
        )

        @variable(sp, trip[1:num_packages, 1:num_packages], Bin)
        @variable(
            sp,
            0 <= counter[1:num_packages, 1:num_packages] <= max_deliveries,
            Int
        )
        @constraint(sp, delivered[1] <= 0)
        @constraint(sp, sum(trip[1, j] for j in 1:num_packages) == 1)
        @constraint(sp, sum(counter[1, j] for j in 1:num_packages) == 1)

        @constraint(sp, sum(trip[i, 1] for i in 1:num_packages) == 1)

        for i in 2:num_packages
            @constraint(
                sp,
                sum(trip[j, i] - trip[i, j] for j in 1:num_packages) == 0
            )
            @constraint(
                sp,
                sum(
                    counter[i, j] - trip[j, i] - counter[j, i] for
                    j in 1:num_packages
                ) == 0
            )
        end

        for i in 2:num_packages
            @constraint(
                sp,
                sum(trip[i, j] for j in 1:num_packages) >= delivered[i]
            )
        end
        @constraint(sp, sum(trip[i, i] for i in 1:num_packages) == 0)

        if typeof(deliveries) == Nothing
            @constraint(
                sp,
                sum(
                    trip[i, j] for i in 1:num_packages for j in 1:num_packages
                ) <= max_deliveries
            )
            @constraint(
                sp,
                sum(
                    counter[i, j] for i in 1:num_packages for
                    j in 1:num_packages
                ) <= max_deliveries * (max_deliveries + 1) / 2
            )
            for i in 1:num_packages
                for j in 1:num_packages
                    @constraint(
                        sp,
                        counter[i, j] <= max_deliveries * trip[i, j]
                    )
                end
            end

        else
            @constraint(
                sp,
                sum(
                    trip[i, j] for i in 1:num_packages for j in 1:num_packages
                ) <= length(deliveries)
            )
            @constraint(
                sp,
                sum(
                    counter[i, j] for i in 1:num_packages for
                    j in 1:num_packages
                ) <= length(deliveries) * (length(deliveries) + 1) / 2
            )
            for i in 1:num_packages
                for j in 1:num_packages
                    @constraint(
                        sp,
                        counter[i, j] <= length(deliveries) * trip[i, j]
                    )
                end
            end
            for i in 2:num_packages
                if i in deliveries
                    @constraint(sp, delivered[i] == 1)
                else
                    @constraint(sp, delivered[i] == 0)
                end
            end
        end

        @objective(
            sp,
            Min,
            sum(
                dist[i, j] * trip[i, j] for i in 1:num_packages for
                j in 1:num_packages
            )
        )

        return sp
    end

    function VRP_branch(jmodel::JuDGEModel, inttol::Float64)
        master = jmodel.master_problem
        subproblems = jmodel.sub_problems
        tree = jmodel.tree

        branches = JuDGE.variable_branch(jmodel, inttol)
        if branches !== nothing
            return branches
        end

        onebranch = zeros(
            length(subproblems[tree][:delivered]),
            length(subproblems[tree][:delivered]),
        )
        map(Main.eval, [:(jmodel = $jmodel)])
        for i in 1:length(subproblems[tree][:delivered])
            for j in i+1:length(subproblems[tree][:delivered])
                for col in master.ext[:columns][tree]
                    if i in keys(col.coeffs[:delivered]) &&
                       j in keys(col.coeffs[:delivered])
                        onebranch[i, j] += value(col.var)
                    end
                end
                onebranch[i, j] =
                    onebranch[i, j] > 1.0 - 10^-6 ? 0 : onebranch[i, j]
            end
        end

        index = 0
        #onebranch=abs.(onebranch-0.5*ones(length(subproblems[tree][:delivered]),length(subproblems[tree][:delivered])))

        # if minimum(onebranch)<0.4999
        # 	index=argmin(onebranch)
        # end

        if maximum(onebranch) > 0.1
            index = argmax(onebranch)
        end

        if index != 0
            println("Branching on: " * string(index))

            function filter1(col::JuDGE.Column)
                return (
                    (index[1] in keys(col.coeffs[:delivered])) ⊻
                    (index[2] in keys(col.coeffs[:delivered]))
                ) ? :ban : nothing
            end
            function filter2(col::JuDGE.Column)
                return (
                    (index[1] in keys(col.coeffs[:delivered])) &
                    (index[2] in keys(col.coeffs[:delivered]))
                ) ? :ban : nothing
            end

            constraint1 = JuDGE.BranchConstraint(
                @expression(
                    subproblems[tree],
                    subproblems[tree][:delivered][index[1]] -
                    subproblems[tree][:delivered][index[2]]
                ),
                :eq,
                0.0,
                tree,
            )
            constraint2 = JuDGE.BranchConstraint(
                @expression(
                    subproblems[tree],
                    subproblems[tree][:delivered][index[1]] +
                    subproblems[tree][:delivered][index[2]]
                ),
                :le,
                1.0,
                tree,
            )

            return [
                JuDGE.Branch(constraint1, filter1),
                JuDGE.Branch(constraint2, filter2),
            ]
        end
    end

    function bp_callback(
        termination::Termination,
        search::Symbol,
        log::Array{Array{JuDGE.ConvergenceState,1},1},
    )
        termination.rlx_reltol = 0.0002
        return termination, search
        # if (termination.allow_frac==:first_fractional && log[end][end].num_frac==0)
        # 	termination.allow_frac=:binary_solve_return_relaxation
        # 	return (termination,:lowestLB)
        # else
        # 	termination.allow_frac=:first_fractional
        #  	return (termination,:depth_first_dive)
        # end
    end

    function heuristic(incumbent::JuDGEModel)
        JuDGE.solve_binary(incumbent, 10^-9, 0.01, true, nothing)
        sequences = get_tours(incumbent, :nodes)
        JuDGE.remove_binary(incumbent)
        changed = true
        improve = 0.0
        while changed
            changed = false
            for n in 0:3
                for i in 1:length(sequences)-1
                    for j in i+1:length(sequences)
                        for k in 1:length(sequences[i])
                            skip = false
                            for p in 0:n
                                if sequences[i][(k+p-1)%length(
                                    sequences[i],
                                )+1] == 1
                                    skip = true
                                    break
                                end
                            end
                            if skip
                                continue
                            end
                            for l in 1:length(sequences[j])
                                skip = false
                                for p in 0:n
                                    if sequences[j][(l+p-1)%length(
                                        sequences[j],
                                    )+1] == 1
                                        skip = true
                                        break
                                    end
                                end
                                if skip
                                    continue
                                end
                                prev_k =
                                    (k - 2 + length(sequences[i])) %
                                    length(sequences[i]) + 1
                                k2 = (k + n - 1) % length(sequences[i]) + 1
                                next_k = (k + n) % length(sequences[i]) + 1
                                prev_l =
                                    (l - 2 + length(sequences[j])) %
                                    length(sequences[j]) + 1
                                l2 = (l + n - 1) % length(sequences[j]) + 1
                                next_l = (l + n) % length(sequences[j]) + 1
                                temp = 0.0
                                temp -=
                                    dist[sequences[i][k2], sequences[i][next_k]]
                                temp -=
                                    dist[sequences[i][prev_k], sequences[i][k]]
                                temp -=
                                    dist[sequences[j][l2], sequences[j][next_l]]
                                temp -=
                                    dist[sequences[j][prev_l], sequences[j][l]]

                                temp2 = 0.0
                                temp2 +=
                                    dist[sequences[i][k2], sequences[j][next_l]]
                                temp2 +=
                                    dist[sequences[j][l2], sequences[i][next_k]]
                                temp2 +=
                                    dist[sequences[i][prev_k], sequences[j][l]]
                                temp2 +=
                                    dist[sequences[j][prev_l], sequences[i][k]]

                                if temp + temp2 >= -10^-6
                                    temp2 = 0.0
                                    temp2 += dist[
                                        sequences[i][k2],
                                        sequences[j][prev_l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][l2],
                                        sequences[i][next_k],
                                    ]
                                    temp2 += dist[
                                        sequences[i][prev_k],
                                        sequences[j][l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][next_l],
                                        sequences[i][k],
                                    ]
                                end

                                if temp + temp2 >= -10^-6
                                    temp2 = 0.0
                                    temp2 += dist[
                                        sequences[i][k2],
                                        sequences[j][prev_l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][l2],
                                        sequences[i][prev_k],
                                    ]
                                    temp2 += dist[
                                        sequences[i][next_k],
                                        sequences[j][l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][next_l],
                                        sequences[i][k],
                                    ]
                                end

                                if temp + temp2 >= -10^-6
                                    temp2 = 0.0
                                    temp2 += dist[
                                        sequences[i][k2],
                                        sequences[j][next_l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][l2],
                                        sequences[i][prev_k],
                                    ]
                                    temp2 += dist[
                                        sequences[i][next_k],
                                        sequences[j][l],
                                    ]
                                    temp2 += dist[
                                        sequences[j][prev_l],
                                        sequences[i][k],
                                    ]
                                end

                                if temp + temp2 < -10^-6
                                    changed = true
                                    improve += temp + temp2
                                    for p in 0:n
                                        store = sequences[i][(k+p-1)%length(
                                            sequences[i],
                                        )+1]
                                        sequences[i][(k+p-1)%length(
                                            sequences[i],
                                        )+1] = sequences[j][(l+p-1)%length(
                                            sequences[j],
                                        )+1]
                                        sequences[j][(l+p-1)%length(
                                            sequences[j],
                                        )+1] = store
                                    end
                                    break
                                end
                            end
                            if changed
                                break
                            end
                        end
                        if changed
                            break
                        end
                    end
                    if changed
                        break
                    end
                end
                if changed
                    break
                end
            end
            if changed
                break
            end
        end
        if improve < 0.0
            println("Heuristic improved solution by " * string(improve))
            for s in sequences
                sp = sub_problems(incumbent.tree, deliveries = s, mipgap = 0.0)
                sp.ext[:form] = :binary
                sp.ext[:objective] = @variable(sp, obj)
                sp.ext[:objective_con] = @constraint(
                    sp,
                    sp.ext[:objective] - objective_function(sp) == 0
                )
                optimize!(sp)
                col = JuDGE.add_column(
                    incumbent.master_problem,
                    sp,
                    incumbent.tree,
                    branches = incumbent.ext[:branches],
                )
            end
        end
        return improve
    end

    function get_tours(jmodel::JuDGEModel, type::Symbol)
        active = JuDGE.get_active_columns(jmodel, inttol = 10^-5)

        sequences = []
        for a in active[mytree]
            deliveries = []
            for d in keys(a[1].coeffs[:delivered])
                push!(deliveries, d)
            end
            if 1 ∉ deliveries
                push!(deliveries, 1)
            end
            sp = sub_problems(mytree, deliveries = deliveries, mipgap = 0.0)
            optimize!(sp)
            if type == :arcs
                arcs = Tuple{Int,Int}[]
                count = length(deliveries)
                i = 1
                while count > 0
                    for j in 1:num_packages
                        if value(sp[:trip][i, j]) > 10^-7
                            push!(arcs, (i, j))
                            i = j
                            count -= 1
                            break
                        end
                    end
                end
                push!(sequences, arcs)
            elseif type == :nodes
                nodes = Int[]
                count = length(deliveries)
                i = 1
                while count > 0
                    for j in 1:num_packages
                        if value(sp[:trip][i, j]) > 10^-7
                            push!(nodes, i)
                            i = j
                            count -= 1
                            break
                        end
                    end
                end
                push!(sequences, nodes)
            end
        end
        return sequences
    end

    function show_vrp(jmodel::JuDGEModel)
        sequences = get_tours(jmodel, :arcs)

        nodetemp = "{id:LABEL,label:\"LABEL\",posX:x,posY:y},"
        temp1 = "var nodes = ["

        for i in 1:num_packages
            temp1 *= replace(
                replace(
                    replace(nodetemp, "LABEL" => string(i)),
                    "x" => positions[i, 1] * 800,
                ),
                "y" => positions[i, 2] * 800,
            )
        end
        temp1 *= "]"

        temp2 = "var edges = ["
        edgetemp = "{from:I,to:J,color:D},"
        d = 0
        colors = [
            "\'#ff0000\'",
            "\'#0000ff\'",
            "\'#00ee00\'",
            "\'#eeee00\'",
            "\'#ee00ee\'",
            "\'#00eeee\'",
            "\'#ee8800\'",
            "\'#88ee00\'",
            "\'#8800ee\'",
            "\'#0088ee\'",
        ]

        for s in sequences
            d += 1
            for e in s
                temp2 *= replace(
                    replace(replace(edgetemp, "I" => e[1]), "J" => e[2]),
                    "D" => colors[d],
                )
            end
        end
        temp2 *= "]"

        s = read(joinpath(dirname(@__DIR__), "visualise", "vrp.html"), String)
        filename = joinpath(
            "vrp" * string(Int(round((time() * 100) % 100000))) * ".html",
        )
        file = open(filename, "w")
        println(file, replace(s, "SOLUTION" => temp1 * "\n" * temp2))
        close(file)

        if Sys.iswindows()
            run(`$(ENV["COMSPEC"]) /c start $(filename)`)
        elseif Sys.isapple()
            run(`open $(filename)`)
        elseif Sys.islinux() || Sys.isbsd()
            run(`xdg-open $(filename)`)
        else
            error(
                "Unable to show plot. Try opening the file $(filename) manually.",
            )
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
    )
    initialise_columns(judy)

    h = nothing
    cb = nothing
    if use_heuristic
        h = heuristic
    end
    if use_branch_callback
        cb = bp_callback
    end
    judy = JuDGE.branch_and_price(
        judy,
        termination = Termination(
            reltol = reltol,
            rlx_reltol = rlx_reltol,
            inttol = 10^-6,
        ),
        bp_callback = cb,
        branch_method = VRP_branch,
        max_no_int = -10,
        heuristic = h,
    )
    show_vrp(judy)

    return JuDGE.get_objval(judy)
end

@time @test VRP(
    seed = 648,
    num_packages = 20,
    num_drivers = 3,
    max_deliveries = 8,
    use_heuristic = false,
) ≈ -103.43628 atol = 1e-3
@time @test VRP(
    seed = 648,
    num_packages = 20,
    num_drivers = 3,
    max_deliveries = 8,
    use_heuristic = true,
) ≈ -103.43628 atol = 1e-3
@time @test VRP(
    seed = 648,
    num_packages = 20,
    num_drivers = 3,
    max_deliveries = 8,
    use_heuristic = true,
    use_branch_callback = true,
    rlx_reltol = 0.15,
) ≈ -103.43628 atol = 1e-4

@test VRP(
    seed = 120,
    num_packages = 30,
    num_drivers = 4,
    max_deliveries = 8,
    use_heuristic = true,
    reltol = 0.005,
    rlx_reltol = 0.002,
) ≈ -130.8 atol = 1e-1
