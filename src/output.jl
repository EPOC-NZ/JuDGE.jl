"""
	print_expansions(jmodel::JuDGEModel;
                     onlynonzero::Bool=true,
                     inttol=10^-9,
                     format=nothing)

Given a solved JuDGE model, this function will write the optimal capacity expansion
decisions to the REPL.

### Required Arguments
`jmodel` is the JuDGE model whose solution we wish to write to a file

### Optional Arguments
`onlynonzero` is a boolean, if set to `true` the function will only print expansions
with a non-zero value.

`inttol` is the integrality tolerance; any expansion variable value less than this
will be treated as 0, and any value greater than 1-`inttol` will be treated as 1

`format` is a function that specifies customised printing of expansion values.
See [Tutorial 2: Formatting output](@ref) for more details.
"""
function print_expansions(
    jmodel::JuDGEModel;
    onlynonzero::Bool = true,
    inttol = 10^-9,
    format = nothing,
)
    function process(x, val)
        if jmodel.sub_problems[jmodel.tree].ext[:options][x][4] == :Con
            return string(val)
        elseif jmodel.sub_problems[jmodel.tree].ext[:options][x][4] == :Bin
            val = val > 1.0 - inttol ? 1.0 : val
            val = val < inttol ? 0.0 : val
            return string(val)
        else
            val = val > ceil(val) - inttol ? ceil(val) : val
            val = val < floor(val) + inttol ? floor(val) : val
            return string(val)
        end
    end

    if termination_status(jmodel.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel.master_problem) != MOI.LOCALLY_SOLVED
        error("You need to first solve the decomposed model.")
    end

    println("\nJuDGE Expansions")

    for node in collect(jmodel.tree)
        for (x, var) in jmodel.master_problem.ext[:expansions][node]
            if typeof(var) <: AbstractArray
                if typeof(format) <: Function
                    if format_output(
                        node,
                        x,
                        format(x, JuMP.value.(var)),
                        onlynonzero,
                        inttol,
                    )
                        continue
                    end
                end
                val = JuMP.value.(var)
                if typeof(var) <: JuMP.Containers.SparseAxisArray
                    val = val.data
                end
                for key in keys(val)
                    if !onlynonzero || val[key] > inttol
                        if typeof(val) <: Array
                            strkey = string(key)
                            strkey = replace(strkey, "CartesianIndex(" => "")
                            strkey = replace(strkey, ")" => "")
                            strkey = replace(strkey, ", " => ",")
                            temp =
                                "Node " *
                                node.name *
                                ": \"" *
                                string(x) *
                                "[" *
                                strkey *
                                "]\" " *
                                process(x, val[key])
                        elseif typeof(val) <: Dict
                            strkey = string(key)
                            strkey = replace(strkey, ")" => "")
                            strkey = replace(strkey, "(" => "")
                            strkey = replace(strkey, ", " => ",")
                            temp =
                                "Node " *
                                node.name *
                                ": \"" *
                                string(x) *
                                "[" *
                                strkey *
                                "]\" " *
                                process(x, val[key])
                        else
                            typeof(val) <: JuMP.Containers.DenseAxisArray
                            temp =
                                "Node " * node.name * ": \"" * string(x) * "["

                            for i in 1:length(val.axes)-1
                                temp *= string(key[i]) * ","
                            end
                            temp *=
                                string(key[length(val.axes)]) *
                                "]\" " *
                                process(x, val[key])
                        end
                        println(temp)
                    end
                end
            else
                if typeof(format) <: Function
                    if format_output(
                        node,
                        x,
                        format(x, JuMP.value(var)),
                        onlynonzero,
                        inttol,
                    )
                        continue
                    end
                end
                if !onlynonzero || JuMP.value(var) > inttol
                    println(
                        "Node " *
                        node.name *
                        ": \"" *
                        string(x) *
                        "\" " *
                        process(x, JuMP.value(var)),
                    )
                end
            end
        end
    end
end

"""
	print_expansions(deteq::DetEqModel;
                     onlynonzero::Bool=true,
                     inttol=10^-9,
                     format=nothing)

Given a solved deterministic equivalent model, this function will write the optimal
capacity expansion decisions to the REPL.

### Required Arguments
`deteq` is the deterministic equivalent model whose solution we wish to write to a file

### Optional Arguments
`onlynonzero` is a boolean, if set to `true` the function will only print expansions
with a non-zero value.

`inttol` is the integrality tolerance; any expansion variable value less than this
will be treated as 0, and any value greater than 1-`inttol` will be treated as 1

`format` is a function that specifies customised printing of expansion values.
See [Tutorial 2: Formatting output](@ref) for more details.
"""
function print_expansions(
    deteq::DetEqModel;
    onlynonzero::Bool = true,
    inttol = 10^-9,
    format = nothing,
)
    function process(var)
        val = JuMP.value(var)
        if is_binary(var)
            return string(val > 1 - inttol ? 1.0 : val)
        else
            return string(val)
        end
    end

    if termination_status(deteq.problem) != MOI.OPTIMAL &&
       termination_status(deteq.problem) != MOI.TIME_LIMIT &&
       termination_status(deteq.problem) != MOI.INTERRUPTED
        error("You need to first solve the deterministic equivalent model.")
    end

    println("\nDeterministic Equivalent Expansions")
    for node in collect(deteq.tree)
        for x in keys(deteq.problem.ext[:master_vars][node])
            var = deteq.problem.ext[:master_vars][node][x]
            if typeof(var) == VariableRef
                if typeof(format) <: Function
                    if format_output(
                        node,
                        x,
                        format(x, JuMP.value(var)),
                        onlynonzero,
                        inttol,
                    )
                        continue
                    end
                end
                if !onlynonzero || JuMP.value(var) > inttol
                    name = deteq.problem.ext[:master_names][node][x]
                    println(
                        "Node " *
                        node.name *
                        ": \"" *
                        name *
                        "\" " *
                        process(var),
                    )
                end
            elseif typeof(var) == Dict{Any,Any}
                if typeof(format) <: Function
                    val = Dict{Any,Float64}()
                    for key in collect(keys(var))
                        val[key] = JuMP.value(var[key])
                    end
                    if format_output(
                        node,
                        x,
                        format(x, val),
                        onlynonzero,
                        inttol,
                    )
                        continue
                    end
                end

                for i in eachindex(var)
                    if !onlynonzero || JuMP.value(var[i]) > inttol
                        name = deteq.problem.ext[:master_names][node][x][i]
                        println(
                            "Node " *
                            node.name *
                            ": \"" *
                            name *
                            "\" " *
                            process(var[i]),
                        )
                    end
                end
            end
        end
    end
end

function format_output(node::AbstractTree, x::Symbol, exps, onlynonzero, inttol)
    if typeof(exps) == Float64 || typeof(exps) == Int
        if !onlynonzero || abs(exps) > inttol
            println(
                "Node " * node.name * ": \"" * string(x) * "\" " * string(exps),
            )
        end
        return true
    elseif typeof(exps) == Dict{AbstractArray,Float64} ||
           typeof(exps) == Dict{AbstractArray,Int}
        for (key, exp) in exps
            if !onlynonzero || abs(exp) > inttol
                println(
                    "Node " *
                    node.name *
                    ": \"" *
                    string(x) *
                    string(key) *
                    "\" " *
                    string(exp),
                )
            end
        end
        return true
    elseif typeof(exps) == Dict{Tuple,Float64} ||
           typeof(exps) == Dict{Tuple,Int}
        for (key, exp) in exps
            if !onlynonzero || abs(exp) > inttol
                s_key = string(key)
                s_key = replace(s_key, "(" => "")
                s_key = replace(s_key, ")" => "")
                s_key = replace(s_key, "\"" => "")
                s_key = replace(s_key, ", " => ",")
                println(
                    "Node " *
                    node.name *
                    ": \"" *
                    string(x) *
                    "[" *
                    s_key *
                    "]\" " *
                    string(exp),
                )
            end
        end
        return true
    elseif typeof(exps) == Dict{Int,Float64} ||
           typeof(exps) == Dict{Symbol,Float64} ||
           typeof(exps) == Dict{String,Float64}
        for (key, exp) in exps
            if !onlynonzero || abs(exp) > inttol
                println(
                    "Node " *
                    node.name *
                    ": \"" *
                    string(x) *
                    "[" *
                    string(key) *
                    "]\" " *
                    string(exp),
                )
            end
        end
        return true
    elseif typeof(exps) != Nothing
        error(
            "Formatting function must return a Float64, Dict{AbstractArray,Float64}, or nothing.",
        )
    end
    return false
end

"""
	write_solution_to_file(model::Union{JuDGEModel,DetEqModel},filename::String)

Given a deterministic equivalent model and a filename, this function writes the
entire solution to a CSV.

### Required Arguments
`model` can be either the `JuDGEModel` or the `DetEqModel` whose solution we wish to write to a file

`filename` is the output filename
"""
function write_solution_to_file(
    model::Union{JuDGEModel,DetEqModel},
    filename::String,
)
    file = open(filename, "w")

    println(file, "node,value,variable,indices")

    solution = solution_to_dictionary(model)

    for node in keys(solution)
        for (name, item) in solution[node]
            if typeof(item) == Float64
                println(
                    file,
                    node.name * "," * string(item) * "," * string(name) * ",",
                )
            else
                for (index, value) in item
                    println(
                        file,
                        node.name *
                        "," *
                        string(value) *
                        "," *
                        string(name) *
                        "," *
                        index,
                    )
                end
            end
        end
    end

    return close(file)
end
