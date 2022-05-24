using DelimitedFiles
using JuMP
using JuDGE
using Dash
using DashWrapper #(] add https://github.com/adow031/DashWrapper.git)

include("solvers/setup_gurobi.jl")

mytree = narytree(5, 2)

function invest_supply_cost(node)
    if node.parent == nothing
        return Dict(zip(supply_nodes, [1.0, 2.0]))
    else
        p = node.parent
        for i in 1:length(p.children)
            if p.children[i] == node
                temp = deepcopy(invest_supply_cost(p))
                for key in keys(temp)
                    temp[key] *= (0.2 * i + 0.8)
                end
                return temp
            end
        end
    end
end

function invest_arc_cost(node)
    if node.parent == nothing
        temp = []
        for i in supply_nodes
            for j in demand_nodes
                push!(temp, (i, j))
            end
        end
        return Dict(zip(temp, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
    else
        p = node.parent
        for i in 1:length(p.children)
            if p.children[i] == node
                temp = deepcopy(invest_arc_cost(p))
                for key in keys(temp)
                    temp[key] *= (0.2 * i + 0.8)
                end
                return temp
            end
        end
    end
end

function demand(node)
    if node.parent == nothing
        return d_dict
    else
        p = node.parent
        for i in 1:length(p.children)
            if p.children[i] == node
                temp = deepcopy(demand(p))
                for key in keys(temp)
                    temp[key] += i * i
                end
                return temp
            end
        end
    end
end

function supply(node)
    return s_dict
end

data_file = "transportation.csv"
data = readdlm(joinpath(@__DIR__, data_file), ',')

supply_nodes = data[3:end, 2]
s = data[3:end, 1]

demand_nodes = collect(data[2, 3:end])
d = collect(data[1, 3:end])

c = data[3:end, 3:end]

# Converting arrays to dictionaries
s_dict = Dict(zip(supply_nodes, s * 2))
d_dict = Dict(zip(demand_nodes, d))

c_dict = Dict()
for i in 1:length(supply_nodes)
    for j in 1:length(demand_nodes)
        c_dict[supply_nodes[i], demand_nodes[j]] = c[i, j]
    end
end

### with judge
function sub_problems(node)
    model = Model(JuDGE_SP_Solver)

    @expansion(model, new_supply[supply_nodes], Bin) #invest in more supply
    @expansion(model, 0 <= new_capacity[supply_nodes, demand_nodes] <= 2, Int) #invest in more arc capacity

    @capitalcosts(
        model,
        sum(invest_supply_cost(node)[i] * new_supply[i] for i in supply_nodes) +
        sum(
            invest_arc_cost(node)[i, j] * new_capacity[i, j] for
            i in supply_nodes for j in demand_nodes
        )
    )

    @variable(model, x[supply_nodes, demand_nodes] >= 0)

    @objective(
        model,
        Min,
        sum(c_dict[i, j] * x[i, j] for i in supply_nodes, j in demand_nodes)
    )
    @constraint(
        model,
        SupplyIncrease[i in supply_nodes],
        sum(x[i, j] for j in demand_nodes) <=
        supply(node)[i] + s_dict[i] * new_supply[i]
    )

    @constraint(
        model,
        CapacityIncrease[i in supply_nodes, j in demand_nodes],
        x[i, j] <= c_dict[i, j] + c_dict[i, j] * new_capacity[i, j]
    )

    @constraint(
        model,
        demandCon[j in demand_nodes],
        sum(x[i, j] for i in supply_nodes) == demand(node)[j]
    )

    return model
end

function format_output(s::Symbol, values)
    if s == :new_capacity
        output = Dict{Tuple,Float64}()
        for i in supply_nodes
            for j in demand_nodes
                output[i, j] = values[i, j] * c_dict[i, j]
            end
        end
        return output
    elseif s == :new_supply
        output = Dict{String,Float64}()
        for i in supply_nodes
            output[i] = values[i] * s_dict[i]
        end
        return output
    end
    return nothing
end

judy = JuDGEModel(
    mytree,
    ConditionallyUniformProbabilities,
    sub_problems,
    JuDGE_MP_Solver,
    discount_factor = 0.9,
)
JuDGE.solve(judy, termination = Termination(inttol = 10^-7), skip_nodes = 2)

println("\nObjective: " * string(objective_value(judy.master_problem)) * "\n")
JuDGE.print_expansions(judy, format = format_output)

println("\nRe-solved Objective: " * string(resolve_subproblems(judy)))

solution = JuDGE.solution_to_dictionary(judy)

function simple_dash()
    # Create the Dash app without any extra blocks or callbacks specified.
    # The path is the folder that the assets directory will be copied into.
    blocks = JuDGE.dash_layout()
    DashWrapper.scale_blocks!(blocks, 800, 1000, 16)

    treedata = JuDGE.export_tree(mytree, data = solution)

    return JuDGE.create_dash(blocks, nothing, path = @__DIR__, data = treedata)
end

function medium_dash()
    blocks = JuDGE.dash_layout()
    push!(blocks, DashWrapper.Block((1, 0), (0.5, 6), "Solution at Node"))
    DashWrapper.scale_blocks!(blocks, 1200, 1000, 16)

    treedata = JuDGE.export_tree(mytree, data = solution)

    # Create the Dash app with the blocks defined above, but no callbacks.
    # The extra file has a javascript function that defines an action to
    # take when a node is selected in the scenario tree.
    return JuDGE.create_dash(
        blocks,
        nothing,
        path = @__DIR__,
        extra_files = [joinpath(@__DIR__, "plot_functions", "medium.js")],
        data = treedata,
    )
end

function advanced_dash()
    # Define the layout of the dashbaord.
    blocks = JuDGE.dash_layout()
    push!(
        blocks,
        DashWrapper.Block((1, 0), (1.2, 3), "Transportation Capacity"),
    )
    push!(blocks, DashWrapper.Block((1, 3), (0.5, 3), "Supply Capacity"))
    push!(blocks, DashWrapper.Block((1.5, 3), (0.7, 3), "Flow"))
    push!(blocks, DashWrapper.Block((2.2, 0), (0.5, 6), "Solution at Node"))
    DashWrapper.scale_blocks!(blocks, 1800, 800, 16)

    treedata = JuDGE.export_tree(mytree, data = solution)

    # Define functions to format the solution data for the plots.
    function capacity_data(node::AbstractTree, type::Symbol, solution::Dict)
        x = []
        y = []
        if type == :existing
            for (n_c, c) in c_dict
                push!(x, n_c[1] * "," * n_c[2])
                push!(y, c)
            end
        elseif type == :built
            for (n_c, c) in solution[node][:new_capacity]
                push!(x, string(n_c))
                push!(y, c - solution[node][:new_capacity_master][n_c])
            end
        else
            type == :justbuilt
            for (n_c, c) in solution[node][:new_capacity_master]
                push!(x, string(n_c))
                push!(y, c)
            end
        end

        return (x = x, y = y)
    end

    function supply_data(node::AbstractTree, type::Symbol, solution::Dict)
        x = []
        y = []
        if type == :existing
            for (n_s, s) in s_dict
                push!(x, n_s)
                push!(y, s)
            end
        elseif type == :built
            for (n_s, s) in solution[node][:new_supply]
                push!(x, string(n_s))
                push!(
                    y,
                    s * s_dict[n_s] -
                    solution[node][:new_supply_master][n_s] * s_dict[n_s],
                )
            end
        else
            type == :justbuilt
            for (n_s, s) in solution[node][:new_supply_master]
                push!(x, string(n_s))
                push!(y, s * s_dict[n_s])
            end
        end

        return (x = x, y = y)
    end

    # Define callback function that creates three plots.
    # These plots depend on the 'node' that is selected.
    function cb_func(arg)
        if arg[1] == 0
            return "", "", ""
        else
            node = collect(mytree)[arg[1]]
        end
        return dcc_graph(
            id = "example-graph-1",
            figure = (
                data = [
                    DashWrapper.join_tuples(
                        capacity_data(node, :existing, solution),
                        (type = "bar", name = "Existing Capacity"),
                    ),
                    DashWrapper.join_tuples(
                        capacity_data(node, :built, solution),
                        (type = "bar", name = "Built Capacity"),
                    ),
                    DashWrapper.join_tuples(
                        capacity_data(node, :justbuilt, solution),
                        (type = "bar", name = "New Capacity"),
                    ),
                ],
                layout = DashWrapper.join_tuples(
                    (
                        barmode = "stack",
                        margin = (l = 50, r = 50, b = 60, t = 25),
                    ),
                    NamedTuple{(:width, :height)}(blocks[3].size),
                ),
            ),
        ),
        dcc_graph(
            id = "example-graph-2",
            figure = (
                data = [
                    DashWrapper.join_tuples(
                        supply_data(node, :existing, solution),
                        (type = "bar", name = "Existing Supply"),
                    ),
                    DashWrapper.join_tuples(
                        supply_data(node, :built, solution),
                        (type = "bar", name = "Built Supply"),
                    ),
                    DashWrapper.join_tuples(
                        supply_data(node, :justbuilt, solution),
                        (type = "bar", name = "New Supply"),
                    ),
                ],
                layout = DashWrapper.join_tuples(
                    (
                        barmode = "stack",
                        margin = (l = 50, r = 50, b = 60, t = 25),
                    ),
                    NamedTuple{(:width, :height)}(blocks[4].size),
                ),
            ),
        ),
        dcc_graph(
            id = "example-graph-3",
            figure = (
                data = [(
                    x = collect(keys(solution[node][:x])),
                    y = collect(values(solution[node][:x])),
                    type = "bar",
                    name = "Flow",
                )],
                layout = DashWrapper.join_tuples(
                    (
                        barmode = "stack",
                        margin = (l = 50, r = 60, b = 100, t = 25),
                    ),
                    NamedTuple{(:width, :height)}(blocks[5].size),
                ),
            ),
        )
    end

    # Defining the callback "default" means that the callback will be called
    # when a node is selected. ['block3','block4','block5'] are the names of
    # the blocks that the three returned objects will be assigned to. There is
    # no custom javascript function that will be called when the callback finishes.
    cb = DashWrapper.DashCallback(
        "default",
        ["block3", "block4", "block5"],
        cb_func,
    )

    # Create the Dash app with the blocks defined above, but no callbacks.
    # The extra file has a javascript function that defines an action to
    # take when a node is selected in the scenario tree.
    return JuDGE.create_dash(
        blocks,
        cb,
        path = @__DIR__,
        extra_files = [joinpath(@__DIR__, "plot_functions", "advanced.js")],
        data = treedata,
    )
end

function callback_dash()
    # Define the layout of the dashbaord
    blocks = JuDGE.dash_layout()
    push!(blocks, DashWrapper.Block((1, 0), (0.5, 6), "Callback output"))
    DashWrapper.scale_blocks!(blocks, 1200, 1000, 16)

    treedata = JuDGE.export_tree(mytree, data = solution)

    # A variable that will be modified by one callback and read by the
    # other.
    public_data = nothing

    # This callback function is attached to the "default" action.
    # This means that it is triggered when the user selects a node.
    # The data returned from this callback is sent to the javascript
    # function js_cb1.
    function cb_func1(arg)
        if arg[1] == 0
            public_data = nothing
            println("No node selected")
        else
            node = collect(mytree)[arg[1]]
            public_data = typeof(node) == Leaf
            println(node)
        end
        return public_data
    end

    # This callback function is attached to "custom", and is triggered when
    # the javascript function 'send_message(...,"custom")' is called.
    # The color returned from this callback is sent to the javascript
    # function js_cb2, and the arg is put into the block4 div.
    function cb_func2(arg)
        if public_data == nothing
            color = "#000000"
        elseif public_data == true
            color = "#DD2222"
        else
            color = "#2222DD"
        end
        return color, arg
    end

    # These functions are in the plot_functions/cb.js file
    cbs = [
        DashWrapper.DashCallback("default", "js_cb1", cb_func1),
        DashWrapper.DashCallback("custom", ["block3"], "js_cb2", cb_func2),
    ]

    # Create the Dash app with the blocks and callbacks specified.
    return JuDGE.create_dash(
        blocks,
        cbs,
        extra_files = [joinpath(@__DIR__, "plot_functions", "cb.js")],
        data = treedata,
    )
end

function run_dash(mode::Symbol)
    if mode == :simple
        app = simple_dash()
    elseif mode == :medium
        app = medium_dash()
    elseif mode == :advanced
        app = advanced_dash()
    elseif mode == :callback
        app = callback_dash()
    else
        error("Mode not defined")
    end

    # Loads the dashboard in the browser, and then runs the server.
    # The first time this is run, the browser may need to be manually
    # refreshed.
    DashWrapper.load_dash(app)
    return DashWrapper.run_server(app, "0.0.0.0", debug = false)
end

run_dash(:simple)
# run_dash(:medium)
# run_dash(:advanced)
# run_dash(:callback)
