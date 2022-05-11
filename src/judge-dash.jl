mutable struct Block
    pos::Tuple{Real,Real}
    size::Tuple{Real,Real}
    title::String
    div_name::String
end

function Block(pos::Tuple{Real,Real}, size::Tuple{Real,Real}, title::String)
    return Block(pos, size, title, "")
end

struct DashCallback
    block_indices::Vector{Real}
    func::Function
    js_fn::String
    name::String
end

function bottom(block::Block)
    return block.pos[2] + block.size[2]
end

function right(block::Block)
    return block.pos[1] + block.size[1]
end

function scale_blocks!(blocks::Vector{Block}, width::Int, height::Int, gap::Int)
    push!(blocks, Block((0, 0), (1, 5), "JuDGE Scenario Tree", "tree"))
    push!(blocks, Block((0, 5), (1, 1), "Settings", "settings"))

    max_x = maximum(right.(blocks))
    max_y = maximum(bottom.(blocks))

    sizes = ((width + gap) / max_x - gap, (height + gap) / max_y - gap)

    for block in blocks
        block.pos =
            (block.pos[1] * (sizes[1] + gap), block.pos[2] * (sizes[2] + gap))
        block.size = (
            block.size[1] * (sizes[1] + gap) - gap,
            block.size[2] * (sizes[2] + gap) - gap,
        )
    end

    return
end

function export_tree(
    some_tree::AbstractTree;
    scale_edges = nothing,
    scale_nodes::Float64 = 0.0,
    max_size::Float64 = 50.0,
    truncate::Int = -1,
    rel_angle::Bool = false,
    style::Symbol = :standard,
    box_size::Int = 800,
    skip_root::Bool = false,
    data::Dict = Dict(),
)
    maxdata = Dict{Symbol,Any}()
    mindata = Dict{Symbol,Any}()
    scale_range = Dict{Symbol,Dict{Symbol,Any}}()

    function data_scale(node)
        #        if data_from_original
        #            node = getID(node)
        #        end
        if skip_root && node.parent == nothing
            return
        end
        for sym in keys(data[node])
            if typeof(data[node][sym]) == Float64
                if !haskey(scale_range, sym)
                    scale_range[sym] = Dict{Symbol,Any}()
                    scale_range[sym][:min] = data[node][sym]
                    scale_range[sym][:max] = data[node][sym]
                else
                    if data[node][sym] < scale_range[sym][:min]
                        scale_range[sym][:min] = data[node][sym]
                    elseif data[node][sym] > scale_range[sym][:max]
                        scale_range[sym][:max] = data[node][sym]
                    end
                end
            elseif typeof(data[node][sym]) <: Dict
                for key in keys(data[node][sym])
                    if !haskey(scale_range, sym)
                        scale_range[sym] = Dict{Symbol,Any}()
                        scale_range[sym][:max] = data[node][sym][key]
                        scale_range[sym][:min] = data[node][sym][key]
                    else
                        if data[node][sym][key] < scale_range[sym][:min]
                            scale_range[sym][:min] = data[node][sym][key]
                        elseif data[node][sym][key] > scale_range[sym][:max]
                            scale_range[sym][:max] = data[node][sym][key]
                        end
                    end
                end
            elseif typeof(data[node][sym]) <: Array
                for val in data[node][sym]
                    if !haskey(scale_range, sym)
                        scale_range[sym] = Dict{Symbol,Any}()
                        scale_range[sym][:max] = val
                        scale_range[sym][:min] = val
                    else
                        if val < scale_range[sym][:min]
                            scale_range[sym][:min] = val
                        elseif val > scale_range[sym][:max]
                            scale_range[sym][:max] = val
                        end
                    end
                end
            end
        end
    end

    function arc_json(node, parent)
        return Dict("from" => get_id[parent] + 1, "to" => get_id[node] + 1)
    end

    get_id = Dict{AbstractTree,Int}()

    angles = Dict{AbstractTree,Float64}()
    position = Dict{AbstractTree,Tuple{Float64,Float64}}()
    position2 = Dict{AbstractTree,Tuple{Float64,Float64}}()

    function setpositions(
        node::AbstractTree,
        rel_angle::Bool,
        l::Float64,
        scale::Float64,
        odd::Bool,
    )
        if typeof(node) == Leaf
            return
        end
        a = -2 * pi / (length(node.children))
        current = 0.0
        if rel_angle
            current = angles[node] + pi + a / 2#0.0
        elseif (length(node.children) % 2) == 1
            current +=
                pi / 2 +
                (length(node.children) - 1) * pi / length(node.children)
        elseif odd
            current += 3 * pi / 2 - pi / length(node.children)
        else
            current += 3 * pi / 2 - 2 * pi / length(node.children)
        end

        for child in node.children
            angles[child] = current
            position[child] = (
                position[node][1] + l * cos(current),
                position[node][2] - l * sin(current),
            )
            if length(node.children) == 2
                if odd
                    setpositions(child, rel_angle, l, scale, !odd)
                else
                    setpositions(child, rel_angle, l * scale^2, scale, !odd)
                end
            else
                setpositions(child, rel_angle, l * scale, scale, !odd)
            end
            current += a
        end
    end

    function setpositions2(node::AbstractTree, leaf_sep::Float64)
        function locate(node::AbstractTree, vert::Float64, horz::Float64)
            if typeof(node) == Leaf ||
               (truncate != -1 && JuDGE.depth(node) >= truncate)
                vert += leaf_sep
                position2[node] = (horz, vert)
                return (vert, vert)
            else
                verts = Float64[]
                for child in node.children
                    pos, vert = locate(child, vert, horz + parch_sep)
                    push!(verts, pos)
                end
                position2[node] = (horz, sum(verts) / length(verts))
                return (position2[node][2], vert)
            end
        end
        num_leaf = length(JuDGE.get_leafnodes(node, truncate = truncate))
        max_depth = JuDGE.depth(collect(node, truncate = truncate)[end])
        parch_sep = 0.8 * leaf_sep * num_leaf / max_depth
        return locate(node, 0.0, 0.0)
    end

    scale_factors = [
        [1.0],
        [1.0, 0.87, 0.83, 0.78, 0.74, 0.71, 0.695, 0.685],
        [1.0, 0.65, 0.52, 0.48],
        [1.0, 0.45, 0.42, 0.42, 0.41, 0.4],
        [1.0, 0.44, 0.37, 0.36],
        [1.0, 0.42, 0.34, 0.33],
        [1.0, 0.35, 0.3],
        [1.0, 0.3, 0.26],
        [1.0, 0.27, 0.23],
        [1.0, 0.24, 0.22],
    ]

    if scale_edges == nothing
        if typeof(some_tree) == Leaf
            scale_edges = 1.0
        else
            dg = length(some_tree.children)
            if dg <= 10
                dp = min(
                    JuDGE.depth(collect(some_tree, truncate = truncate)[end]),
                    length(scale_factors[dg]),
                )
                scale_edges = scale_factors[dg][dp]
            else
                dp = JuDGE.depth(collect(some_tree, truncate = truncate)[end])
                scale_edges = 0.22 * 0.98^(dg^dp - 11)
            end
        end
    end

    angles[some_tree] = rel_angle ? 0.0 : -pi / 2
    position[some_tree] = (0.0, 0.0)

    setpositions(some_tree, rel_angle, 700.0, scale_edges, true)
    setpositions2(some_tree, 70.0)
    index = 0
    if skip_root
        index -= 1
    end

    arcs = Dict[]
    nodes = Dict[]

    for node in collect(some_tree, truncate = truncate)
        get_id[node] = index
        parent = node.parent
        if parent != nothing && (!skip_root || parent.parent != nothing)
            push!(arcs, arc_json(node, parent))
        end
        index += 1
    end

    for node in collect(some_tree, truncate = truncate)
        if !skip_root || node != some_tree
            data_scale(node)
            temp = Dict{Symbol,Any}()
            temp[:id] = get_id[node] + 1
            temp[:label] = node.name
            temp[:level] = JuDGE.depth(node)
            temp[:leaf] =
                typeof(node) == Leaf ||
                (truncate != -1 && JuDGE.depth(node) >= truncate)
            temp[:posX] = position[node][1]
            temp[:posY] = position[node][2]
            temp[:posX2] = position2[node][1]
            temp[:posY2] = position2[node][2]
            if haskey(data, node)
                temp[:data] = data[node]
            end
            push!(nodes, temp)
        end
    end

    if scale_nodes == 0.0
        scale_nodes = min(
            1.0,
            exp(
                log(
                    (
                        400 *
                        scale_edges^JuDGE.depth(
                            collect(some_tree, truncate = truncate)[end],
                        )
                    ) / max_size,
                ) /
                JuDGE.depth(collect(some_tree, truncate = truncate)[end]),
            ),
        )
    end
    min_size = 25

    return Dict(
        "nodes" => nodes,
        "arcs" => arcs,
        "node_scale" => scale_nodes,
        "min_size" => min_size,
        "max_size" => max_size,
        "scale" => scale_range,
    )
end

function make_window(
    ID::String,
    pos::Tuple{Real,Real},
    size::Tuple{Real,Real};
    title::String = "",
    content = nothing,
)
    inside = []
    if title != ""
        push!(inside, html_div(className = "topbar grey") do
            return title
        end)
    end

    if content == nothing
        push!(inside, html_div(id = ID * "_content", className = "content"))
    else
        push!(inside, html_div(id = ID * "_content", className = "content") do
            return content
        end)
    end

    return html_div(
        id = ID * "_frame",
        className = "frame",
        style = Dict(
            "left" => string(pos[1]) * "px",
            "top" => string(pos[2]) * "px",
            "width" => string(size[1]) * "px",
            "height" => string(size[2]) * "px",
        ),
    ) do
        return inside
    end
end

function load_dash()
    return JuDGE.run_file(
        joinpath(dirname(@__DIR__), "dash", "judge-dash.html"),
    )
end

function create_assets(
    path::String,
    extra_files::Vector{String},
    callbacks::Vector{DashCallback},
)
    srcdir = joinpath(dirname(@__DIR__), "dash")
    dstdir = joinpath(path, "assets")

    try
        rm(dstdir; force = true, recursive = true)
    catch
    end

    if !isdir(dstdir)
        mkdir(dstdir)
    end

    cp(
        joinpath(srcdir, "colors.js"),
        joinpath(dstdir, "colors.js"),
        force = true,
    )
    cp(joinpath(srcdir, "tree.js"), joinpath(dstdir, "tree.js"), force = true)
    cp(
        joinpath(srcdir, "judge-dash.js"),
        joinpath(dstdir, "judge-dash.js"),
        force = true,
    )
    cp(
        joinpath(srcdir, "jdash.css"),
        joinpath(dstdir, "jdash.css"),
        force = true,
    )
    cp(
        joinpath(srcdir, "jsonview", "jsonview.bundle.js"),
        joinpath(dstdir, "jsonview.bundle.js"),
        force = true,
    )
    cp(
        joinpath(srcdir, "jsonview", "jsonview.bundle.css"),
        joinpath(dstdir, "jsonview.bundle.css"),
        force = true,
    )
    cp(
        joinpath(srcdir, "svg-pan-zoom", "svg-pan-zoom.min.js"),
        joinpath(dstdir, "svg-pan-zoom.min.js"),
        force = true,
    )

    for file in extra_files
        f = file[findlast('\\', file)+1:end]
        println(file)
        println(f)
        println(joinpath(dstdir, f))
        cp(file, joinpath(dstdir, f), force = true)
    end

    if length(callbacks) != 0
        file = open(joinpath(dstdir, "callbacks.js"), "w")

        println(file, "function callback_controller(trigger,data) {")
        println(file, "\tswitch (trigger) {")
        for cb in callbacks
            println(file, "\t\tcase \"" * cb.name * "\":")
            println(file, "\t\t\t" * cb.js_fn * "(data);")
            #println(file, "\t\t\tconsole.log(data);")
            println(file, "\t\t\tbreak;")
        end
        println(file, "\t\tdefault:")
        println(file, "\t\t\tconsole.log(\"Unknown trigger\");")
        println(file, "\t}")
        println(file, "}")
        close(file)
    end
end

function create_dash(
    tree::AbstractTree,
    solution::Dict,
    blocks::Union{Nothing,Vector{Block}},
    callbacks::Union{Nothing,DashCallback,Vector{DashCallback}};
    path::String = pwd(),
    extra_files::Vector{String} = String[],
)
    if typeof(callbacks) == DashCallback
        callbacks = [callbacks]
    elseif callbacks == nothing
        callbacks = DashCallback[]
    end

    if blocks == nothing
        blocks = Block[]
        scale_blocks!(blocks, 800, 1000, 16)
    end

    windows = []

    i = 0

    if length(blocks) > 0
        for block in blocks
            i += 1
            if block.div_name == ""
                push!(
                    windows,
                    make_window(
                        "block" * string(i),
                        block.pos,
                        block.size,
                        title = block.title,
                    ),
                )
            else
                push!(
                    windows,
                    make_window(
                        block.div_name,
                        block.pos,
                        block.size,
                        title = block.title,
                    ),
                )
            end
        end
    end
    treedata = export_tree(tree, data = solution)

    treedata["callbacks"] = String[]

    cb_default = 0
    js_cb = false
    for cb in callbacks
        if cb.name != "default"
            push!(
                windows,
                html_button(
                    id = "cb-" * cb.name * "-button",
                    children = cb.name,
                    n_clicks = 0,
                    style = Dict("display" => "none"),
                ),
            )
            push!(
                windows,
                dcc_store(
                    id = "msg-" * cb.name,
                    storage_type = "session",
                    data = "",
                ),
            )
            if cb.js_fn != ""
                js_cb = true
                push!(
                    windows,
                    dcc_store(
                        id = "data-" * cb.name,
                        storage_type = "session",
                        data = "",
                    ),
                )
                push!(treedata["callbacks"], cb.name)
                push!(
                    windows,
                    html_div(
                        id = "cb-" * cb.name * "_content",
                        style = Dict("display" => "none"),
                    ),
                )
            end
        else
            cb_default += 1
            if cb.js_fn != ""
                js_cb = true
                push!(treedata["callbacks"], "default")
                push!(
                    windows,
                    dcc_store(
                        id = "data-" * cb.name,
                        storage_type = "session",
                        data = "",
                    ),
                )
            end
        end
    end

    push!(
        windows,
        dcc_store(id = "tree_data", storage_type = "session", data = treedata),
    )
    push!(
        windows,
        make_window(
            "cb-default",
            (
                blocks[end-1].pos[1] + 16,
                blocks[end-1].pos[2] + blocks[end-1].size[2] - 28,
            ),
            (200, 28),
        ),
    )
    push!(
        windows,
        dcc_store(id = "msg-default", storage_type = "session", data = "0;;"),
    )
    push!(
        windows,
        html_button(
            id = "cb-default-button",
            children = "default",
            n_clicks = 0,
            style = Dict("display" => "none"),
        ),
    )

    if cb_default > 1
        error("Can only specify one default callback.")
    end

    offset = maximum(right.(blocks)) * 0.15

    create_assets(path, extra_files, js_cb ? callbacks : DashCallback[])

    backup_dir = pwd()

    cd(path)
    app = dash(
        external_scripts = [
            "https://cdnjs.cloudflare.com/ajax/libs/jquery/2.1.3/jquery.min.js",
        ],
        external_stylesheets = [],
    )
    cd(backup_dir)

    app.title = "JuDGE Solution Dashboard"

    app.layout =
        app.layout = html_div(
            id = "demo",
            style = Dict(
                "width" => string(maximum(right.(blocks))) * "px",
                "padding-top" =>
                    "clamp(" *
                    string(offset) *
                    "px,15%," *
                    string(2 * offset) *
                    "px)",
                "margin-top" => "-" * string(offset) * "px",
                "margin-left" => "auto",
                "margin-right" => "auto",
            ),
        ) do
            html_div(className = "full") do
                return windows
            end
        end

    if length(callbacks) == 0
        callback!(
            app,
            Output("cb-default_content", "children"),
            Input("cb-default-button", "n_clicks"),
            State("msg-default", "data"),
            prevent_initial_call = true,
        ) do input_1, input_2
            index, arg = process_cb_input(input_2)
            nodes = collect(tree)
            if index == 0
                return "No node selected"
            else
                return "Selected: " * nodes[index].name
            end
        end
    else
        for cb in callbacks
            if cb.name == "default"
                outputs = [Output("cb-default_content", "children")]

                if cb.js_fn != ""
                    push!(outputs, Output("data-default", "data"))
                end

                for index in cb.block_indices
                    push!(
                        outputs,
                        Output(
                            "block" * string(index) * "_content",
                            "children",
                        ),
                    )
                end

                if length(outputs) == 1
                    error("This callback doesn't do anything. Delete it.")
                end

                callback!(
                    app,
                    outputs,
                    Input("cb-default-button", "n_clicks"),
                    State("msg-default", "data"),
                    prevent_initial_call = true,
                ) do input_1, input_2
                    index, arg = process_cb_input(input_2)
                    nodes = collect(tree)

                    if index == 0
                        out = ("No node selected",)
                        out2 = cb.func(nothing, arg)
                    else
                        out = ("Selected: " * nodes[index].name,)
                        out2 = cb.func(nodes[index], arg)
                    end
                    if !(typeof(out2) <: Tuple)
                        return join_tuples(out, (out2,))
                    else
                        return join_tuples(out, out2)
                    end
                end
            else
                outputs = Output{String}[]
                if cb.js_fn != ""
                    push!(
                        outputs,
                        Output("cb-" * cb.name * "_content", "children"),
                    )
                    push!(outputs, Output("data-" * cb.name, "data"))
                end
                for index in cb.block_indices
                    push!(
                        outputs,
                        Output(
                            "block" * string(index) * "_content",
                            "children",
                        ),
                    )
                end

                if length(outputs) == 0
                    error("This callback doesn't do anything. Delete it.")
                end

                callback!(
                    app,
                    outputs,
                    Input("cb-" * cb.name * "-button", "n_clicks"),
                    State("msg-" * cb.name, "data"),
                    #prevent_initial_call=true,
                ) do input_1, input_2
                    print(cb.name * ": ")
                    println(input_2)

                    out = cb.func(input_2)

                    if !(typeof(out) <: Tuple)
                        out = (out,)
                    end

                    if cb.js_fn != ""
                        return join_tuples((input_1,), out)
                    else
                        return out
                    end
                end
            end
        end
    end
    #app.root_path=path
    return app
end

function process_cb_input(input)
    println("Message received: " * string(input))
    temp = split(input, ";")
    index = parse(Int, temp[1])
    deleteat!(temp, 1)
    return index, string.(temp)
end
