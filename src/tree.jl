using JSON

# trait for iterablility
abstract type IterableTrait end
struct Iterable{T} <: IterableTrait end
struct NotIterable{T} <: IterableTrait end
IterableTrait(::Type) = NotIterable{Any}()
IterableTrait(::Type{A}) where {A<:AbstractArray{T}} where {T} = Iterable{T}()
function IterableTrait(::Type{A}) where {A<:Dict{U,T}} where {T} where {U}
    return Iterable{T}()
end

# definition of a tree
abstract type AbstractTree end

mutable struct Leaf <: AbstractTree
    ID::Union{Nothing,AbstractTree}
    name::String
    parent::Union{Nothing,AbstractTree}
    ext::Dict{Symbol,Any}
    function Leaf()
        return new(nothing, "", nothing, Dict{Symbol,Any}())
    end
    function Leaf(parent::AbstractTree)
        return new(nothing, "", parent, Dict{Symbol,Any}())
    end
    function Leaf(name::String)
        return new(nothing, name, nothing, Dict{Symbol,Any}())
    end
    function Leaf(name::String, parent::AbstractTree)
        return new(nothing, name, parent, Dict{Symbol,Any}())
    end
end

mutable struct Tree <: AbstractTree
    children::Union{Vector{AbstractTree},Vector{Leaf},Vector{Tree}}
    ID::Union{Nothing,AbstractTree}
    name::String
    parent::Union{Nothing,AbstractTree}
    ext::Dict{Symbol,Any}
    Tree(children) = Tree(IterableTrait(typeof(children)), children)
    function Tree(::Iterable{<:AbstractTree}, children)
        return new(children, nothing, "", nothing, Dict{Symbol,Any}())
    end
    function Tree()
        return new(
            Vector{AbstractTree}(),
            nothing,
            "",
            nothing,
            Dict{Symbol,Any}(),
        )
    end
    function Tree(parent::AbstractTree)
        return new(
            Vector{AbstractTree}(),
            nothing,
            "",
            parent,
            Dict{Symbol,Any}(),
        )
    end
    function Tree(name::String)
        return new(
            Vector{AbstractTree}(),
            nothing,
            name,
            nothing,
            Dict{Symbol,Any}(),
        )
    end
    function Tree(name::String, parent::AbstractTree)
        return new(
            Vector{AbstractTree}(),
            nothing,
            name,
            parent,
            Dict{Symbol,Any}(),
        )
    end
end

function Base.copy(leaf::Leaf)
    return Leaf()
end

function Base.copy(tree::Tree)
    return Tree(map(copy, tree.children))
end

function Base.map(f, leaf::Leaf)
    return f(leaf)
end

function Base.map(f, tree::Tree)
    return Tree(map(x -> map(f, x), tree.children))
end

function Base.map(f, dict::Dict)
    return Dict(key => f(dict[key]) for key in keys(dict))
end

function Base.show(io::IO, tree::AbstractTree)
    if typeof(tree) == Tree
        if tree.name != ""
            print(
                io,
                "Subtree rooted at node " *
                tree.name *
                " containing " *
                string(count(tree)) *
                " nodes",
            )
        else
            print(io, "Subtree containing " * string(count(tree)) * " nodes")
        end
    else
        if tree.name != ""
            print(io, "Leaf node " * tree.name)
        else
            print(io, "Leaf node")
        end
    end
end

function getID(node::AbstractTree)
    if node.ID == nothing
        return node
    else
        return node.ID
    end
end

function count(tree::AbstractTree)
    i = 1
    if typeof(tree) == Tree
        for child in tree.children
            i += count(child)
        end
    end
    return i
end

function print_tree(some_tree::AbstractTree)
    function helper(tree::Tree, depth)
        println("  "^depth * "--" * tree.name)
        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        return println("  "^depth * "--" * leaf.name)
    end
    helper(some_tree, 0)
    return nothing
end

"""
	print_tree(some_tree::AbstractTree, data::Dict{AbstractTree,Any})

Given `some_tree`, this function prints a simple representation of the tree to the REPL.

### Required Arguments
`some_tree` is the tree we wish to visualise

### Optional Arguments
`data` is a dictionary indexed by each node in `some_tree`
"""
function print_tree(
    some_tree::AbstractTree,
    data::Dict{AbstractTree,T} where {T<:Any},
)
    function helper(tree::Tree, depth)
        if tree in keys(data)
            println(
                "  "^depth * "--" * tree.name * " (" * string(data[tree]) * ")",
            )
        elseif getID(tree) in keys(data)
            println(
                "  "^depth *
                "--" *
                tree.name *
                " (" *
                string(data[getID(tree)]) *
                ")",
            )
        else
            println("  "^depth * "--" * tree.name)
        end

        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        if leaf in keys(data)
            println(
                "  "^depth * "--" * leaf.name * " (" * string(data[leaf]) * ")",
            )
        elseif getID(leaf) in keys(data)
            println(
                "  "^depth *
                "--" *
                leaf.name *
                " (" *
                string(data[getID(leaf)]) *
                ")",
            )
        else
            println("  "^depth * "--" * leaf.name)
        end
    end
    helper(some_tree, 0)

    return nothing
end

"""
	visualize_tree(some_tree::AbstractTree,
        data::Dict{AbstractTree,Dict{Symbol,Any}};
        scale_edges=nothing,
        scale_all=1.0)

Given `some_tree`, this function generates a html/js visualization of the tree.

### Required Arguments
`some_tree` is the tree we wish to visualise.

`data` is a dictionary of the data we wish to display, each element is another
dictionary indexed by the nodes of the tree.

### Optional Arguments
`scale_edges` this is the scale factor for edges as the network gets deeper.

`scale_nodes` this is the scale factor for nodes as the network gets deeper.

`max_size` this is the size of the root node.
"""
function visualize_tree(
    some_tree::AbstractTree,
    data::Dict{AbstractTree,Dict{Symbol,Any}};
    scale_edges = nothing,
    scale_nodes = 0.0,
    max_size = 50.0,
    custom::Union{Nothing,Dict{Symbol,Tuple{String,String,String}}} = nothing,
)
    maxdata = Dict{Symbol,Any}()
    mindata = Dict{Symbol,Any}()

    function node_json(node)
        first = true
        for sym in keys(data[node])
            if sym == :graph_data || sym == :custom_data
                continue
            end
            if typeof(data[node][sym]) == Float64
                if sym ∉ keys(maxdata)
                    maxdata[sym] = data[node][sym]
                    mindata[sym] = data[node][sym]
                else
                    if data[node][sym] < mindata[sym]
                        mindata[sym] = data[node][sym]
                    elseif data[node][sym] > maxdata[sym]
                        maxdata[sym] = data[node][sym]
                    end
                end
            elseif typeof(data[node][sym]) <: Dict
                for key in keys(data[node][sym])
                    if sym ∉ keys(maxdata)
                        maxdata[sym] = data[node][sym][key]
                        mindata[sym] = data[node][sym][key]
                    else
                        if data[node][sym][key] < mindata[sym]
                            mindata[sym] = data[node][sym][key]
                        elseif data[node][sym][key] > maxdata[sym]
                            maxdata[sym] = data[node][sym][key]
                        end
                    end
                end
            end
        end
    end

    function arc_json(node, parent)
        temp = "{"
        temp *= "from:"
        temp *= string(get_id[parent])
        temp *= ",to:"
        temp *= string(get_id[node])
        return temp *= ",},"
    end

    index = 1
    nodes = "var nodes = ["
    arcs = "var edges = ["
    get_id = Dict{AbstractTree,Int}()

    angles = Dict{AbstractTree,Float64}()
    position = Dict{AbstractTree,Tuple{Float64,Float64}}()
    position2 = Dict{AbstractTree,Tuple{Float64,Float64}}()

    function setpositions(
        node::AbstractTree,
        alpha::Float64,
        l::Float64,
        scale::Float64,
        odd::Bool,
    )
        if typeof(node) == Leaf
            return
        end
        a = 0.0
        current = 0.0
        if (length(node.children) % 2) == 1
            a = -2 * pi / (length(node.children))
            current =
                pi / 2 +
                (length(node.children) - 1) * pi / length(node.children)
        elseif odd
            a = -2 * pi / (length(node.children))
            current = 3 * pi / 2 - pi / length(node.children)
        else
            a = -2 * pi / (length(node.children))
            current = 3 * pi / 2 - 2 * pi / length(node.children)
        end

        for child in node.children
            angles[child] = current
            position[child] = (
                position[node][1] + l * cos(current),
                position[node][2] - l * sin(current),
            )
            if length(node.children) == 2
                if odd
                    setpositions(child, alpha, l, scale, !odd)
                else
                    setpositions(child, alpha, l * scale^2, scale, !odd)
                end
            else
                setpositions(child, alpha, l * scale, scale, !odd)
            end
            current += a
        end
    end

    function setpositions2(node::AbstractTree, leaf_sep::Float64)
        function locate(node::AbstractTree, vert::Float64, horz::Float64)
            if typeof(node) == Leaf
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
        num_leaf = length(get_leafnodes(node))
        max_depth = depth(collect(node)[end])
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
                    depth(collect(some_tree)[end]),
                    length(scale_factors[dg]),
                )
                scale_edges = scale_factors[dg][dp]
            else
                scale_edges = 0.22 * 0.91^(dg - 10)
            end
        end
    end

    angles[some_tree] = -pi / 2
    position[some_tree] = (0.0, 0.0)

    setpositions(some_tree, 0.5, 700.0, scale_edges, true)
    setpositions2(some_tree, 80.0)

    for node in collect(some_tree)
        get_id[node] = index
        node_json(node)
        parent = node.parent
        if parent != nothing
            arcs *= arc_json(node, parent)
        end
        index += 1
    end

    nodejson = Dict{AbstractTree,Dict{Symbol,Any}}()
    for node in collect(some_tree)
        nodejson[node] = Dict{Symbol,Any}()
        nodejson[node][:id] = get_id[node]
        nodejson[node][:label] = node.name
        nodejson[node][:level] = depth(node)
        nodejson[node][:posX] = position[node][1]
        nodejson[node][:posY] = position[node][2]
        nodejson[node][:posX2] = position2[node][1]
        nodejson[node][:posY2] = position2[node][2]
        nodejson[node][:data] = data[node]
    end

    for node in collect(some_tree)
        nodes *= JSON.json(nodejson[node])
        nodes *= ","
    end

    nodes *= "];"
    arcs *= "];"

    scale = "var scale = {"

    for sym in keys(maxdata)
        if length(scale) != 13
            scale *= ","
        end
        scale *= string(sym)
        scale *= ":"
        scale *= "{min:"
        scale *= string(mindata[sym])
        scale *= ",max:"
        if maxdata[sym] == mindata[sym]
            maxdata[sym] += 1
        end
        scale *= string(maxdata[sym])
        scale *= "}"
    end
    scale *= "};"
    if scale_nodes == 0.0
        scale_nodes = min(
            1.0,
            exp(
                log(
                    (400 * scale_edges^depth(collect(some_tree)[end])) /
                    max_size,
                ) / depth(collect(some_tree)[end]),
            ),
        )
    end
    min_size = 25
    adata =
        nodes *
        "\n" *
        arcs *
        "\n" *
        scale *
        "\n" *
        "var node_scale=" *
        string(scale_nodes) *
        ";" *
        "var min_size=" *
        string(min_size) *
        ";" *
        "var max_size=" *
        string(max_size) *
        ";"

    s = read(joinpath(dirname(@__DIR__), "visualise", "visualize.html"), String)
    c = read(joinpath(dirname(@__DIR__), "visualise", "colors.js"), String)
    filename =
        joinpath("vis" * string(Int(round((time() * 100) % 100000))) * ".html")
    file = open(filename, "w")
    s = replace(replace(s, "#DATA#" => adata), "#COLORS#" => c)
    if custom != nothing
        temp = ""
        for cf in values(custom)
            temp *= read(cf[3], String) * "\n"
        end
        s = replace(s, "#FUNCTIONS#" => temp)

        temp = "if ('custom_data' in nodes[selected-1].data) {\n"
        for (key, cf) in custom
            temp *=
                "if ('" *
                string(key) *
                "' in nodes[selected-1].data.custom_data) {\n"
            temp *=
                cf[2] *
                "(nodes[selected-1].data.custom_data." *
                string(key) *
                ");\n"
            temp *= "}\n"
            temp *= "else {\n"
            temp *= cf[2] * "(null);\n"
            temp *= "}\n"
        end
        temp *= "}\n"
        temp *= "else {\n"
        for (key, cf) in custom
            temp *= cf[2] * "(null);\n"
        end
        temp *= "}\n"
        s = replace(s, "#FUNCTION_CALLS#" => temp)

        temp = ""
        for cf in values(custom)
            temp *= cf[1]
        end
        s = replace(s, "#DIVS#" => temp)
    else
        s = replace(s, "#FUNCTIONS#" => "")
        s = replace(s, "#FUNCTION_CALLS#" => "")
        s = replace(s, "#DIVS#" => "")
    end
    println(file, s)
    close(file)

    if Sys.iswindows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    elseif Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
    else
        error("Unable to show plot. Try opening the file $(filename) manually.")
    end

    return
end

"""
	collect(tree::Tree;order=:depth)

Given `tree`, this function returns an array of corresponding nodes. By default this will be in a depth-first order.

### Required Arguments
`tree` is the tree from which we wish to collect the nodes

### Optional Arguments
`order` can be set to `:depth` or `:breadth` to specify the order that the nodes are listed in the array.

### Examples
    nodes = collect(tree) #gets an array of nodes from tree in depth-first order
    nodes = collect(tree,order=:breadth) #gets an array of nodes from tree in breadth-first order
"""
function Base.collect(tree::Tree; order = :depth)
    index = 1
    collection = Vector{AbstractTree}()
    push!(collection, tree)
    while index <= length(collection)
        if typeof(collection[index]) == Tree
            if order == :depth
                for i in eachindex(collection[index].children)
                    insert!(
                        collection,
                        index + i,
                        collection[index].children[i],
                    )
                end
            elseif order == :breadth
                append!(collection, collection[index].children)
            else
                error("\'order\' should be set to :depth or :breadth")
            end
        end
        index += 1
    end
    return collection
end

function Base.collect(leaf::Leaf; order = :depth)
    return [leaf]
end

function label_nodes(tree::AbstractTree)
    for t in collect(tree)
        t.name = history2(t)
    end
    return tree
end

"""
	get_leafnodes(tree::AbstractTree)

Given `tree`, this function returns an array of corresponding `Leaf` nodes.

### Required Arguments
`tree` is the tree from which we wish to collect leaf nodes

### Example
    leafnodes = JuDGE.get_leafnodes(tree)
"""
function get_leafnodes(tree::AbstractTree)
    function helper(leaf::Leaf, collection)
        return push!(collection, leaf)
    end
    function helper(someTree::Tree, collection)
        for child in someTree.children
            helper(child, collection)
        end
    end
    result = Vector{Leaf}()
    helper(tree, result)
    return result
end

"""
	ConditionallyUniformProbabilities(tree::AbstractTree)

Given a tree, this function returns a dictionary which maps nodes of the tree to probabilities, given that there are conditionally uniform probabilities over the children of any node.

### Required Arguments
`tree` is the tree for which the probability distribution will be generated

### Example
    probs = ConditionallyUniformProbabilities(tree)
"""
function ConditionallyUniformProbabilities(tree::T where {T<:AbstractTree})
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return 1.0
        else
            parent = subtree.parent
            return result(parent) / length(parent.children)
        end
    end
    prob = Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node] = result(node)
    end
    return prob
end

"""
	UniformLeafProbabilities(tree::AbstractTree)

Given a tree, this function returns a dictionary which maps nodes of the tree to probabilities, given that there is a uniform distribution over the leaf nodes

### Required Arguments
`tree` is the tree for which the probability distribution will be generated

### Example
    probs = UniformLeafProbabilities(tree)
"""
function UniformLeafProbabilities(tree::T where {T<:AbstractTree})
    function result(subtree::T where {T<:AbstractTree})
        if typeof(subtree) == Leaf
            return p
        else
            pr = 0
            for child in subtree.children
                pr += result(child)
            end
            return pr
        end
    end

    leafnodes = get_leafnodes(tree)
    p = 1 / length(leafnodes)
    prob = Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node] = result(node)
    end
    return prob
end

"""
	convert_probabilities(tree::AbstractTree, probabilities::Dict{NodeID,Float64})

Given a dictionary of conditional probabilities for each node in `tree`, this function
returns a dictionary that maps the `NodeID` of each node of `tree` to the corresponding
unconditional probability.

### Required Arguments
`tree` is the tree that the probabilities pertain to

`probabilities` is a dictionary of condition probabilities for each node in `tree`

### Example
    probs = JuDGE.convert_probabilities(tree,probabilities)
"""
function convert_probabilities(
    tree::T where {T<:AbstractTree},
    probabilities::Dict{AbstractTree,Float64},
)
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return probabilities[getID(subtree)]
        else
            parent = subtree.parent
            return result(parent) * probabilities[getID(subtree)]
        end
    end
    prob = Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node] = result(node)
    end
    return prob
end

"""
	depth(tree::AbstractTree)

Given `tree`, this function returns the depth.
The root node has a depth of 0.

### Required Arguments
`tree` is the node we wish to find the depth for.

### Example
    depth = JuDGE.depth(tree) #returns the depth of a node in a tree
"""
function depth(tree::AbstractTree)
    return length(history(tree)) - 1
end

"""
	history(tree::AbstractTree)

Given `tree`, this function returns the history back up to the root node of `tree`.

### Required Arguments
`tree` is the node that we wish to find the history for.

### Example
    history = JuDGE.history(tree) #get a vector of nodes that precede tree
"""
function history(tree::T where {T<:AbstractTree})
    function helper(state, subtree)
        if subtree.parent == nothing
            push!(state, subtree)
            return state
        else
            push!(state, subtree)
            helper(state, subtree.parent)
        end
        return state
    end
    return helper(Vector{AbstractTree}(), tree)
end

function history2(tree::T where {T<:AbstractTree})
    function helper(state, subtree)
        if subtree.parent == nothing
            state = "1" * state
            return state
        else
            parent = subtree.parent
            for i in 1:length(parent.children)
                if subtree == parent.children[i]
                    state = string(i) * state
                    break
                end
            end
            state = helper(state, subtree.parent)
        end
        return state
    end
    return helper("", tree)
end

"""
	narytree(depth::Int, degree::Int)

Given the `depth` and `degree`, this function returns an N-ary tree. Note that a depth of 0 return a single `Leaf` node (which is also the root node of the tree).

### Required Arguments
`depth` is the maximum number of arcs from the root node any `Leaf` node

`degree` is the number of children of all nodes, other than the `Leaf` nodes

### Example
    tree = narytree(2,2)
"""
function narytree(depth::Int, degree::Int)
    if depth == 0
        return Leaf("1")
    end

    tree = Tree("1")
    layer = [tree]
    for dp in 1:depth
        nextlayer = AbstractTree[]
        for t in layer
            for dg in 1:degree
                if dp != depth
                    nt = Tree(t.name * string(dg), t)
                    push!(t.children, nt)
                    push!(nextlayer, nt)
                else
                    push!(t.children, Leaf(t.name * string(dg), t))
                end
            end
        end
        layer = nextlayer
    end
    #label_nodes(tree)
    return tree
end

"""
	get_node(tree::AbstractTree, indices::Vector{Int})

Given a `tree`, and an array of `indices`, this function returns the corresponding node in the tree.

### Required Arguments
`tree` is the tree from which we are finding the node

`indicies` is an array of integer indices identifying a node within `tree`.

### Examples
    node = get_node(tree,[1]) #get the root node
    node = get_node(tree,[1,1]) #get the first child of the root node
    node = get_node(tree,[1,2]) #get the second child of the root node
"""
function get_node(tree::AbstractTree, indices::Vector{Int})
    node = tree
    for i in 2:length(indices)
        node = node.children[indices[i]]
    end
    return node
end

"""
	get_groups(tree::AbstractTree; combine=0)

Given a `tree`, this function will split it up into an array of subtrees. These can
be provided as `blocks` for the `JuDGE.solve()` function to perform partial pricing.

### Required Arguments
`tree` is the tree from which we are finding the node

### Optional Arguments
`combine` this parameter determines the size of the subtrees (higher creates larger subtrees).
If set to 0, the subtrees will be the sets of paths to the leaf nodes.

### Examples
    blocks = get_groups(tree)
"""
function get_groups(tree::AbstractTree; combine = 0)
    leafnodes = get_leafnodes(tree)

    groups = Vector{AbstractTree}[]
    nodes = AbstractTree[]

    for leaf in leafnodes
        node = history(leaf)[combine+1]
        if node ∉ nodes
            push!(nodes, node)
        end
    end

    for node in nodes
        scenario = history(node)
        children = collect(node)
        popfirst!(children)
        scenario = [scenario; children]
        push!(groups, scenario)
    end

    return groups
end

function get_scenarios(tree::AbstractTree)
    scenarios = get_groups(tree, combine = 0)

    trees = AbstractTree[]

    for i in 1:length(scenarios)
        t = Leaf()
        t.ID = scenarios[i][1]
        t.name = t.ID.name
        for j in 2:length(scenarios[i])
            t = Tree([t])
            t.ID = scenarios[i][j]
            t.name = t.ID.name
        end
        push!(trees, t)
    end

    return trees
end

mutable struct Node
    children::Vector{Node}
    pr::Float64
    name::String
end

function save_tree_to_file(tree::AbstractTree, filename::String)
    nodes = collect(tree, order = :breadth)

    f = open(filename, "w")
    println(f, "n,p")
    for n in nodes
        if parent(n) == nothing
            println(f, n.name * "," * "-")
        else
            println(f, n.name * "," * n.parent.name)
        end
    end
    return close(f)
end

"""
	tree_from_leaves(leafnodes::Vector{Vector{Int}}, probs::Vector{Float64})

Construct tree from Array of leaf nodes, and (optionally) the corresponding probabilities

### Required Arguments
`leafnodes` is an array of arrays defining the set of leaf nodes

### Optional Arguments
`probs` is an array of probabilities for the leaf nodes

### Example
    (tree,prob) = tree_from_leaves([[1,1,1],[1,1,2],[1,2,1],[1,2,2]],[0.25,0.25,0.25,0.25])
    tree = tree_from_leaves([[1,1,1],[1,1,2],[1,2,1],[1,2,2]])
"""#
function tree_from_leaves(
    leafnodes::Vector{Vector{Int}},
    probs::Vector{Float64},
)
    prob = Dict{AbstractTree,Float64}()

    function groupnode(node::Node)
        if length(node.children) == 0
            output = Leaf()
        else
            v = Vector{AbstractTree}()
            for i in 1:length(node.children)
                push!(v, groupnode(node.children[i]))
                node.pr += node.children[i].pr
            end
            output = Tree(v)
            for c in output.children
                c.parent = output
            end
        end
        prob[output] = node.pr
        return output
    end

    root = Node(Vector{Node}(), 0.0, "")
    for i in 1:length(leafnodes)
        if leafnodes[i][1] != 1
            error("there must be a unique root node")
        end
        parent = root
        for j in 2:length(leafnodes[i])
            while leafnodes[i][j] > length(parent.children)
                n = Node(Vector{Node}(), 0.0, "")
                push!(parent.children, n)
            end
            parent = parent.children[leafnodes[i][j]]
            if j == length(leafnodes[i])
                parent.pr = probs[i]
            end
        end
    end
    tree = groupnode(root)
    label_nodes(tree)
    return (tree, prob)
end

# Construct tree from Array of leaf nodes, without probabilities
function tree_from_leaves(leafnodes::Vector{Vector{Int}})
    return tree_from_leaves(leafnodes, zeros(length(leafnodes)))[1]
end

# Construct tree from nested Vector, the first element of each vector is the
# probability of a node, remaining elements are vectors representing child nodes
function tree_from_nodes(nodes::Vector{Any})
    prob = Dict{AbstractTree,Float64}()

    function groupnode(node::Any, prob, prev)
        if length(node) == 1
            output = Leaf()
        else
            v = Vector{AbstractTree}()
            for i in 2:length(node)
                push!(v, groupnode(node[i], prob, prev * node[1]))
            end
            output = Tree(v)
            for c in output.children
                c.parent = output
            end
        end
        prob[output] = node[1] * prev
        return output
    end
    tree = groupnode(nodes, prob, 1.0)
    label_nodes(tree)
    return (tree, prob)
end

"""
	tree_from_file(filename::String)

Construct tree from a file, each line in the file is of the form B,A,...
representing an arc in the tree, from node "A" to node "B". The total
number of columns is arbitrary. The first row of the file should be
n,p,... these column headers are converted into symbols used to index
the `data`. Each column itself converted into a dictionary, indexed by
the node.

### Required Arguments
`string` is the full path of the file containing the tree

### Example
    tree, data = tree_from_file(joinpath(@__DIR__,"tree.csv"))
"""#
function tree_from_file(filename::String)
    data = Dict{Symbol,Dict{Node,Float64}}()
    data2 = Dict{Symbol,Dict{AbstractTree,Float64}}()
    nodes = Dict{String,Node}()
    count = Dict{Node,Int}()
    first = true
    f = open(filename)
    headers = Vector{Symbol}()
    for l in eachline(f)
        if first
            first = false
            a = split(l, ",")
            if a[1] != "n" || a[2] != "p"
                error("First row of tree file should be \"n,p,...\"")
            end
            for i in 3:length(a)
                push!(headers, Symbol(a[i]))
                data[Symbol(a[i])] = Dict{Node,Float64}()
                data2[Symbol(a[i])] = Dict{AbstractTree,Float64}()
            end
        else
            a = split(l, ",")
            if length(a) == 0
                @warn("Found blank line")
            elseif length(a) != length(headers) + 2
                @warn("Found line with wrong number of elements: \"" * l * "\"")
                continue
            end

            if !((string)(a[2]) in keys(nodes))
                if (string)(a[2]) != "-"
                    nodes[(string)(a[2])] =
                        Node(Vector{Node}(), 1.0, (string)(a[2]))
                    count[nodes[(string)(a[2])]] = 0
                end
            end
            if !((string)(a[1]) in keys(nodes))
                nodes[(string)(a[1])] =
                    Node(Vector{Node}(), 1.0, (string)(a[1]))
                count[nodes[(string)(a[1])]] = 0
            end
            if (string)(a[2]) != "-"
                push!(nodes[(string)(a[2])].children, nodes[(string)(a[1])])
            end
            for i in 3:length(a)
                data[headers[i-2]][nodes[(string)(a[1])]] = parse(Float64, a[i])
            end
        end
    end

    close(f)

    for (nn, n) in nodes
        for i in n.children
            count[i] += 1
            if count[i] > 1
                error("Node " * i.name * " defined multiple times")
            end
        end
    end

    found = 0
    root = nothing
    for n in keys(count)
        if count[n] == 0
            root = n
            found += 1
        end
    end

    if found == 0
        error("No root node found")
    elseif found > 1
        error("Multiple root nodes found")
    end

    function groupnode(node::Node)
        if length(node.children) == 0
            output = Leaf(node.name)
        else
            v = Vector{AbstractTree}()
            for n in node.children
                push!(v, groupnode(n))
            end
            output = Tree(v)
            output.name = node.name
            for c in output.children
                c.parent = output
            end
        end
        for h in headers
            data2[h][output] = data[h][node]
        end
        return output
    end

    tree = groupnode(root)
    return (tree, data2)
end
