using JuDGE

function TreeA()
    leafnodesA = [
        [1, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 2, 1, 1, 1],
        [1, 2, 1, 1, 2],
        [1, 2, 1, 2],
        [1, 2, 1, 3],
    ]

    probsA = [0.2, 0.3, 0.1, 0.2, 0.05, 0.15]

    treeA, probabilityA = tree_from_leaves(leafnodesA, probsA)

    JuDGE.print_tree(treeA)

    nodesA = collect(treeA)

    return probabilityA[nodesA[10]] == 0.1 &&
           probabilityA[treeA[leafnodesA[3]]] == 0.1 &&
           probabilityA[treeA[1, 2, 1, 1, 1]] == 0.1
end

function TreeB()
    nodeprobsB = [
        1.0,
        [0.5, [0.4, [1.0]], [0.6, [1.0]]],
        [0.5, [1.0, [0.6, [1 / 3], [2 / 3]], [0.1], [0.3]]],
    ]

    treeB, probabilityB = tree_from_nodes(nodeprobsB)

    JuDGE.print_tree(treeB)

    nodesB = collect(treeB)
    return probabilityB[treeB.children[2].children[1].children[1].children[1]] ≈
           0.1 && probabilityB[nodesB[10]] ≈ 0.1
end

function TreeC(; visualise = false)
    treeC, data = tree_from_file(joinpath(@__DIR__, "treeC.tree"))

    nodesC = collect(treeC)
    JuDGE.print_tree(treeC, data[:pr])
    probabilityC = JuDGE.convert_probabilities(treeC, data[:pr])
    JuDGE.print_tree(treeC, probabilityC)

    if visualise
        vis_data = Dict{AbstractTree,Dict{Symbol,Any}}()
        JuDGE.add_to_dictionary!(vis_data, data[:pr], :pr_conditional)
        JuDGE.add_to_dictionary!(vis_data, probabilityC, :pr_absolute)
        JuDGE.visualize_tree(
            treeC,
            vis_data,
            scale_edges = 0.9,
            filename = "treeC",
            rel_angle = true,
        )
    end

    return probabilityC[treeC[1, 2, 1, 1, 1]] ≈ 0.1
end

function TreeD()
    tree = narytree(3, 3)
    scenarios = JuDGE.get_scenarios(tree)
    return scenarios[16].children[1] != tree.children[2] &&
           JuDGE.getID(scenarios[16].children[1]) ==
           JuDGE.getID(tree.children[2])
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    TreeA()
    TreeB()
    TreeC(visualise = true)
    TreeD()
end
