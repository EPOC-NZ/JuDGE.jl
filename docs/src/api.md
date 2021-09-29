# API Reference

## AbstractTree Functions

### Defining Trees
```@docs
JuDGE.narytree
JuDGE.tree_from_file
JuDGE.tree_from_leaves
JuDGE.print_tree(::AbstractTree, ::Dict{AbstractTree,T} where T <: Any)
```

### Nodes of Trees
```@docs
Base.collect
JuDGE.get_leafnodes
JuDGE.get_node
```

### Tree Probabilities
```@docs
JuDGE.convert_probabilities
JuDGE.ConditionallyUniformProbabilities
JuDGE.UniformLeafProbabilities
```

### Other Tree functions
```@docs
JuDGE.depth
JuDGE.history
JuDGE.visualize_tree
JuDGE.get_groups
```

## JuDGE Functions

### JuDGE Solving Functions
```@docs
JuDGE.JuDGEModel
JuDGE.solve(::JuDGEModel)
JuDGE.branch_and_price
JuDGE.Termination
JuDGE.variable_branch
JuDGE.resolve_subproblems
JuDGE.set_policy
```

### JuDGE Macros for Subproblems
```@docs
JuDGE.@expansion
JuDGE.@shutdown
JuDGE.@enforced
JuDGE.@state
JuDGE.@capitalcosts
JuDGE.@ongoingcosts
```

### JuDGE Solutions / Output
```@docs
JuDGE.solution_to_dictionary(::JuDGEModel)
JuDGE.get_active_columns(::JuDGEModel)
JuDGE.write_solution_to_file(::JuDGEModel,::String)
JuDGE.print_expansions(::JuDGEModel)
```

## Deterministic Equivalent

### Define and Solve DetEq Model
```@docs
JuDGE.DetEqModel
JuDGE.solve(::DetEqModel)
JuDGE.set_starting_solution!
```

### DetEq Solutions / Output
```@docs
JuDGE.solution_to_dictionary(::DetEqModel)
JuDGE.write_solution_to_file(::DetEqModel,::String)
JuDGE.print_expansions(::DetEqModel)
```

## Risk
```@docs
JuDGE.RiskNeutral()
JuDGE.Risk(::Float64,::Float64;::Union{Dict{Leaf,Float64},Nothing},::Union{Float64,Nothing},::Float64)
JuDGE.Risk(::Float64;::Union{Dict{Leaf,Float64},Nothing},::Union{Float64,Nothing},::Float64)
```

## Other functions
```@docs
JuDGE.remove_from_dictionary!
JuDGE.add_to_dictionary!
```
