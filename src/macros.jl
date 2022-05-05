macro judge_var(model, variable, class, aargs, aakws)
    lag = 0
    span = 1000
    initial = 0.0
    ub = nothing
    lb = nothing
    penalty = nothing
    state_name = :nothing

    vartype = :Con

    if length(aargs) == 0
        push!(aargs, :Con)
    end

    if length(aargs) != 1
        ex = quote
            error(
                "Master + subproblem variable macro takes at most three positional arguments",
            )
        end
        return ex
    elseif aargs[1] ∉ [:Con, :Bin, :Int]
        ex = quote
            error("Optional third positional argument must be \'Bin\' or \'Int\'")
        end
        return ex
    end

    for (a, b) in aakws
        if a == :lag
            lag = b
        elseif a == :duration
            span = b
        elseif a == :initial
            initial = b
        elseif a == :lb
            lb = b
        elseif a == :ub
            ub = b
        elseif a == :penalty
            if typeof(b) <: Array &&
               length(b) == 2 &&
               b[1] >= 0 &&
               b[2] >= 0 &&
               b[1] + b[2] > 0
                penalty = b
            else
                error("Invalid penality vector.")
            end
        elseif a == :state_name
            state_name = b
        else
            ex = quote
                error(
                    "Invalid keyword '" *
                    string($(Meta.quot(a))) *
                    "' for master + subproblem variable macro",
                )
            end
            return ex
        end
    end

    if class == :(:state)
        if lag != 0
            @warn("'lag' keyword has been ignored")
            lag = 0
        end
        if span != 1000
            @warn("'duration' keyword has been ignored")
        end
        span = 1
    elseif state_name != :nothing
        @warn("'state_name' keyword has been ignored")
        state_name = :nothing
    end

    tmp = nothing

    ex = quote
        if !haskey($model.ext, :expansions)
            $model.ext[:expansions] = Dict{Symbol,Any}()
            $model.ext[:options] = Dict{Symbol,Tuple}()
        end
        if length($aargs) == 1
            if $aargs[1] == :Con
                tmp = @variable($model, $variable)
            elseif $aargs[1] == :Bin
                tmp = @variable($model, $variable, Bin)
            else
                tmp = @variable($model, $variable, Int)
            end
        end

        if $(Meta.quot(state_name)) == :nothing
            sym = [k for (k, v) in $model.obj_dict if v === tmp][1]
            $model.ext[:expansions][sym] = tmp
            $model.ext[:options][sym] = (
                $class,
                $lag,
                $span,
                $aargs[1],
                $lb,
                $ub,
                $initial,
                $penalty,
            )
        else
            $model.ext[:expansions][$(Meta.quot(state_name))] = tmp
            $model.ext[:options][$(Meta.quot(state_name))] = (
                $class,
                $lag,
                $span,
                $aargs[1],
                $lb,
                $ub,
                $initial,
                $penalty,
            )
            $state_name = tmp
        end
    end
    return esc(ex)
end

"""
	expansion(model, variable, args...)

Defines an expansion variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of expansion variables.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the expansion variable to

`variable` is the name of the variable being created, this will be continuous by default; follows JuMP syntax if defining a set of variables.

### Optional Arguments
This macro has a third, unnamed, argument which can be set to Con, Bin, or Int, similar to the `@variable` macro.

`lag` is the number of nodes in the scenario between an expansion being decided, and it becoming available.

`duration` is the number of consecutive nodes in the scenario over which an expansion is available.

`lb` is the lower bound for this variable in the master problem (typically omitted).

`ub` is the upper bound for this variable in the master problem (typically omitted).

### Examples
    @expansion(model, expand[1:5], Bin) #defines an array of 5 binary variables with no lag, and unlimited lifespan
    @expansion(model, expand[1:5,1:2]>=0, lag=1) #defines a matrix of 10 continuous variables with a lag of 1, and unlimited duration
    @expansion(model, 0<=expand<=10, Int, duration=2) #defines a single integer variable with a lag of 0, and a duration of 2
"""
macro expansion(model, variable, args...)
    aargs = []
    aakws = Pair{Symbol,Any}[]
    custom_penalty = false
    for el in args
        if Meta.isexpr(el, :(=))
            push!(aakws, Pair(el.args...))
            if el.args[1] == :penalty
                custom_penalty = true
            end
        else
            push!(aargs, el)
        end
    end

    if !custom_penalty
        push!(aakws, :penalty => [0, Inf])
    end

    ex = quote
        JuDGE.@judge_var($model, $variable, :cumulative, $aargs, $aakws)
    end

    return esc(ex)
end

"""
	shutdown(model, variable, args...)

Defines an shutdown variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of shutdown variables.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the shutdown variable to

`variable` is the name of the variable being created, this will be continuous by default; follows JuMP syntax if defining a set of variables.

### Optional Arguments
This macro has a third, unnamed, argument which can be set to Con, Bin, or Int, similar to the `@variable` macro.

`lag` is the number of nodes in the scenario between an shutdown being decided, and it becoming unavailable.

`duration` is the number of consecutive nodes in the scenario over which the shutdown will last.

`lb` is the lower bound for this variable in the master problem (typically omitted).

`ub` is the upper bound for this variable in the master problem (typically omitted).

### Examples
    @shutdown(model, shut[1:5], Bin) #defines an array of 5 binary variables with no lag, and unlimited lifespan
    @shutdown(model, shut[1:5,1:2]>=0, lag=1) #defines a matrix of 10 continuous variables with a lag of 1, and unlimited duration
    @shutdown(model, 0<=shut<=10, Int, duration=2) #defines a single integer variable with a lag of 0, and a duration of 2
"""
macro shutdown(model, variable, args...)
    aargs = []
    aakws = Pair{Symbol,Any}[]
    custom_penalty = false
    for el in args
        if Meta.isexpr(el, :(=))
            push!(aakws, Pair(el.args...))
            if el.args[1] == :penalty
                custom_penalty = true
            end
        else
            push!(aargs, el)
        end
    end

    if !custom_penalty
        push!(aakws, :penalty => [Inf, 0])
    end

    ex = quote
        JuDGE.@judge_var($model, $variable, :cumulative, $aargs, $aakws)
    end

    return esc(ex)
end

"""
	enforced(model, variable, args...)

Defines an enforced variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of enforced variables. These
variables can be used as either expansion or shutdown variables, but since the constraint in the master problem is an equality, convergence is
more difficult since there is less flexibility when solving the master problem.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the expansion variable to

`variable` is the name of the variable being created, this will be continuous by default; follows JuMP syntax if defining a set of variables.

### Optional Arguments
This macro has a third, unnamed, argument which can be set to Con, Bin, or Int, similar to the `@variable` macro.

`lag` is the number of nodes in the scenario between an expansion being decided, and it becoming available.

`duration` is the number of consecutive nodes in the scenario over which an expansion is available.

`lb` is the lower bound for this variable in the master problem (typically omitted).

`ub` is the upper bound for this variable in the master problem (typically omitted).

`penalty` is a placeholder for a future feature, which may allow the violation of master/subproblem equality constraint, at a cost.

### Examples
    @expansion(model, forced[1:5], Bin) #defines an array of 5 binary variables with no lag, and unlimited lifespan
    @expansion(model, forced[1:5,1:2]>=0, lag=1) #defines a matrix of 10 continuous variables with a lag of 1, and unlimited duration
    @expansion(model, 0<=forced<=10, Int, duration=2) #defines a single integer variable with a lag of 0, and a duration of 2
"""
macro enforced(model, variable, args...)
    aargs = []
    aakws = Pair{Symbol,Any}[]
    custom_penalty = false
    for el in args
        if Meta.isexpr(el, :(=))
            push!(aakws, Pair(el.args...))
            if el.args[1] == :penalty
                custom_penalty = true
            end
        else
            push!(aargs, el)
        end
    end

    if !custom_penalty
        push!(aakws, :penalty => [Inf, Inf])
    end

    ex = quote
        JuDGE.@judge_var($model, $variable, :cumulative, $aargs, $aakws)
    end

    return esc(ex)
end

"""
	state(model, variable, args...)

Defines a state variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of state variables. These
variables can be used to model inventory that is carried forward between the subproblems.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the expansion variable to

`variable` is the name of the variable being created in the subproblem, this will be continuous by default; follows JuMP syntax if defining a set of variables.
The subproblem variable corresponds to the change in the state.

### Optional Arguments
This macro has a third, unnamed, argument which can be set to Con, Bin, or Int, similar to the `@variable` macro.

`state_name` is the name for the state variable in the master problem. If omitted, the name of the master problem variable
will match the subproblem variable (but this may cause confusion, since only the master problem variable is the state).
See the `inventory.jl` example to see how this should be implemented.

`initial` is the initial value for the master problem's state variable at the root node.

`lb` is the lower bound for the variable in the master problem (typically omitted).

`ub` is the upper bound for the variable in the master problem (typically omitted).

### Examples
    @state(sp, -50<=Δstock<=50, state_name=stock, lb=0, ub=200, initial=0) #defines a state variable called stock in the master
                                                                           #(starting at 0, and able to take values 0 to 200),
                                                                           #and Δstock in the subproblem (able to change the stock level by ±50).
"""
macro state(model, variable, args...)
    aargs = []
    aakws = Pair{Symbol,Any}[]
    custom_penalty = false
    free_disposal = false
    for el in args
        if Meta.isexpr(el, :(=))
            if el.args[1] == :penalty
                custom_penalty = true
                push!(aakws, Pair(el.args[1], eval(el.args[2])))
                continue
            elseif el.args[1] == :free_disposal
                if el.args[2]
                    free_disposal = true
                end
                continue
            end
            push!(aakws, Pair(el.args...))
        else
            push!(aargs, el)
        end
    end

    if !custom_penalty
        if free_disposal
            push!(aakws, :penalty => [Inf, 0])
        else
            push!(aakws, :penalty => [Inf, Inf])
        end
    elseif free_disposal
        error("Cannot set both free_disposal and penalty for @state variable.")
    end

    ex = quote
        JuDGE.@judge_var($model, $variable, :state, $aargs, $aakws)
    end

    return esc(ex)
end

"""
	capitalcosts(model, expr)

Defines a linear expression specifying the capital cost of expansions and shutdowns at the current node

### Required Arguments
`model` is the JuDGE subproblem corresponding to the node in the scenario tree that we are adding specifying the costs for

`expr` is an `AffExpr` which gives the total cost of choosing expansion and shutdown variables at the current node

### Example
    @capitalcosts(model, sum(expand[i]*cost[node][i] for i in 1:5))
"""
macro capitalcosts(model, expr)
    ex = quote
        #$model.ext[:capitalcosts] = @expression($model, $expr)
        $model.ext[:capitalcosts] = Dict()
        if typeof($expr) == AffExpr
            for (term, coef) in $expr.terms
                $model.ext[:capitalcosts][term] = coef
            end
            $model.ext[:capitalcosts][:constant] = $expr.constant
        elseif typeof($expr) == VariableRef
            $model.ext[:capitalcosts][$expr] = 1.0
            $model.ext[:capitalcosts][:constant] = 0.0
        elseif typeof($expr) == Float64
            $model.ext[:capitalcosts][:constant] = $expr
        end
    end
    return esc(ex)
end

"""
	ongoingcosts(model, expr)

Defines a linear expression specifying the ongoing costs of expansions and shutdowns available at the current node

### Required Arguments
`model` is the JuDGE subproblem corresponding to the node in the scenario tree that we are adding specifying the costs for

`expr` is an `AffExpr` which gives the ongoing cost of expansions and shutdowns available at the current node

### Example
    @ongoingcosts(model, sum(expand[i]*ongoingcosts[node][i] for i in 1:5))
"""
macro ongoingcosts(model, expr)
    ex = quote
        #      $model.ext[:ongoingcosts] = @expression($model, $expr)
        $model.ext[:ongoingcosts] = Dict()
        if typeof($expr) == AffExpr
            for (term, coef) in $expr.terms
                $model.ext[:ongoingcosts][term] = coef
            end
            $model.ext[:ongoingcosts][:constant] = $expr.constant
        elseif typeof($expr) == VariableRef
            $model.ext[:ongoingcosts][$expr] = 1.0
            $model.ext[:ongoingcosts][:constant] = 0.0
        elseif typeof($expr) == Float64
            $model.ext[:ongoingcosts][:constant] = $expr
        end
    end
    return esc(ex)
end
