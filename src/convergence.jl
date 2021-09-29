struct ConvergenceState
    obj::Float64
    ub::Float64
    lb::Float64
    int_abs::Float64
    int_rel::Float64
    rlx_abs::Float64
    rlx_rel::Float64
    time::Float64
    iter::Int
    num_frac::Int
    function ConvergenceState(
        obj::Float64,
        ub::Float64,
        lb::Float64,
        int_abs::Float64,
        int_rel::Float64,
        rlx_abs::Float64,
        rlx_rel::Float64,
        time::Float64,
        iter::Int,
        num_frac::Int,
    )
        return new(
            obj,
            ub,
            lb,
            int_abs,
            int_rel,
            rlx_abs,
            rlx_rel,
            time,
            iter,
            num_frac,
        )
    end
    function ConvergenceState(
        obj::Float64,
        ub::Float64,
        lb::Float64,
        time::Float64,
        iter::Int,
        num_frac::Int,
    )
        return new(
            obj,
            ub,
            lb,
            ub - lb,
            (ub - lb) / abs(lb),
            obj - lb,
            (obj - lb) / abs(lb),
            time,
            iter,
            num_frac,
        )
    end
end

"""
Termination(;abstol::Float64=10^-9,
            reltol::Float64=10^-9,
            rlx_abstol::Float64=10^-9,
            rlx_reltol::Float64=10^-9,
            time_limit::Float64=Inf,
            max_iter::Int=typemax(Int),
            inttol::Float64=10^-8,
            allow_frac::Symbol=:binary_solve)

Define the stopping conditions for `JuDGE.solve()` / `JuDGE.branch_and_price()`.

### Optional Arguments
`abstol` is the absolute tolerance for the best integer-feasible objective value and the lower bound.

`reltol` is the relative tolerance for the best integer-feasible objective value and the lower bound.

`rlx_abstol` is the absolute tolerance for the relaxed master objective value and the lower bound.

`rlx_reltol` is the relative tolerance for the relaxed master objective value and the lower bound.

`time_limit` is the maximum duration in seconds.

`max_iter` is the maximum number of iterations.

`inttol` is the maximum deviation from 0 or 1 for any binary/integer variable for integer feasible solutions.

`allow_frac` indicates whether a fractional solution will be returned; possible values are:
    `:binary_solve` a binary solve of master will be performed (if needed) prior to the solution being returned;
    `:binary_solve_return_relaxation` a binary solve of master will be performed (if needed), updating the upper bound,
    but the master problem relation will be returned;
    `:first_fractional` will return the first fractional master solution found;
    `:no_binary_solve` will simply return the solution to the relaxed master problem when terminated.

### Examples
   Termination(rlx_abstol=10^-6)
   Termination(abstol=10^-3,inttol=10^-6)
"""
mutable struct Termination
    abstol::Float64
    reltol::Float64
    rlx_abstol::Float64
    rlx_reltol::Float64
    time_limit::Float64
    max_iter::Int
    inttol::Float64
    allow_frac::Symbol

    function Termination(;
        abstol::Float64 = 10^-9,
        reltol::Float64 = 10^-9,
        rlx_abstol::Float64 = 10^-9,
        rlx_reltol::Float64 = 10^-9,
        time_limit::Float64 = Inf,
        max_iter::Int = typemax(Int),
        inttol::Float64 = 10^-8,
        allow_frac::Symbol = :binary_solve,
    )
        if allow_frac ∉ [
            :binary_solve,
            :binary_solve_return_relaxation,
            :first_fractional,
            :no_binary_solve,
        ]
            error(
                "'allow_frac' can take values: :binary_solve, :binary_solve_return_relaxation, :first_fractional or :no_binary_solve",
            )
        end
        return new(
            abstol,
            reltol,
            rlx_abstol,
            rlx_reltol,
            time_limit,
            max_iter,
            inttol,
            allow_frac,
        )
    end
end

function display(termination::Termination)
    println("---------------------------------------------")
    println("Termination          Absolute       Relative")
    println("---------------------------------------------")
    Printf.@printf(
        "Binary/Integer:     %5.4e     %5.4e \n",
        termination.abstol,
        termination.reltol
    )
    Printf.@printf(
        "Relaxation:         %5.4e     %5.4e \n",
        termination.rlx_abstol,
        termination.rlx_reltol
    )
    println("---------------------------------------------")
    Printf.@printf("Integer tolerance:  %5.4e\n", termination.inttol)
    Printf.@printf("Time-limit:         %5.4e\n", termination.time_limit)
    if termination.max_iter == typemax(Int)
        print("Max iterations:     Unlimited\n")
    else
        Printf.@printf("Max iterations:     %10d\n", termination.max_iter)
    end

    Printf.@printf("Allow fractional:   %s\n", string(termination.allow_frac))
    return println("---------------------------------------------")
end

function display(cs::ConvergenceState; relaxation::Bool = true)
    print("")
    if cs.obj == Inf
        print("          ")
    elseif cs.obj >= 0
        print(" ")
    end
    Printf.@printf(" %e  | ", cs.obj)

    if cs.ub == Inf
        print("          ")
    elseif cs.ub >= 0
        print(" ")
    end
    Printf.@printf("%e ", cs.ub)

    if cs.lb == -Inf
        print("         ")
    elseif cs.lb >= 0
        print(" ")
    end
    Printf.@printf("%e  |  ", cs.lb)

    if relaxation
        temp1 = cs.rlx_abs
        temp2 = cs.rlx_rel
    else
        temp1 = cs.int_abs
        temp2 = cs.int_rel
    end

    if temp1 == Inf
        print("          ")
    elseif temp1 >= 0
        print(" ")
    end
    Printf.@printf("%e   ", temp1)

    if temp2 == Inf || isnan(temp2)
        print("          ")
    elseif temp2 >= 0
        print(" ")
    end
    Printf.@printf("%e  |   ", temp2)

    Printf.@printf("%9d  | %9.3f  %7d", cs.num_frac, cs.time, cs.iter)
    if relaxation
        println("")
    else
        println("*")
    end
end

function InitialConvergenceState()
    return ConvergenceState(Inf, Inf, -Inf, Inf, Inf, Inf, Inf, 0.0, 0, 0)
end

function has_converged(b::Termination, a::ConvergenceState)
    if a.int_abs <= b.abstol ||
       a.int_rel <= b.reltol ||
       a.rlx_abs <= b.rlx_abstol ||
       a.rlx_rel <= b.rlx_reltol ||
       a.time > b.time_limit ||
       a.iter >= b.max_iter
        true
    else
        false
    end
end
