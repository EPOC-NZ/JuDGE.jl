import HiGHS
using JuMP

JuDGE_SP_Solver = optimizer_with_attributes(
    HiGHS.Optimizer,
    "output_flag" => false,
    "mip_rel_gap" => 0.0,
)

JuDGE_DE_Solver = optimizer_with_attributes(
    HiGHS.Optimizer,
    "output_flag" => true,
    "mip_rel_gap" => 0.0,
)

JuDGE_MP_Solver = (
    optimizer_with_attributes(
        HiGHS.Optimizer,
        #    "solver" => "ipm", HiGHS' interior point method solver doesn't always work
        "output_flag" => false,
        "mip_rel_gap" => 0.0,
    ),
    JuDGE_SP_Solver,
)
