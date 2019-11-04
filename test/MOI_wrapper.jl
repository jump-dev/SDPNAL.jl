using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SDPNAL
const OPTIMIZER = SDPNAL.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "SDPNAL"
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(OPTIMIZER, false)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, true)
end

# UniversalFallback is needed for starting values, even if they are ignored by SDPNAL
const CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
const CACHED = MOIU.CachingOptimizer(CACHE, OPTIMIZER)
const BRIDGED = MOIB.full_bridge_optimizer(CACHED, Float64)
const CONFIG = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Options" begin
    optimizer = SDPNAL.Optimizer(printlevel = 1)
    @test MOI.get(optimizer, MOI.RawParameter("printlevel")) == 1

    param = MOI.RawParameter("bad_option")
    err = MOI.UnsupportedAttribute(param)
    @test_throws err SDPNAL.Optimizer(bad_option = 1)
end

@testset "Unit" begin
    MOIT.unittest(BRIDGED, CONFIG, [
        #  Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #   Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        #  Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #   Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        #  Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #   Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        #  Expression: MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        #   Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        "solve_with_lowerbound",
        # Unbounded problem not supported
        "solve_unbounded_model",
        # `TimeLimitSec` not supported.
        "time_limit_sec",
        # `NumberOfThreads` not supported.
        "number_threads",
        # Quadratic functions are not supported
        "solve_qcp_edge_cases", "solve_qp_edge_cases",
        # Integer and ZeroOne sets are not supported
        "solve_integer_edge_cases", "solve_objbound_edge_cases",
        "solve_zero_one_with_bounds_1",
        "solve_zero_one_with_bounds_2",
        "solve_zero_one_with_bounds_3"])
end
@testset "Continuous Linear" begin
    # See explanation in `MOI/test/Bridges/lazy_bridge_OPTIMIZER.jl`.
    # This is to avoid `Variable.VectorizeBridge` which does not support
    # `ConstraintSet` modification.
    MOIB.remove_bridge(BRIDGED, MOIB.Constraint.ScalarSlackBridge{Float64})
    MOIT.contlineartest(BRIDGED, CONFIG, String[
        # Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #  Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        # Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        # Expression: ≈(MOI.get(model, MOI.ObjectiveValue()), 4, atol=atol, rtol=rtol)
        #  Evaluated: 3.990384638642647 ≈ 4 (atol=0.0001, rtol=0.0001)
        # Expression: ≈(MOI.get(model, MOI.VariablePrimal(), [x, y]), [T(4), zero(T)], atol=atol, rtol=rtol)
        #  Evaluated: [3.9838363372959242, 0.006548301346723018] ≈ [4.0, 0.0] (atol=0.0001, rtol=0.0001)
        "linear5",
        # Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #  Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        # Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        # Expression: MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        "linear10b",
        # Throws error: total dimension of C should be > length(b)
        #"linear15",
        # Infeasible problems not supported
        "linear8a", "linear8b", "linear8c", "linear12",
        "partial_start"
    ])
end
@testset "Continuous Conic" begin
    MOIT.contconictest(BRIDGED, CONFIG, [
        # Infeasible problems not supported
        "lin3", "norminf2", "normone2", "rotatedsoc2",
        # Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #  Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        # Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NO_SOLUTION == MathOptInterface.FEASIBLE_POINT
        # Expression: MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NO_SOLUTION == MathOptInterface.FEASIBLE_POINT
        "soc2n", "soc2p",
        # MethodError: no method matching symvec(::Array{Float64,1})
        "soc3",
        # Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #  Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        # Expression: MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        "geomean1f", "geomean1v", "lin4",
        # `ExponentialCone` and `PowerCone` not supported.
        "exp", "dualexp", "pow", "dualpow", "logdet",
        # `RootDetConeSquare` -> `RootDetConeTriangle` bridge missing.
        "rootdets"
    ])
end
