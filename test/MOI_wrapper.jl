using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SDPNAL
const OPTIMIZER = SDPNAL.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

# TODO move to MOI.Test
@testset "Upper bound for $set" for set in [MOI.Nonnegatives(1), MOI.PositiveSemidefiniteConeTriangle(1)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, set)
    fx = MOI.SingleVariable(x[1])
    c = MOI.add_constraint(OPTIMIZER, fx, MOI.LessThan(1.0))
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = 1.0fx
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), c) ≈ 1.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), c) ≈ -1.0 atol = 1e-5
end
@testset "Lower bound for $set" for set in [MOI.Nonnegatives(1), MOI.PositiveSemidefiniteConeTriangle(1)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, set)
    fx = MOI.SingleVariable(x[1])
    c = MOI.add_constraint(OPTIMIZER, fx, MOI.GreaterThan(1.0))
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = 1.0fx
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), c) ≈ 1.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), c) ≈ 1.0 atol = 1e-5
end
@testset "Off-diagonal bounds for $set" for set in [MOI.Interval(-2.0, 2.0), MOI.Interval(-2.0, 3.0), MOI.Interval(-3.0, 2.0)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, MOI.PositiveSemidefiniteConeTriangle(2))
    a, b, c = x
    fa = MOI.SingleVariable(a)
    fb = MOI.SingleVariable(b)
    fc = MOI.SingleVariable(c)
    con_a = MOI.add_constraint(OPTIMIZER, fa, set)
    con_b = MOI.add_constraint(OPTIMIZER, fb, set)
    con_c = MOI.add_constraint(OPTIMIZER, fc, set)
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = 1.0fa + 1.0fc
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ObjectiveValue()) ≈ 2set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), a) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), b) ≈ 0.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), c) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), cx) ≈ zeros(3) atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_a) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_a) ≈ -1.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_b) ≈ 0.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_b) ≈ 0.0 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_c) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_c) ≈ -1.0 atol = 1e-5
end
@testset "Off-diagonal bounds for $set" for set in [MOI.Interval(-2.0, 2.0), MOI.Interval(-2.0, 3.0), MOI.Interval(-3.0, 2.0)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, MOI.PositiveSemidefiniteConeTriangle(2))
    a, b, c = x
    fa = MOI.SingleVariable(a)
    fb = MOI.SingleVariable(b)
    fc = MOI.SingleVariable(c)
    con_a = MOI.add_constraint(OPTIMIZER, fa, set)
    con_b = MOI.add_constraint(OPTIMIZER, fb, set)
    con_c = MOI.add_constraint(OPTIMIZER, fc, set)
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = 1.0fb
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ObjectiveValue()) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), a) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), b) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), c) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), cx) ≈ [1/4, -1/4, 1/4] atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_a) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_a) ≈ -1/4 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_b) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_b) ≈ -1/4 atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_c) ≈ set.upper atol = 1e-5
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_c) ≈ -1/4 atol = 1e-5
end

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "SDPNAL"
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(OPTIMIZER, false)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, true)
end

# UniversalFallback is needed for starting values, even if they are ignored by SDPNAL
const CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
MOI.empty!(OPTIMIZER)
const CACHED = MOIU.CachingOptimizer(CACHE, OPTIMIZER)
const BRIDGED = MOIB.full_bridge_optimizer(CACHED, Float64)
MOIB.remove_bridge(BRIDGED, MOIB.Constraint.ScalarSlackBridge{Float64})
const CONFIG = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

#MOIT.linear9test(BRIDGED, CONFIG)

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
        "solve_with_lowerbound", "solve_affine_lessthan", "solve_affine_deletion_edge_cases", "solve_duplicate_terms_scalar_affine", "solve_with_upperbound", "solve_duplicate_terms_vector_affine",
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
    #MOIB.remove_bridge(BRIDGED, MOIB.Constraint.ScalarSlackBridge{Float64}) # Already removed above
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
