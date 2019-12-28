using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SDPNAL
# With the default tolerance (`tol=1e-6`), the termination status is often
# `NUMERICAL_ERROR` when there are inequality constraints.
const OPTIMIZER = SDPNAL.Optimizer(tol=1e-4)
MOI.set(OPTIMIZER, MOI.Silent(), true)
const CONFIG = MOIT.TestConfig(atol=1e-3, rtol=1e-3)

# TODO move to MOI.Test
@testset "Upper bound for $set" for set in [MOI.Nonnegatives(1), MOI.PositiveSemidefiniteConeTriangle(1)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, set)
    fx = MOI.SingleVariable(x[1])
    c = MOI.add_constraint(OPTIMIZER, fx, MOI.LessThan(1.0))
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = 1.0fx
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), c) ≈ 1.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), c) ≈ -1.0 atol = CONFIG.atol
end
@testset "Lower bound for $set" for set in [MOI.Nonnegatives(1), MOI.PositiveSemidefiniteConeTriangle(1)]
    x, cx = MOI.add_constrained_variables(OPTIMIZER, set)
    fx = MOI.SingleVariable(x[1])
    c = MOI.add_constraint(OPTIMIZER, fx, MOI.GreaterThan(1.0))
    MOI.set(OPTIMIZER, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = 1.0fx
    MOI.set(OPTIMIZER, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.optimize!(OPTIMIZER)
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), c) ≈ 1.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), c) ≈ 1.0 atol = CONFIG.atol
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
    @test MOI.get(OPTIMIZER, MOI.ObjectiveValue()) ≈ 2set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), a) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), b) ≈ 0.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), c) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), cx) ≈ zeros(3) atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_a) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_a) ≈ -1.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_b) ≈ 0.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_b) ≈ 0.0 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_c) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_c) ≈ -1.0 atol = CONFIG.atol
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
    @test MOI.get(OPTIMIZER, MOI.ObjectiveValue()) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), a) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), b) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.VariablePrimal(), c) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), cx) ≈ [1/4, -1/4, 1/4] atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_a) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_a) ≈ -1/4 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_b) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_b) ≈ -1/4 atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintPrimal(), con_c) ≈ set.upper atol = CONFIG.atol
    @test MOI.get(OPTIMIZER, MOI.ConstraintDual(), con_c) ≈ -1/4 atol = CONFIG.atol
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

@testset "Options" begin
    optimizer = SDPNAL.Optimizer(printlevel = 1)
    @test MOI.get(optimizer, MOI.RawParameter("printlevel")) == 1

    param = MOI.RawParameter("bad_option")
    err = MOI.UnsupportedAttribute(param)
    @test_throws err SDPNAL.Optimizer(bad_option = 1)
end

@testset "Unit" begin
    MOIT.unittest(BRIDGED, CONFIG, [
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
        # Error using mexMatvec
        # mexMatvec: dimension of 4TH input is incorrect
        # Error in matvecSSN
        # Error in psqmrAB
        # Error in SSNfun
        # Error in sdpnalplus_main_default
        # Error in sdpnalplus
        # >>
        # Intel MKL ERROR: Parameter 6 was incorrect on entry to DGEMV .
        # Error using save
        # Variable 'jx_sdpnalplus_arg_out_1' not found.
        "linear9",
        # Infeasible problems not supported
        "linear8a", "linear8b", "linear8c", "linear12",
        "partial_start"
    ])
end
@testset "Continuous Conic" begin
    MOIT.contconictest(BRIDGED, CONFIG, [
        # Infeasible problems not supported
        "lin3", "norminf2", "normone2", "rotatedsoc2", "soc3",
        # Error using mexMatvec
        # mexMatvec: dimension of 4TH input is incorrect
        # Error in matvecSSN
        # Error in psqmrAB
        # Error in SSNfun
        # Error in sdpnalplus_main_default
        # Error in sdpnalplus
        # >>
        # Intel MKL ERROR: Parameter 6 was incorrect on entry to DGEMV .
        # Error using save
        # Variable 'jx_sdpnalplus_arg_out_1' not found.
        "lin4",
        # >> s2      [iter=1753,(b) reset ADM.maxiter to 200]
        # s2      [iter=1966,(b) reset ADM.maxiter to 300]
        #       [iter=2281,(b) reset ADM.maxiter to 400]
        # s2s2s2s2s2s2s2s2s2s2s2s2s2s2s2s2s2      [iter=2694,(b) reset ADM.maxiter to 500]
        #       [iter=3210,(b) reset ADM.maxiter to 600]
        # s2      [iter=3824,(b) reset ADM.maxiter to 700]
        # s2s2s2s2s2s2s2s2s2      [iter=4537,(b) reset ADM.maxiter to 800]
        # Test Failed at /home/blegat/.julia/packages/MathOptInterface/izyVf/src/Test/contconic.jl:1262
        #   Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #    Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
        "geomean1f", "geomean1v",
        # `ExponentialCone` and `PowerCone` not supported.
        "exp", "dualexp", "pow", "dualpow", "logdet",
        # `RootDetConeSquare` -> `RootDetConeTriangle` bridge missing.
        "rootdets"
    ])
end
