module TestSDPNAL

using Test
using MathOptInterface
import SDPNAL

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_solver_name()
    @test MOI.get(SDPNAL.Optimizer(), MOI.SolverName()) == "SDPNAL"
end

function test_supports_incremental_interface()
    @test MOI.supports_incremental_interface(SDPNAL.Optimizer())
end

function test_options()
    optimizer = SDPNAL.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("printlevel"), 1)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("printlevel")) == 1

    param = MOI.RawOptimizerAttribute("bad_option")
    err = MOI.UnsupportedAttribute(param)
    @test_throws err MOI.set(optimizer, MOI.RawOptimizerAttribute("bad_option"), 1)
end

function test_runtests()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(
            MOI.Utilities.CachingOptimizer(
                MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
                SDPNAL.Optimizer(),
            ),
            Float64,
        ),
        # This does not work as with some modifications, the bridges with try
        # getting `ConstraintFunction` which is not supported by SDPNAL
        #MOI.instantiate(SDPNAL.Optimizer, with_bridge_type=Float64),
    )
    # `Variable.ZerosBridge` makes dual needed by some tests fail.
    MOI.Bridges.remove_bridge(model.optimizer, MathOptInterface.Bridges.Variable.ZerosBridge{Float64})
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            rtol = 1e-1,
            atol = 1e-1,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ObjectiveBound,
                MOI.SolverVersion,
            ],
        ),
        exclude = String[
            # Unable to bridge RotatedSecondOrderCone to PSD because the dimension is too small: got 2, expected >= 3.
            "test_conic_SecondOrderCone_INFEASIBLE",
            "test_constraint_PrimalStart_DualStart_SecondOrderCone",
            # Error using mexMatvec
            # mexMatvec: dimension of 4TH input is incorrect
            #
            # Error in matvecSSN
            #
            # Error in psqmrAB
            #
            # Error in SSNfun
            #
            # Error in sdpnalplus_main_default
            #
            # Error in sdpnalplus
            #
            # >>
            # Intel MKL ERROR: Parameter 6 was incorrect on entry to DGEMV .
            # Error using save
            # Variable 'jx_sdpnalplus_arg_out_1' not found.
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
            "test_conic_linear_INFEASIBLE",
            "test_conic_linear_INFEASIBLE_2",
            "test_infeasible_MAX_SENSE",
            "test_infeasible_MAX_SENSE_offset",
            "test_infeasible_MIN_SENSE",
            "test_infeasible_MIN_SENSE_offset",
            "test_infeasible_affine_MAX_SENSE",
            "test_infeasible_affine_MAX_SENSE_offset",
            "test_infeasible_affine_MIN_SENSE",
            "test_infeasible_affine_MIN_SENSE_offset",
            "test_linear_DUAL_INFEASIBLE",
            "test_linear_DUAL_INFEASIBLE_2",
            "test_linear_FEASIBILITY_SENSE",
            "test_linear_INFEASIBLE",
            "test_linear_INFEASIBLE_2",
            "test_linear_add_constraints",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_GreaterThan",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_lower",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_LessThan",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_VariableIndex_LessThan",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_VariableIndex_LessThan_max",
            "test_unbounded_MAX_SENSE",
            "test_unbounded_MAX_SENSE_offset",
            "test_unbounded_MIN_SENSE",
            "test_unbounded_MIN_SENSE_offset",
            # Incorrect result
            "test_conic_GeometricMeanCone_VectorOfVariables",
            "test_objective_FEASIBILITY_SENSE_clears_objective",
            "test_objective_qp_ObjectiveFunction_edge_cases",
            "test_objective_qp_ObjectiveFunction_zero_ofdiag",
            "test_variable_solve_with_lowerbound",
            #  Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
            #   Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
            "test_conic_GeometricMeanCone_VectorAffineFunction",
            "test_conic_GeometricMeanCone_VectorAffineFunction_3",
            "test_conic_GeometricMeanCone_VectorOfVariables_3",
            "test_conic_NormInfinityCone_3",
            "test_conic_NormOneCone",
            "test_conic_RootDetConeTriangle_VectorAffineFunction",
            "test_conic_RootDetConeTriangle_VectorOfVariables",
            "test_conic_RotatedSecondOrderCone_INFEASIBLE",
            "test_conic_SecondOrderCone_Nonnegatives",
            "test_conic_SecondOrderCone_Nonpositives",
            "test_conic_SecondOrderCone_negative_post_bound_2",
            "test_conic_SecondOrderCone_negative_post_bound_3",
            "test_conic_SecondOrderCone_no_initial_bound",
            "test_conic_linear_VectorAffineFunction",
            "test_constraint_ScalarAffineFunction_LessThan",
            "test_constraint_ScalarAffineFunction_duplicate",
            "test_constraint_VectorAffineFunction_duplicate",
            "test_linear_LessThan_and_GreaterThan",
            "test_linear_VectorAffineFunction",
            "test_linear_inactive_bounds",
            "test_linear_integration_delete_variables",
            "test_linear_integration_modification",
            "test_linear_modify_GreaterThan_and_LessThan_constraints",
            "test_linear_transform",
            "test_modification_affine_deletion_edge_cases",
            "test_modification_coef_scalaraffine_lessthan",
            "test_modification_const_vectoraffine_nonpos",
            "test_modification_delete_variables_in_a_batch",
            "test_modification_func_scalaraffine_lessthan",
            "test_modification_func_vectoraffine_nonneg",
            "test_modification_multirow_vectoraffine_nonpos",
            "test_quadratic_nonhomogeneous",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_EqualTo_lower",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_EqualTo_upper",
            "test_variable_solve_with_upperbound",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
        ],
    )
    return
end

end  # module

TestSDPNAL.runtests()
