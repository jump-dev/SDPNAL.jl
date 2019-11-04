using LinearAlgebra # For rmul!

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const AFFEQ = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Cdouble}, MOI.EqualTo{Cdouble}}

@enum VariableType NNEG PSD

struct VariableInfo
    variable_type::VariableType
    cone_index::Int
    index_in_cone::Int
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    b::Vector{Float64}

    variable_info::Vector{VariableInfo}

    # PSDCone variables
    psdc_dims::Vector{Int}
    psdc_Cvar::Vector{Vector{Int}}
    psdc_Cval::Vector{Vector{Float64}}
    psdc_Avar::Vector{Vector{Int}}
    psdc_Acon::Vector{Vector{Int}}
    psdc_Aval::Vector{Vector{Float64}}

    # NonNEGatives variables
    num_nneg::Int
    nneg_info::Vector{Int} # Similar to `info` field of `MOI.Bridges.Variable.Map`.
    nneg_Cvar::Vector{Int}
    nneg_Cval::Vector{Float64}
    nneg_Avar::Vector{Int}
    nneg_Acon::Vector{Int}
    nneg_Aval::Vector{Float64}

    objective_sense::MOI.OptimizationSense
    objective_constant::Float64

    primal_objective_value::Float64
    dual_objective_value::Float64
    nneg_X::Vector{Float64}
    psdc_X::Vector{Vector{Float64}}
    y::Vector{Float64}
    nneg_Z::Vector{Float64}
    psdc_Z::Vector{Vector{Float64}}
    info::Dict{String, Any}
    runhist::Dict{String, Any}
    status::Union{Nothing, Int}
    solve_time::Float64

    silent::Bool
    options::Dict{Symbol, Any}
    function Optimizer(; kwargs...)
        optimizer = new(
            Float64[], VariableInfo[],
            Int[], Vector{Int}[], Vector{Float64}[], Vector{Int}[], Vector{Int}[], Vector{Float64}[],
            0, Int[], Int[], Float64[], Int[], Int[], Float64[],
            MOI.FEASIBILITY_SENSE, 0.0,
            NaN, NaN,
            Float64[], Vector{Float64}[], # X
            Float64[], # y
            Float64[], Vector{Float64}[], # Z
            Dict{String, Any}(), Dict{String, Any}(),
            nothing, NaN,
            false, Dict{Symbol, Any}())
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawParameter(string(key)), value)
        end
        return optimizer
    end
end

function MOI.supports(optimizer::Optimizer, param::MOI.RawParameter)
    return param.name in ALLOWED_OPTIONS
end
function MOI.set(optimizer::Optimizer, param::MOI.RawParameter, value)
    if !MOI.supports(optimizer, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    optimizer.options[Symbol(param.name)] = value
end
function MOI.get(optimizer::Optimizer, param::MOI.RawParameter)
    # TODO: This gives a poor error message if the name of the parameter is invalid.
    return optimizer.options[Symbol(param.name)]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

MOI.get(::Optimizer, ::MOI.SolverName) = "SDPNAL"
function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
    return optimizer.solve_time
end

function MOI.is_empty(optimizer::Optimizer)
    return isempty(optimizer.b) &&
        isempty(optimizer.variable_info) &&
        iszero(optimizer.psdc_dims) &&
        iszero(optimizer.num_nneg) &&
        optimizer.objective_sense == MOI.FEASIBILITY_SENSE &&
        iszero(optimizer.objective_constant)
end

function MOI.empty!(optimizer::Optimizer)
    empty!(optimizer.b)

    empty!(optimizer.variable_info)

    empty!(optimizer.psdc_dims)
    empty!(optimizer.psdc_Cvar)
    empty!(optimizer.psdc_Cval)
    empty!(optimizer.psdc_Avar)
    empty!(optimizer.psdc_Acon)
    empty!(optimizer.psdc_Aval)

    optimizer.num_nneg = 0
    empty!(optimizer.nneg_info)
    empty!(optimizer.nneg_Cvar)
    empty!(optimizer.nneg_Cval)
    empty!(optimizer.nneg_Avar)
    empty!(optimizer.nneg_Acon)
    empty!(optimizer.nneg_Aval)

    optimizer.objective_sense = MOI.FEASIBILITY_SENSE
    optimizer.objective_constant = 0.0

    optimizer.primal_objective_value = NaN
    optimizer.dual_objective_value = NaN
    empty!(optimizer.nneg_X)
    empty!(optimizer.psdc_X)
    empty!(optimizer.y)
    empty!(optimizer.nneg_Z)
    empty!(optimizer.psdc_Z)
    empty!(optimizer.info)
    empty!(optimizer.runhist)
    optimizer.status = nothing
    optimizer.solve_time = NaN
end

function MOI.supports(
    optimizer::Optimizer,
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.VectorOfVariables}, ::Type{MOI.Reals})
    return false
end
const SupportedSets = Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.VectorOfVariables},
    ::Type{<:SupportedSets})
    return true
end
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Cdouble}},
    ::Type{MOI.EqualTo{Cdouble}})
    return true
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end
MOIU.supports_default_copy_to(::Optimizer, copy_names::Bool) = !copy_names

# Variables
function _add_nonneg_variable(optimizer::Optimizer)
    optimizer.num_nneg += 1
    push!(optimizer.variable_info, VariableInfo(NNEG, 1, optimizer.num_nneg))
    return MOI.VariableIndex(length(optimizer.variable_info))
end
function _add_psdc_variable(optimizer::Optimizer, index_in_cone)
    push!(optimizer.variable_info, VariableInfo(PSD, length(optimizer.psdc_dims), index_in_cone))
    return MOI.VariableIndex(length(optimizer.variable_info))
end
function _add_constrained_variables(optimizer::Optimizer, set::MOI.Nonnegatives)
    push!(optimizer.nneg_info, MOI.dimension(set))
    for i in 2:MOI.dimension(set)
        push!(optimizer.nneg_info, i)
    end
    return [_add_nonneg_variable(optimizer) for i in 1:MOI.dimension(set)]
end
function _add_constrained_variables(optimizer::Optimizer, set::MOI.PositiveSemidefiniteConeTriangle)
    push!(optimizer.psdc_dims, MOI.side_dimension(set))
    push!(optimizer.psdc_Cvar, Int[])
    push!(optimizer.psdc_Cval, Float64[])
    push!(optimizer.psdc_Avar, Int[])
    push!(optimizer.psdc_Acon, Int[])
    push!(optimizer.psdc_Aval, Float64[])
    return [_add_psdc_variable(optimizer, i) for i in 1:MOI.dimension(set)]
end
function MOI.add_constrained_variables(optimizer::Optimizer, set::SupportedSets)
    vis = _add_constrained_variables(optimizer, set)
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(set)}(first(vis).value)
    return vis, ci
end

# Objective
function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveSense)
    return optimizer.objective_sense
end
sense_to_sign(sense::MOI.OptimizationSense) = sense == MOI.MAX_SENSE ? -1 : 1
function MOI.set(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    if sense != optimizer.objective_sense
        sign = sense_to_sign(sense)
        rmul!(optimizer.nneg_Cval, -1)
        for i in eachindex(optimizer.psdc_dims)
            rmul!(optimizer.psdc_Cval[i], -1)
        end
    end
    optimizer.objective_sense = sense
end
function MOI.set(optimizer::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                 func::MOI.ScalarAffineFunction{Float64})
    optimizer.objective_constant = MOI.constant(func)
    empty!(optimizer.nneg_Cvar)
    empty!(optimizer.nneg_Cval)
    for i in eachindex(optimizer.psdc_dims)
        empty!(optimizer.psdc_Cvar[i])
        empty!(optimizer.psdc_Cval[i])
    end
    sign = sense_to_sign(optimizer.objective_sense)
    for term in func.terms
        info = optimizer.variable_info[term.variable_index.value]
        if info.variable_type == NNEG
            push!(optimizer.nneg_Cvar, info.index_in_cone)
            push!(optimizer.nneg_Cval, sign * term.coefficient)
        else
            @assert info.variable_type == PSD
            push!(optimizer.psdc_Cvar[info.cone_index], info.index_in_cone)
            push!(optimizer.psdc_Cval[info.cone_index], sign * term.coefficient)
        end
    end
end

# Constraints
function is_diagonal_index(k)
    # See https://www.juliaopt.org/MathOptInterface.jl/v0.9.3/apireference/#MathOptInterface.AbstractSymmetricMatrixSetTriangle
    i = div(1 + isqrt(8k - 7), 2)
    j = k - div((i - 1) * i, 2)
    return i == j
end
function MOI.add_constraint(optimizer::Optimizer, func::MOI.ScalarAffineFunction{Float64}, set::MOI.EqualTo{Float64})
    if !iszero(MOI.constant(func))
        # We use the fact that the initial function constant was zero to
        # implement getters for `MOI.ConstraintPrimal`.
        throw(MOI.ScalarFunctionConstantNotZero{
             Float64, typeof(func), typeof(set)}(MOI.constant(func)))
    end
    push!(optimizer.b, MOI.constant(set))
    con = length(optimizer.b)
    for term in func.terms
        info = optimizer.variable_info[term.variable_index.value]
        if info.variable_type == NNEG
            push!(optimizer.nneg_Avar, info.index_in_cone)
            push!(optimizer.nneg_Acon, con)
            push!(optimizer.nneg_Aval, term.coefficient)
        else
            @assert info.variable_type == PSD
            push!(optimizer.psdc_Avar[info.cone_index], info.index_in_cone)
            push!(optimizer.psdc_Acon[info.cone_index], con)
            coef = is_diagonal_index(info.index_in_cone) ? term.coefficient : term.coefficient / √2
            push!(optimizer.psdc_Aval[info.cone_index], coef)
        end
    end
    return AFFEQ(con)
end

# TODO could do something more efficient here
#      `SparseMatrixCSC` is returned in SumOfSquares.jl test `sos_horn`
symvec(Q::SparseMatrixCSC) = symvec(Matrix(Q))
function symvec(Q::Matrix)
    n = LinearAlgebra.checksquare(Q)
    vec_dim = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    q = Vector{eltype(Q)}(undef, vec_dim)
    k = 0
    for j in 1:n
        for i in 1:j
            k += 1
            q[k] = Q[i, j]
        end
    end
    @assert k == length(q)
    return q
end

function MOI.optimize!(optimizer::Optimizer)
    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:printlevel] = 0
    end

    blk1 = Any[]
    blk2 = Any[]
    m = length(optimizer.b)
    At = SparseArrays.SparseMatrixCSC{Float64,Int}[]
    # FIXME I get a strange failure with sparse vectors, need to investigate
    C = Union{Matrix{Float64}, Vector{Float64}}[]
    if !isempty(optimizer.psdc_dims)
        for (i, dim) in enumerate(optimizer.psdc_dims)
            vec_dim = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(dim))
            push!(At, sparse(optimizer.psdc_Avar[i], optimizer.psdc_Acon[i], optimizer.psdc_Aval[i], vec_dim, m))
            c = Vector(sparsevec(optimizer.psdc_Cvar[i], optimizer.psdc_Cval[i], vec_dim))
            Ci = zeros(dim, dim)
            k = 0
            for col in 1:dim
                for row in 1:(col - 1)
                    k += 1
                    Ci[row, col] = c[k] / 2
                    Ci[col, row] = c[k] / 2
                end
                k += 1
                Ci[col, col] = c[k]
            end
            push!(C, Ci)
            push!(blk1, "s")
            push!(blk2, Float64(dim))
        end
    end
    if !iszero(optimizer.num_nneg)
        push!(At, sparse(optimizer.nneg_Avar, optimizer.nneg_Acon, optimizer.nneg_Aval, optimizer.num_nneg, m))
        push!(C, Vector(sparsevec(optimizer.nneg_Cvar, optimizer.nneg_Cval, optimizer.num_nneg)))
        push!(blk1, "l")
        push!(blk2, Float64(optimizer.num_nneg))
    end
    blk = [blk1 blk2]

    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:printlevel] = 0
    end

    obj, X, s, y, Z1, Z2, y2, v, optimizer.info, optimizer.runhist = sdpnalplus(
        blk, At, C, optimizer.b; options...)

    optimizer.primal_objective_value, optimizer.dual_objective_value = obj
    k = 0
    optimizer.psdc_X = symvec.(X[k .+ eachindex(optimizer.psdc_dims)])
    optimizer.psdc_Z = symvec.(Z1[k .+ eachindex(optimizer.psdc_dims)])
    k += length(optimizer.psdc_dims)
    if iszero(optimizer.num_nneg)
        empty!(optimizer.nneg_X)
        empty!(optimizer.nneg_Z)
    else
        k += 1
        optimizer.nneg_X = X[k]
        optimizer.nneg_Z = Z1[k]
    end
    optimizer.y = y
    optimizer.status = optimizer.info["termcode"]
    optimizer.solve_time = optimizer.info["totaltime"]
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    if optimizer.status === nothing
        throw(MOI.OptimizeNotCalled())
    else
        # See `solver_main_default/sdpnalplus.m`
        if optimizer.status == 0
            return "SDP is solved to the required accuracy"
        elseif optimizer.status == 1
            return "SDP is not solved successfully"
        else
            if !(-1 <= optimizer.status <= -3)
                error("Unrecognized SDPNAL termcode `$(optimizer.status)`.")
            end
            stoptol = get(optimizer.options, "stoptol", 1e-6)
            prim = optimizer.status == -2 ? 1.05stoptol : 10stoptol
            dual = optimizer.status == -1 ? 1.05stoptol : (optimizer.status == -3 ? 10stoptol : stoptol)
            return "SDP is partially solved successfully: primfeasorg < $prim, dualfeasorg < $dual"
        end
    end
end

const TERMINATION_STATUS = [
    MOI.NUMERICAL_ERROR,
    MOI.NUMERICAL_ERROR,
    MOI.NUMERICAL_ERROR,
    MOI.OPTIMAL,
    MOI.NUMERICAL_ERROR
]
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    if optimizer.status === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    else
        return TERMINATION_STATUS[4 + optimizer.status]
    end
end

const PRIMAL_STATUS = [
    MOI.NEARLY_FEASIBLE_POINT,
    MOI.NEARLY_FEASIBLE_POINT,
    MOI.NEARLY_FEASIBLE_POINT,
    MOI.FEASIBLE_POINT,
    MOI.NO_SOLUTION
]
function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return PRIMAL_STATUS[4 + optimizer.status]
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return PRIMAL_STATUS[4 + optimizer.status]
end

MOI.get(::Optimizer, ::MOI.ResultCount) = 1
function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    sign = sense_to_sign(optimizer.objective_sense)
    return sign * optimizer.primal_objective_value + optimizer.objective_constant
end
function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    sign = sense_to_sign(optimizer.objective_sense)
    return sign * optimizer.dual_objective_value + optimizer.objective_constant
end

function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(optimizer, attr)
    info = optimizer.variable_info[vi.value]
    if info.variable_type == NNEG
        return optimizer.nneg_X[info.index_in_cone]
    else
        @assert info.variable_type == PSD
        return optimizer.psdc_X[info.cone_index][info.index_in_cone]
    end
end

function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, S}) where S<:SupportedSets
    MOI.check_result_index_bounds(optimizer, attr)
    info = optimizer.variable_info[ci.value]
    if info.variable_type == NNEG
        dim = optimizer.nneg_info[info.index_in_cone]
        return optimizer.nneg_X[(info.index_in_cone - 1) .+ (1:dim)]
    else
        @assert info.variable_type == PSD
        return optimizer.psdc_X[info.cone_index]
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal, ci::AFFEQ)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.b[ci.value]
end

function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, S}) where S<:SupportedSets
    MOI.check_result_index_bounds(optimizer, attr)
    info = optimizer.variable_info[ci.value]
    if info.variable_type == NNEG
        dim = optimizer.nneg_info[info.index_in_cone]
        return optimizer.nneg_Z[(info.index_in_cone - 1) .+ (1:dim)]
    else
        @assert info.variable_type == PSD
        return optimizer.psdc_Z[info.cone_index]
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual, ci::AFFEQ)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.y[ci.value]
end
