using LinearAlgebra # For rmul!

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const AFFEQ = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}
const BoundsSet = Union{MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64}}

@enum VariableType NNEG PSD

struct VariableInfo
    variable_type::VariableType
    cone_index::Int
    index_in_cone::Int
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    b::Vector{Float64}
    l::Vector{Float64}
    u::Vector{Float64}

    variable_info::Vector{VariableInfo}
    single_variable_mask::Vector{UInt8}

    # PSDCone variables
    psdc_dims::Vector{Int}
    psdc_Cvar::Vector{Vector{Int}}
    psdc_Cval::Vector{Vector{Float64}}
    psdc_Avar::Vector{Vector{Int}}
    psdc_Acon::Vector{Vector{Int}}
    psdc_Aval::Vector{Vector{Float64}}
    psdc_Bvar::Vector{Vector{Int}}
    psdc_Bcon::Vector{Vector{Int}}
    psdc_Bval::Vector{Vector{Float64}}
    psdc_L::Vector{Matrix{Float64}}
    psdc_U::Vector{Matrix{Float64}}

    # NonNEGatives variables
    num_nneg::Int
    nneg_info::Vector{Int} # Similar to `info` field of `MOI.Bridges.Variable.Map`.
    nneg_Cvar::Vector{Int}
    nneg_Cval::Vector{Float64}
    nneg_Avar::Vector{Int}
    nneg_Acon::Vector{Int}
    nneg_Aval::Vector{Float64}
    nneg_Bvar::Vector{Int}
    nneg_Bcon::Vector{Int}
    nneg_Bval::Vector{Float64}
    nneg_L::Vector{Float64}
    nneg_U::Vector{Float64}

    objective_sense::MOI.OptimizationSense
    objective_constant::Float64

    primal_objective_value::Float64
    dual_objective_value::Float64
    nneg_X::Vector{Float64}
    psdc_X::Vector{Vector{Float64}}
    s::Vector{Float64}
    y::Vector{Float64}
    y2::Vector{Float64}
    nneg_Z1::Vector{Float64}
    nneg_Z2::Vector{Float64}
    psdc_Z1::Vector{Vector{Float64}}
    psdc_Z2::Vector{Vector{Float64}}
    info::Dict{String, Any}
    runhist::Dict{String, Any}
    status::Union{Nothing, Int}
    solve_time::Float64

    silent::Bool
    options::Dict{Symbol, Any}
    function Optimizer(; kwargs...)
        optimizer = new(
            Float64[], Float64[], Float64[], VariableInfo[], UInt8[],
            Int[], Vector{Int}[], Vector{Float64}[], Vector{Int}[], Vector{Int}[], Vector{Float64}[], Vector{Int}[], Vector{Int}[], Vector{Float64}[], Matrix{Float64}[], Matrix{Float64}[],
            0, Int[], Int[], Float64[], Int[], Int[], Float64[], Int[], Int[], Float64[], Float64[], Float64[],
            MOI.FEASIBILITY_SENSE, 0.0,
            NaN, NaN,
            Float64[], Vector{Float64}[], # X
            Float64[], # s
            Float64[], Float64[], # y, y2 ≈ v
            Float64[], Float64[], Vector{Float64}[], Vector{Float64}[], # Z
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
    return isempty(optimizer.b) && isempty(optimizer.l) && isempty(optimizer.u) &&
        isempty(optimizer.variable_info) &&
        iszero(optimizer.psdc_dims) &&
        iszero(optimizer.num_nneg) &&
        optimizer.objective_sense == MOI.FEASIBILITY_SENSE &&
        iszero(optimizer.objective_constant)
end

function MOI.empty!(optimizer::Optimizer)
    empty!(optimizer.b)
    empty!(optimizer.l)
    empty!(optimizer.u)

    empty!(optimizer.variable_info)
    empty!(optimizer.single_variable_mask)

    empty!(optimizer.psdc_dims)
    empty!(optimizer.psdc_Cvar)
    empty!(optimizer.psdc_Cval)
    empty!(optimizer.psdc_Avar)
    empty!(optimizer.psdc_Acon)
    empty!(optimizer.psdc_Aval)
    empty!(optimizer.psdc_Bvar)
    empty!(optimizer.psdc_Bcon)
    empty!(optimizer.psdc_Bval)
    empty!(optimizer.psdc_L)
    empty!(optimizer.psdc_U)

    optimizer.num_nneg = 0
    empty!(optimizer.nneg_info)
    empty!(optimizer.nneg_Cvar)
    empty!(optimizer.nneg_Cval)
    empty!(optimizer.nneg_Avar)
    empty!(optimizer.nneg_Acon)
    empty!(optimizer.nneg_Aval)
    empty!(optimizer.nneg_Bvar)
    empty!(optimizer.nneg_Bcon)
    empty!(optimizer.nneg_Bval)
    empty!(optimizer.nneg_L)
    empty!(optimizer.nneg_U)

    optimizer.objective_sense = MOI.FEASIBILITY_SENSE
    optimizer.objective_constant = 0.0

    optimizer.primal_objective_value = NaN
    optimizer.dual_objective_value = NaN
    empty!(optimizer.nneg_X)
    empty!(optimizer.psdc_X)
    empty!(optimizer.s)
    empty!(optimizer.y)
    empty!(optimizer.y2)
    empty!(optimizer.nneg_Z1)
    empty!(optimizer.nneg_Z2)
    empty!(optimizer.psdc_Z1)
    empty!(optimizer.psdc_Z2)
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

MOI.supports_add_constrained_variables(::Optimizer, ::Type{MOI.Reals}) = false

const SupportedSets = Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}
MOI.supports_add_constrained_variables(::Optimizer, ::Type{<:SupportedSets}) = true
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:Union{MOI.EqualTo{Float64}, BoundsSet}})
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
    push!(optimizer.single_variable_mask, 0x0)
    return MOI.VariableIndex(length(optimizer.variable_info))
end
function _add_psdc_variable(optimizer::Optimizer, index_in_cone)
    push!(optimizer.variable_info, VariableInfo(PSD, length(optimizer.psdc_dims), index_in_cone))
    push!(optimizer.single_variable_mask, 0x0)
    return MOI.VariableIndex(length(optimizer.variable_info))
end
function _add_constrained_variables(optimizer::Optimizer, set::MOI.Nonnegatives)
    push!(optimizer.nneg_info, MOI.dimension(set))
    for i in 2:MOI.dimension(set)
        push!(optimizer.nneg_info, i)
    end
    for i in 1:MOI.dimension(set)
        push!(optimizer.nneg_L, -Inf)
        push!(optimizer.nneg_U, Inf)
    end
    return [_add_nonneg_variable(optimizer) for i in 1:MOI.dimension(set)]
end
function _add_constrained_variables(optimizer::Optimizer, set::MOI.PositiveSemidefiniteConeTriangle)
    d = MOI.side_dimension(set)
    push!(optimizer.psdc_dims, d)
    push!(optimizer.psdc_Cvar, Int[])
    push!(optimizer.psdc_Cval, Float64[])
    push!(optimizer.psdc_Avar, Int[])
    push!(optimizer.psdc_Acon, Int[])
    push!(optimizer.psdc_Aval, Float64[])
    push!(optimizer.psdc_Bvar, Int[])
    push!(optimizer.psdc_Bcon, Int[])
    push!(optimizer.psdc_Bval, Float64[])
    push!(optimizer.psdc_L, fill(-Inf,  d, d))
    push!(optimizer.psdc_U, fill(Inf,  d, d))
    return [_add_psdc_variable(optimizer, i) for i in 1:MOI.dimension(set)]
end
function MOI.add_constrained_variables(optimizer::Optimizer, set::SupportedSets)
    vis = _add_constrained_variables(optimizer, set)
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(set)}(first(vis).value)
    return vis, ci
end

# Variable bounds
function MOI.supports_constraint(
    optimizer::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<:BoundsSet})
    return true
end
function MOI.add_constraint(
    optimizer::Optimizer, f::MOI.SingleVariable, set::BoundsSet)
    vi = f.variable
    flag = MOIU.single_variable_flag(typeof(set))
    mask = optimizer.single_variable_mask[vi.value]
    MOIU.throw_if_lower_bound_set(vi, typeof(set), mask, Float64)
    MOIU.throw_if_upper_bound_set(vi, typeof(set), mask, Float64)
    info = optimizer.variable_info[vi.value]
    if !iszero(flag & MOIU.LOWER_BOUND_MASK)
        if info.variable_type == NNEG
            L = optimizer.nneg_L
        else
            L = optimizer.psdc_L[info.cone_index]
        end
        L[info.index_in_cone] = MOIU.extract_lower_bound(set)
    end
    if !iszero(flag & MOIU.UPPER_BOUND_MASK)
        if info.variable_type == NNEG
            U = optimizer.nneg_U
        else
            U = optimizer.psdc_U[info.cone_index]
        end
        U[info.index_in_cone] = MOIU.extract_upper_bound(set)
    end
    optimizer.single_variable_mask[vi.value] = mask | flag
    return MOI.ConstraintIndex{MOI.SingleVariable, typeof(set)}(vi.value)
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
function MOI.add_constraint(optimizer::Optimizer, func::MOI.ScalarAffineFunction{Float64}, set::BoundsSet)
    if !iszero(MOI.constant(func))
        # We use the fact that the initial function constant was zero to
        # implement getters for `MOI.ConstraintPrimal`.
        throw(MOI.ScalarFunctionConstantNotZero{
             Float64, typeof(func), typeof(set)}(MOI.constant(func)))
    end
    flag = MOIU.single_variable_flag(typeof(set))
    push!(optimizer.l, iszero(flag & MOIU.LOWER_BOUND_MASK) ? -Inf : MOIU.extract_lower_bound(set))
    push!(optimizer.u, iszero(flag & MOIU.UPPER_BOUND_MASK) ? Inf : MOIU.extract_upper_bound(set))
    con = length(optimizer.l)
    for term in func.terms
        info = optimizer.variable_info[term.variable_index.value]
        if info.variable_type == NNEG
            push!(optimizer.nneg_Bvar, info.index_in_cone)
            push!(optimizer.nneg_Bcon, con)
            push!(optimizer.nneg_Bval, term.coefficient)
        else
            @assert info.variable_type == PSD
            push!(optimizer.psdc_Bvar[info.cone_index], info.index_in_cone)
            push!(optimizer.psdc_Bcon[info.cone_index], con)
            coef = is_diagonal_index(info.index_in_cone) ? term.coefficient : term.coefficient / √2
            push!(optimizer.psdc_Bval[info.cone_index], coef)
        end
    end
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(con)
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
function symvec(Q::Vector)
    if length(Q) > 1
        error("Expected square matrix from SDPNAL but got a vector of length $(length(Q))")
    end
    return Q
end
function _vec(x::SparseMatrixCSC)
    if size(x, 2) > 1
        error("Expected a vector from SDPNAL but got a sparse matrix of size $(size(x))")
    end
    return Vector(x[:, 1])
end
_vec(x::Float64) = [x]
_vec(x::Vector) = x

function MOI.optimize!(optimizer::Optimizer)
    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:printlevel] = 0
    end

    blk1 = Any[]
    blk2 = Any[]
    m = length(optimizer.b)
    p = length(optimizer.l)
    At = SparseArrays.SparseMatrixCSC{Float64,Int}[]
    Bt = SparseArrays.SparseMatrixCSC{Float64,Int}[]
    # FIXME I get a strange failure with sparse vectors, need to investigate
    C = Union{Matrix{Float64}, Vector{Float64}}[]
    L = Vector{Union{Matrix{Float64}, Vector{Float64}}}(undef, 0)
    U = Vector{Union{Matrix{Float64}, Vector{Float64}}}(undef, 0)
    if !isempty(optimizer.psdc_dims)
        for (i, dim) in enumerate(optimizer.psdc_dims)
            vec_dim = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(dim))
            push!(At, sparse(optimizer.psdc_Avar[i], optimizer.psdc_Acon[i], optimizer.psdc_Aval[i], vec_dim, m))
            push!(Bt, sparse(optimizer.psdc_Bvar[i], optimizer.psdc_Bcon[i], optimizer.psdc_Bval[i], vec_dim, p))
            c = Vector(sparsevec(optimizer.psdc_Cvar[i], optimizer.psdc_Cval[i], vec_dim))
            Ci = zeros(dim, dim)
            Li_vec = optimizer.psdc_L[i]
            Li = zeros(dim, dim)
            Ui_vec = optimizer.psdc_U[i]
            Ui = zeros(dim, dim)
            k = 0
            for col in 1:dim
                for row in 1:(col - 1)
                    k += 1
                    Ci[row, col] = c[k] / 2
                    Ci[col, row] = c[k] / 2
                    Li[row, col] = Li_vec[k]
                    Li[col, row] = Li_vec[k]
                    Ui[row, col] = Ui_vec[k]
                    Ui[col, row] = Ui_vec[k]
                end
                k += 1
                Ci[col, col] = c[k]
                Li[col, col] = Li_vec[k]
                Ui[col, col] = Ui_vec[k]
            end
            push!(C, Ci)
            push!(L, Li)
            push!(U, Ui)
            push!(blk1, "s")
            push!(blk2, Float64(dim))
        end
    end
    if !iszero(optimizer.num_nneg)
        push!(At, sparse(optimizer.nneg_Avar, optimizer.nneg_Acon, optimizer.nneg_Aval, optimizer.num_nneg, m))
        push!(Bt, sparse(optimizer.nneg_Bvar, optimizer.nneg_Bcon, optimizer.nneg_Bval, optimizer.num_nneg, p))
        push!(C, Vector(sparsevec(optimizer.nneg_Cvar, optimizer.nneg_Cval, optimizer.num_nneg)))
        push!(blk1, "l")
        push!(blk2, Float64(optimizer.num_nneg))
        push!(L, optimizer.nneg_L)
        push!(U, optimizer.nneg_U)
    end
    blk = [blk1 blk2]

    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:printlevel] = 0
    end

    obj, X, s, y, Z1, Z2, y2, v, optimizer.info, optimizer.runhist = sdpnalplus(
        blk, At, C, optimizer.b, L, U, Bt, optimizer.l, optimizer.u; options...)

    optimizer.primal_objective_value, optimizer.dual_objective_value = obj
    k = 0
    optimizer.psdc_X = symvec.(X[k .+ eachindex(optimizer.psdc_dims)])
    optimizer.psdc_Z1 = symvec.(Z1[k .+ eachindex(optimizer.psdc_dims)])
    optimizer.psdc_Z2 = symvec.(Z2[k .+ eachindex(optimizer.psdc_dims)])
    k += length(optimizer.psdc_dims)
    if iszero(optimizer.num_nneg)
        empty!(optimizer.nneg_X)
        empty!(optimizer.nneg_Z1)
        empty!(optimizer.nneg_Z2)
    else
        k += 1
        optimizer.nneg_X = X[k]
        optimizer.nneg_Z1 = Z1[k]
        optimizer.nneg_Z2 = _vec(Z2[k])
    end
    optimizer.s = s
    optimizer.y = y
    optimizer.y2 = _vec(y2) # Should we used `v` instead ? They should be approximately equal
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
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, <:BoundsSet})
    MOI.check_result_index_bounds(optimizer, attr)
    return MOI.get(optimizer, MOI.VariablePrimal(attr.N),
                   MOI.VariableIndex(ci.value))
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
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         <:BoundsSet})
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.s[ci.value]
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, <:SupportedSets})
    MOI.check_result_index_bounds(optimizer, attr)
    info = optimizer.variable_info[ci.value]
    if info.variable_type == NNEG
        dim = optimizer.nneg_info[info.index_in_cone]
        return optimizer.nneg_Z1[(info.index_in_cone - 1) .+ (1:dim)]
    else
        @assert info.variable_type == PSD
        return optimizer.psdc_Z1[info.cone_index]
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, <:BoundsSet})
    MOI.check_result_index_bounds(optimizer, attr)
    info = optimizer.variable_info[ci.value]
    if info.variable_type == NNEG
        Z2 = optimizer.nneg_Z2
    else
        @assert info.variable_type == PSD
        Z2 = optimizer.psdc_Z2[info.cone_index]
    end
    return Z2[info.index_in_cone]
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual, ci::AFFEQ)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.y[ci.value]
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         <:BoundsSet})
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.y2[ci.value]
end
