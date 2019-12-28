module SDPNAL

using SparseArrays
using MATLAB

export sdpnalplus

# See solver_main_default/sdpnalplus.m
const ALLOWED_OPTIONS = [
    "tol",
    "sigma",
    "maxiter",
    "maxtime",
    "printlevel",
    "rescale",
    "stopoption",
    "AATsolve",
    "BBTsolve",
    "beta",
    "SSNprecond",
    "ADMmaxiter",
    "ADMtol"
]

_array(x::AbstractMatrix) = x
_array(x::Vector) = x
_array(x::Float64) = [x]

# TODO log in objective, OPTION, initial iterates X0, y0, Z0
# Solve the primal/dual pair
# min c'x,      max b'y
# s.t. Ax = b,   c - A'x ∈ K
#       x ∈ K
function sdpnalplus(
    blk::Matrix,
    At::Vector{<:Union{Matrix{Float64}, SparseMatrixCSC{Float64}}},
    C::Vector{<:Union{Matrix{Float64}, Vector{Float64}}},
    b::Vector{Float64},
    L = [], U = [],
    Bt = [], l = [], u = []; kws...)

    #C::Vector{<:Union{Vector{Float64}, SparseVector{Float64}}}, b::Vector{Float64})

    options = Dict{String, Any}(string(key) => value for (key, value) in kws)
    @assert all(i -> size(At[i], 2) == length(b), 1:length(At))
    @assert length(At) == size(blk, 1)
    #@assert all(i -> size(A[i], 1) == dim(blk[i, 1], blk[i, 2]), 1:length(A))
    #@assert all(i -> length(C[i], 1) == dim(blk[i, 1], blk[i, 2]), 1:length(A))
    # There are 10 output arguments so we use `10` below
    obj, X, s, y, Z1, Z2, y2, v, info, runhist = mxcall(:sdpnalplus, 10, blk, At, C, b, L, U, Bt, l, u, options)
    return obj, _array.(X), _array(s), _array(y), _array.(Z1), _array.(Z2), y2, v, info, runhist
end

include("MOI_wrapper.jl")

end # module
