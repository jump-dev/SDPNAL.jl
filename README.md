# SDPNAL

`SDPNAL.jl` is an interface to the **[SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)**
solver. It exports the `sdpt3` function that is a thin wrapper on top of the
`sdpt3` MATLAB function and use it to define the `SDPNAL.Optimizer` object that
implements the solver-independent `MathOptInterface` API.

To use it with JuMP, simply do
```julia
using JuMP
using SDPNAL
model = Model(with_optimizer(SDPNAL.Optimizer))
```
To suppress output, do
```julia
model = Model(with_optimizer(SDPNAL.Optimizer, printlevel=0))
```

Note that contrary to implementation of other solver-independent interfaces,
using SDPNAL from JuMP/MOI allows to fully exploit the particular structures of the SDPNAL interface
and does not create superfluous slack variables and equality constraints discussed in [the SDPNAL guide](https://arxiv.org/pdf/1710.10604.pdf):

> A new interface is necessary to facilitate the modeling of an SDP problem for Sdpnal+ because of latter’s flexibility to directly accept inequality constraints of the form “l ≤ B(X) ≤ u”,
> and bound constraints of the form “L ≤ X ≤ U”.
> The flexibility can significantly simplify the generation of the data in the Sdpnal+ format as compared
> to what need to be done in CVX or YALMIP to reformulate them as equality constraints through introducing extra variables.
> In addition, the final number of equality constraints present in the data input to Sdpnal+ can also be substantially fewer than those present in CVX or YALMIP.
> It is important to note here that the number of equality constraints present in the generated problem data can greatly affect the computational efficiency
> of the solvers, especially for interior-point based solvers.

## Installation

You can install `SDPNAL.jl` through the Julia package manager:
```julia
] add https://github.com/JuliaOpt/SDPNAL.jl.git
```
but you first need to make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the SDPNAL software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.
