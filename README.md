# SDPNAL

`SDPNAL.jl` is an interface to the **[SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)**
solver. It exports the `sdpt3` function that is a thin wrapper on top of the
`sdpt3` MATLAB function and use it to define the `SDPNAL.Optimizer` object that
implements the solver-independent `MathOptInterface` API.

To use it with JuMP, simply do
```julia
using JuMP
using SDPNAL
model = Model(SDPNAL.Optimizer)
```
To suppress output, do `set_silent(model)` or
```julia
model = Model(optimizer_with_attribute(SDPNAL.Optimizer, "printlevel" => 0))
set_silent(model)
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
] add SDPNAL
```
but you first need to make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the SDPNAL software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.

### Troubleshooting

There is a `startup.m` file at the root of the SDPNAL folder.
This adds all subfolder recursively when MATLAB starts.
However, the `interface` subfolder constains a .git subfolder which contains a very large tree of subfolders.
Because of this MATLAB crashes if SDPNAL is in its path because the `startup.m` requests MATLAB to try to parse all the files in the `.git` folder.
To resolve this, delete the `startup.m` file and `.git` folder and add the subfolders manually your `toolbox/local/pathdef.m` file as follows:
```
function p = pathdef

% (...)

p = [...
%%% BEGIN ENTRIES %%%
'/path/to/SDPNALv1.0:', ...
'/path/to/SDPNALv1.0/interface:', ...
'/path/to/SDPNALv1.0/mexfun:', ...
'/path/to/SDPNALv1.0/solver:', ...
'/path/to/SDPNALv1.0/solver_main_default:', ...
'/path/to/SDPNALv1.0/util:', ...
% (...)
```

If you have [SDPT3](https://github.com/jump-dev/SDPT3.jl) in addition to SDPNAL in the MATLAB's path (i.e. the `toolbox/local/pathdef.m` file) then you might have issues.
As both solvers define a `validate` function, this might make SDPNAL call SDPT3's `validate` function instead of SDPT3's `validate` function.
