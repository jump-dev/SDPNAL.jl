# SDPNAL.jl

[SDPNAL.jl](https://github.com/jump-dev/SDPNAL.jl) is wrapper for the
[SDPNALplus](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus) solver.

The wrapper has two components:

 * an exported `sdpnalplus` function that is a thin wrapper on top of the
   `sdpnalplus` MATLAB function
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is maintained by the JuMP community and is not an official wrapper
of SDPNALplus.

## License

`SDPNAL.jl` is licensed under the [MIT License](https://github.com/jump-dev/SDPNAL.jl/blob/master/LICENSE.md).

The underlying solver, [SDPNALplus](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)
is licensed under the [Creative Commons Attribution-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-sa/4.0/).

In addition, SDPNAL requires an installation of MATLAB, which is a closed-source
commercial product for which you must [obtain a license](https://www.mathworks.com/products/matlab.html).

## Use with JuMP

To use SDPNAL with [JuMP](https://github.com/jump-dev/JuMP.jl), do:
```julia
using JuMP, SDPNAL
model = Model(SDPNAL.Optimizer)
set_attribute(model, "printlevel", 0)
```

Note that, contrary to implementation of other solver-independent interfaces,
using SDPNAL from JuMP or MOI fully exploits the particular structures of the
SDPNAL interface and does not create superfluous slack variables and equality
constraints as discussed in [the SDPNAL guide](https://arxiv.org/pdf/1710.10604.pdf):

> A new interface is necessary to facilitate the modeling of an SDP problem for
> SDPNAL+ because of latter’s flexibility to directly accept inequality
> constraints of the form “l ≤ B(X) ≤ u”, and bound constraints of the form
> “L ≤ X ≤ U”. The flexibility can significantly simplify the generation of the
> data in the SDPNAL+ format as compared to what need to be done in CVX or
> YALMIP to reformulate them as equality constraints through introducing extra
> variables. In addition, the final number of equality constraints present in
> the data input to SDPNAL+ can also be substantially fewer than those present
> in CVX or YALMIP. It is important to note here that the number of equality
> constraints present in the generated problem data can greatly affect the
> computational efficiency of the solvers, especially for interior-point based
> solvers.

## Installation

First, make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package, and that
the SDPNALplus software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.

Then, install `SDPNAL.jl` using `Pkg.add`:
```julia
import Pkg
Pkg.add("SDPNAL")
```

There is a `startup.m` file at the root of the SDPNAL folder. This adds all
subdirectories recursively when MATLAB starts. However, the `interface` directory
contains a `.git` subdirectory which contains a very large number of files.
Because of this, MATLAB crashes if SDPNAL is in its path because the `startup.m`
requests MATLAB to try to parse all the files in the `.git` folder. To resolve
this problem, delete the `startup.m` file and `.git` folder, and add the
subdirectories manually your `toolbox/local/pathdef.m` file as follows:
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

If you have [SDPT3](https://github.com/jump-dev/SDPT3.jl) in addition to SDPNAL
in the MATLAB path (that is, the `toolbox/local/pathdef.m` file) then you
might have issues because both solvers define a `validate` function, and this
might make SDPNAL call SDPT3's `validate` function instead of SDPT3's `validate`
function.
