
# CombinatorialMultigrid.jl
Implements the Combinatorial Multigrid Preconditioner


In order to run CMG we present a quick example. Lets load an example matrix ```X``` and build the ```b``` side. 

```
## load example matrix

X = wtedChimera(100_000);
LX = lap(X);
b = rand(Float64, size(X, 1));
b = b .- sum(b)/length(b);
```
CMG needs to be built before solving a linear system. To build CMG we provide two wrapper functions: ```cmg_preconditioner_lap```, which requires the user to provide the laplacian matrix and ```cmg_preconditioner_adj``` which requires the user to provide the adjacent matrix. 

```
## build cmg preconditioner 
t = @elapsed (pfunc, h) = cmg_preconditioner_lap(LX);
@info "Time Required to build CMG Solver: $(t) seconds"
t = @elapsed x = pfunc(b);
@info "Time Required to find x: $(t) seconds"

## solve with pcg and cmg preconditioner
f1 = pcgSolver(LX,pfunc);
t = @elapsed x = f1(b, maxits=40, tol=1e-6,verbose=true);
@info "Time Required to solve system: $(t) seconds"

## solve with approxchol_lap from laplacians
solver = approxchol_lap(X; tol=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams());
@info "Time Required to build Lap Solver: $(t) seconds"
t = @elapsed x = solver(b);
@info "Time Required to find x: $(t) seconds"
```
