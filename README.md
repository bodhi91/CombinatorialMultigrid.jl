
# CombinatorialMultigrid.jl
This package implements the Combinatorial Multigrid Preconditioner. Refer to [1] for details on the algorithm. 


In order to run CMG we present a quick example. Lets load a 1Mx1M example matrix ```X``` and build the ```b``` side. 

```
## load example matrix
file = matopen("../example/X.mat"); X = read(file, "X"); close(file)
LX = lap(X);
b = rand(Float64, size(X, 1));
b = b .- sum(b)/length(b);
```
CMG needs to be built before solving a linear system. To build CMG we provide two wrapper functions: ```cmg_preconditioner_lap```, which requires the user to provide the laplacian matrix and ```cmg_preconditioner_adj``` which requires the user to provide the adjacent matrix. 
We run the following script:

```
t = @elapsed (pfunc, h) = cmg_preconditioner_lap(LX);
@info "Time Required to build CMG Solver: $(t) seconds"
t = @elapsed x = pfunc(b);
@info "Time Required to find x: $(t) seconds"
```
Both ```cmg_preconditioner_lap``` and ```cmg_preconditioner_adj``` returns two parameters: the solver function and the hierarchy.
The above script generates the following output: 
```
[ Info: Time Required to build CMG Solver: 3.258120718 seconds
[ Info: Time Required to find x: 0.194163587 seconds
```
We try to solve a linear system using CMG. For this purpose we leverage ```pcg``` from the ```Laplacians``` package. We run the following script: 
```
f = pcgSolver(LX,pfunc);
t = @elapsed x = f(b, maxits=40, tol=1e-6,verbose=true);
@info "Time Required to solve system: $(t) seconds"
```
This generates the following output: 
```
PCG stopped after: 7.828 seconds and 29 iterations with relative error 8.429320186485909e-7.
[ Info: Time Required to solve system: 8.034237971 seconds
```

For comparison we run ```approxchol_lap``` which is the fastest solver from the ```Laplacians``` package. We run the following script: 
```
solver = approxchol_lap(X; tol=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams());
@info "Time Required to build Lap Solver: $(t) seconds"
t = @elapsed x = solver(b);
@info "Time Required to find x: $(t) seconds"
```
This generates the following output: 
```
[ Info: Time Required to build Lap Solver: 30.104803977 seconds
[ Info: Time Required to find x: 12.226966502 seconds
```

```CMG``` builds the solver in ```3.26 seconds``` compared to ```30 seconds``` with ```approxchol_lap``` and solves ```x``` in ```0.19 seconds``` compared to ```12.23 seconds```.


**Citations:***
*[1] Ioannis Koutis, Gary L. Miller, David Tolliver, Combinatorial preconditioners and multilevel solvers for problems in computer vision and image processing, Computer Vision and Image Understanding, Volume 115, Issue 12, 2011, Pages 1638-1646, ISSN 1077-3142, https://doi.org/10.1016/j.cviu.2011.05.013.*
