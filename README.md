
# CombinatorialMultigrid.jl
Implements the Combinatorial Multigrid Preconditioner


To run cmg do the following: 
```
using Laplacians
include("driver.jl")
import Main.cmg: solve_cmg

A = wtedChimera(10000)
LA = lap(A) 
b = rand(Float64, size(A, 1))

(pfunc, h) = solve_cmg(LA)
v = pfunc(b)
```
solve_cmg returns the functor and the hierarchy. In order to integrate CMG with an iterative solver like LOBPCG, wrap the functor with lPreconditioner and pass it as a preconditioner to LOBPCG. 
For example: 
```
lobpcg(LA, largest=false, nev=1, maxiter=55, P=lPreconditioner(pfunc))
```
