using Test; Laplacians; using LinearAlgebra; using BenchmarkTools
include("../src/CombinatorialMultigrid.jl");

using CombinatorialMultigrid

## load microchip
X = wtedChimera(100_000)
LX = lap(X);
b1 = rand(Float64, size(X, 1));
b1 = b1 .- sum(b1)/length(b1);

@time (pfunc, h) = cmg_preconditioner_lap(LX);
x= pfunc(b1);
@time x= pfunc(b1);
