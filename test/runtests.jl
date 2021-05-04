using Test; Laplacians; using MAT; using LinearAlgebra; using BenchmarkTools
include("../src/CombinatorialMultigrid.jl");

## load microchip
file = matopen("X.mat"); X = read(file, "X"); close(file)
LX = lap(X);
b1 = rand(Float64, size(X, 1));
b1 = b1 .- sum(b1)/length(b1);

@time (pfunc, h) = solve_cmg(LX);
x= pfunc(b1);
@time x= pfunc(b1);
