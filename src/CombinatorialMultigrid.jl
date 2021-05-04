module CombinatorialMultigrid
    using SparseArrays
    using LinearAlgebra
    using LDLFactorizations

    include("cmgAlg.jl")
    export solve_cmg
end
