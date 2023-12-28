module LatticeJl
    using LinearAlgebra

    mutable struct Lattice
        L :: Int64
        spins :: Matrix{Vector{Float64}}

        function Lattice(L::Int64)
            new(L, [normalize([1,1]-2*rand(2)) for _ in 1:L, _ in 1:L])
        end
    end

    function get_neighbours(i::Int64, j::Int64, lattice::Lattice)
        
        neighbours = [
            lattice.spins[1+mod(i-2, lattice.L), j],
            lattice.spins[1+mod(i, lattice.L), j],
            lattice.spins[i, 1+mod(j-2, lattice.L)],
            lattice.spins[i, 1+mod(j, lattice.L)]
        ] 

        return neighbours
    end
end
    
