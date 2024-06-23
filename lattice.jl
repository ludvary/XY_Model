module LatticeJl
    using LinearAlgebra

    mutable struct Lattice
        L :: Int64
        spins :: Matrix{Vector{Float64}}
        E :: Float64
        M :: Vector{Float64}

        function Lattice(L::Int64)

            new(L, [normalize([1,1]-2*rand(2)) for _ in 1:L, _ in 1:L], 0.0, [0.0, 0.0])        # [normalize(...) for ...] creates a matrix of vectors [[x1,y1], [x2, y2], ....] which
                                                                                                # are all normalized
        end
    end

    function get_neighbours(i::Int64, j::Int64, lattice::Lattice)                               # returns a list containing all four nearest neighbours
        
        neighbours = [
            lattice.spins[1+mod(i-2, lattice.L), j],
            lattice.spins[1+mod(i, lattice.L), j],
            lattice.spins[i, 1+mod(j-2, lattice.L)],
            lattice.spins[i, 1+mod(j, lattice.L)]
        ] 

        return neighbours
    end
end
    
