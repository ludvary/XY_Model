module LatticeJl
using LinearAlgebra

mutable struct Lattice
    L::Int64
    coords::Vector{Vector{Float64}}
    spins::Vector{Vector{Float64}} # this contains the angle in radians

    function Lattice(L)
        new(L, [[i, j] for i in 1:L for j in 1:L], [normalize(randn(2)) for _ in 1:L^2]) # normalize(randn(2)) makes randomnly oriented unit Vectors
    end
end


# This can be made faster

function get_neighbour_indices(lattice_point::Int64, lattice::Lattice)

    neighbour_indices = [lattice_point + 1, lattice_point - 1, lattice_point + lattice.L, lattice_point - lattice.L]

    if lattice_point <= lattice.L  # means it is the top row
        if (lattice_point) == lattice.L  # means it is top right
            neighbour_indices = [1, lattice.L, lattice_point + lattice.L, lattice_point - 1]

        elseif (lattice_point) == 1  # means it is top left
            neighbour_indices = [2, lattice_point + lattice.L, lattice.L, (lattice.L)^2 - lattice.L + 1]
        else
            neighbour_indices = [lattice_point + 1, lattice_point - 1, lattice_point + lattice.L, (lattice.L)^2 - (lattice.L - lattice_point)]  # means it is in top row but not on the extremes
        end
    end

    if lattice_point > (lattice.L)^2 - lattice.L  # means it is bottom row
        if lattice_point == (lattice.L)^2  # means it is bottom right
            neighbour_indices = [lattice.L, lattice_point - lattice.L + 1, lattice.L, (lattice_point) - 1]

        elseif mod(lattice_point, lattice.L) == 1  # means it is bottom left
            neighbour_indices = [lattice_point + 1, lattice_point - lattice.L, (lattice.L)^2, 1]

        else
            neighbour_indices = [lattice_point + 1, lattice_point - 1, lattice_point - lattice.L, lattice.L - ((lattice.L)^2 - lattice_point)] # means it is in the bottom rown not in the extremes
        end
    end

    if mod(lattice_point, lattice.L) == 1 && lattice_point != 1 && lattice_point != (lattice.L)^2 - (lattice.L - 1)  # means at the leftmost column but not at top or bottom
        neighbour_indices = [lattice_point + 1, lattice_point + lattice.L, lattice_point - lattice.L, lattice_point + (lattice.L - 1)]
    end

    if mod(lattice_point, lattice.L) == 0 && lattice_point != (lattice.L)^2 && lattice_point != lattice.L # means at the rightmost column but not at top or bottom
        neighbour_indices = [lattice_point - 1, lattice_point + lattice.L, lattice_point - lattice.L, lattice_point - (lattice.L - 1)]
    end

    return neighbour_indices
end

#=
for every point in the lattice make a lookup_table as soon as the code starts and when you need to find neighbour_indices (which you need to do like sim_steps*eqm_steps of times, just use the lookup table)
=#

function make_lookup_table(lattice::Lattice) 
    N = lattice.L^2

    lookup_table = Dict{Int64, Vector{Int64}}()

    for i in 1:N
        neighbour_indices = get_neighbour_indices(i, lattice)
        lookup_table[i] = neighbour_indices
    end

    return lookup_table
end

end

