include("lattice.jl")
using .LatticeJl
using Statistics
using LinearAlgebra
using PyPlot
using ProgressMeter
using PolyFit


const L = 300
const eqm_steps = 10_000_000
const sim_steps = 1_000_000
const N = L^2
const J = 1
const B_mag = 0
const B = [B_mag, 0]
const δ = 0.1 * pi
const T_init = 0.1
const T_final = 6.5
const T_step = 0.01

function main()
    lattice = LatticeJl.Lattice(L)

    function find_mag(lattice :: LatticeJl.Lattice)
        return sum(lattice.spins)
    end

    function find_energy(lattice :: LatticeJl.Lattice)
        E = 0
        for i in 1:L
            for j in 1:L
                neighbours = LatticeJl.get_neighbours(i, j, lattice)
                E -= J * dot(lattice.spins[i, j], sum(neighbours)) + dot(B, lattice.spins[i, j])
            end
        end
        return E/2
    end

    function rotation_matrix(theta)

        R = [cos(theta) -sin(theta);
            sin(theta) cos(theta)]

        return R
    end


    function do_eqm(lattice :: LatticeJl.Lattice, Temperature :: Float64)

        E = find_energy(lattice)
        M = find_mag(lattice)

        for _ in 1:eqm_steps
            
            i, j = rand(1:L), rand(1:L)
            Δδ = 2 * δ * rand() - δ

            sum_neighbours = sum(LatticeJl.get_neighbours(i, j, lattice))
            rotation_mat = rotation_matrix(Δδ)

            # ΔE = J * dot(sum_neighbours, ( (lattice.spins[i, j])' * (I(2) - rotation_mat)) )  # doesnt work
            ΔE = J * dot(sum_neighbours, (I(2) - rotation_mat) * lattice.spins[i, j])           # works

            if ΔE < 0
                E += ΔE
                M += (rotation_mat - I(2)) * lattice.spins[i, j]
                lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]

            else
                prob = exp(-ΔE / Temperature)

                if rand() < prob
                    E += ΔE
                    M += (rotation_mat - I(2)) * lattice.spins[i, j]
                    lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]
                end
            end

        end
        return lattice
    end

    function do_sim(lattice :: LatticeJl.Lattice, Temperature :: Float64)

        E = find_energy(lattice)
        M = find_mag(lattice)

        E_arr = zeros(Float64, eqm_steps)
        M_arr = zeros(Float64, eqm_steps)

        for k in 1:eqm_steps
            
            i, j = rand(1:L), rand(1:L)

            Δδ = 2 * δ * rand() - δ

            sum_neighbours = sum(LatticeJl.get_neighbours(i, j, lattice))
            rotation_mat = rotation_matrix(Δδ)

            # ΔE = J * dot(sum_neighbours, ( (lattice.spins[i, j])' * (I(2) - rotation_mat)) )  # doesnt work
            ΔE = J * dot(sum_neighbours, (I(2) - rotation_mat) * lattice.spins[i, j])           # works

            if ΔE < 0
                E += ΔE
                M += (rotation_mat - I(2)) * lattice.spins[i, j]
                lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]

            else
                prob = exp(-ΔE / Temperature)

                if rand() < prob
                    E += ΔE
                    M += (rotation_mat - I(2)) * lattice.spins[i, j]
                    lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]
                end
            end

            E_arr[k] = E/N
            M_arr[k] = norm(M)/N
        end
        return lattice, E_arr, M_arr
    end

    number_of_temps = Int(abs(T_final - T_init)/T_step)

    T_arr = zeros(Float64, number_of_temps)
    E_avg_arr = zeros(Float64, number_of_temps)
    M_avg_arr = zeros(Float64, number_of_temps)
    var_E_at_T_arr = zeros(Float64, number_of_temps)
    var_M_at_T_arr = zeros(Float64, number_of_temps)
    
    @showprogress dt=1 desc="Computing...." for i in 1:number_of_temps
        T = T_init + i*T_step
        lattice = do_eqm(lattice, T)
        lattice, E_arr, M_arr = do_sim(lattice, T)

        T_arr[i] = T
        E_avg_arr[i] = mean(E_arr)
        M_avg_arr[i] = mean(M_arr)
        var_E_at_T_arr[i] = var(E_arr)/T^2
        var_M_at_T_arr[i] = var(M_arr)/T
    end

    
    order = 15
    poly_Cv = polyfit(T_arr, var_E_at_T_arr, order)
    poly_Chi = polyfit(T_arr, var_M_at_T_arr, order)
    poly_E = polyfit(T_arr, E_avg_arr, order)
    poly_M = polyfit(T_arr, M_avg_arr, order)

    function fitting(poly, x)
        value = 0
        for i in 0:order
            value += poly[i]*x^i
        end
        return value
    end

    # plot(T_arr, E_avg_arr, label="<E>", lw=0.8, color="navy")
    # plot(T_arr, M_avg_arr, label="<M>", lw=0.8, color="maroon")
    # legend()
    # grid(true)
    # show()


    plot(T_arr, [fitting(poly_E, x) for x in T_arr], color="navy", label="Interpolation(10th order Poly)")
    plot(T_arr, E_avg_arr, lw=0.4, color="navy", alpha=0.6, label="<E>")
    xlabel("T", fontsize=20)
    ylabel("<E> per spin", fontsize=20)
    legend()
    show()


    plot(T_arr, [fitting(poly_M, x) for x in T_arr], color="maroon", label="Interpolation (10th order Poly)")
    plot(T_arr, M_avg_arr, lw=0.4, color="maroon", alpha=0.6, label="<M>")
    xlabel("T", fontsize=20)
    ylabel("<M> per spin", fontsize=20)
    legend()
    show()

    plot(T_arr, [fitting(poly_Cv, x) for x in T_arr], color="grey", label="Interpolation (10th order Poly)")
    plot(T_arr, var_E_at_T_arr, lw=0.4, color="grey", alpha=0.6, label="Cv")
    xlabel("T", fontsize=20)
    ylabel(L"$C_v$", fontsize=20)
    legend()
    show()

    plot(T_arr, [fitting(poly_Chi, x) for x in T_arr], color="green", label="Interpolation (10th order Poly)")
    plot(T_arr, var_M_at_T_arr, lw=0.4, color="green", alpha=0.6, label=L"\chi")
    xlabel("T", fontsize=20)
    ylabel(L"\chi", fontsize=20)
    legend()
    show()
end
main()
