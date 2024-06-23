include("even_faster_lattice.jl")
using .LatticeJl
using Statistics
using LinearAlgebra
using PyPlot
using ProgressMeter
using PolyFit


const L = 10                        # side length of square lattice
const eqm_steps = 100_000           # number of iterations (not mc steps) in one call of do_eqm() function
const sim_steps = 20_000            # number of iterations (not mc steps) in one call of do_sim() function
const config_avs = 2000
const N = L^2                       # total number of spins
const J = 1                         # greater this number, greater the interaction strength
const B_mag = 0                     # magnetic field
const B = [B_mag, 0]                # magnetic field
const δ = 0.1 * pi                  # the spins get perturbed by this angle 
const T_init = 0.01                 # intial temparature
const T_final = 4.0                 # final temparature
const T_step = 0.03                 # temparature step


function main()


    # function to find the magnetization of the whole lattice
    function find_mag(lattice :: LatticeJl.Lattice)
        return sum(lattice.spins)
    end


    # function to find the energy of the whole lattice
    function find_energy(lattice :: LatticeJl.Lattice)
        E = 0
        for i in 1:L
            for j in 1:L
                neighbours = LatticeJl.get_neighbours(i, j, lattice)
                E -= J * dot(lattice.spins[i, j], sum(neighbours)) + dot(B, lattice.spins[i, j])
            end
        end
        return E/2      # due to double counting
    end


    # rotation matrix in 2D
    function rotation_matrix(theta)

        R = [cos(theta) -sin(theta);
            sin(theta) cos(theta)]

        return R
    end


    
    # function to carry out the equillibrium steps
    function do_eqm(lattice :: LatticeJl.Lattice, Temperature :: Float64)

        for _ in 1:eqm_steps
            
            # randomly choose a lattice point
            i, j = rand(1:L), rand(1:L)

            # find the pertubation angle
            Δδ = 2 * δ * rand() - δ

            # get the sum of spins of nearest neighbours
            sum_neighbours = sum(LatticeJl.get_neighbours(i, j, lattice))

            # find the rotation matrix
            rotation_mat = rotation_matrix(Δδ)

            # find the change in local energy if the pertubation was to be accepted
            ΔE = J * dot(sum_neighbours, (I(2) - rotation_mat) * lattice.spins[i, j])

            # if change is -ve, accept the pertubation
            if ΔE < 0
                # update the values E, M, and spin 
                lattice.E += ΔE
                lattice.M += (rotation_mat - I(2)) * lattice.spins[i, j]
                lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]

            # else accept with some probability (kb is 1)
            else
                prob = exp(-ΔE / Temperature)

                if rand() < prob
                    # update the values E, M, and spin 
                    lattice.E += ΔE
                    lattice.M += (rotation_mat - I(2)) * lattice.spins[i, j]
                    lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]
                end
            end

        end
        # return the equillibriated(just invented a new word?) lattice 
        return lattice
    end

    # same as do_eqm() but this time we actually collect data
    function do_sim(lattice :: LatticeJl.Lattice, Temperature :: Float64)

        # init arrays for collecting data
        E_arr = zeros(Float64, eqm_steps)
        
        M_arr = zeros(Float64, eqm_steps)

        for k in 1:eqm_steps
            
            i, j = rand(1:L), rand(1:L)

            Δδ = 2 * δ * rand() - δ

            sum_neighbours = sum(LatticeJl.get_neighbours(i, j, lattice))
            rotation_mat = rotation_matrix(Δδ)

            ΔE = J * dot(sum_neighbours, (I(2) - rotation_mat) * lattice.spins[i, j])

            if ΔE < 0
                lattice.E += ΔE
                lattice.M += (rotation_mat - I(2)) * lattice.spins[i, j]
                lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]

            else
                prob = exp(-ΔE / Temperature)

                if rand() < prob
                    lattice.E += ΔE
                    lattice.M += (rotation_mat - I(2)) * lattice.spins[i, j]
                    lattice.spins[i, j] = rotation_mat * lattice.spins[i, j]
                end
            end

            # write E per spin and M per spin into the array (direction doesnt matter, we just want to see how aligned the spins are, hence the norm())
            E_arr[k] = lattice.E/N
            M_arr[k] = norm(lattice.M)/N
        end
        
        # return the lattice and the arrays
        return lattice, E_arr, M_arr
    end

    # number of times the temparature is updated
    number_of_temps = Int(floor(abs(T_final - T_init)/T_step))

    # init arrays
    T_arr = zeros(Float64, number_of_temps)
    E_avg_arr = zeros(Float64, number_of_temps)
    M_avg_arr = zeros(Float64, number_of_temps)
    Cv_arr = zeros(Float64, number_of_temps)
    Chi_arr = zeros(Float64, number_of_temps)

    # init matrices
    E_avg_mat = [zeros(Float64, number_of_temps) for _ in 1:config_avs]
    M_avg_mat = [zeros(Float64, number_of_temps) for _ in 1:config_avs]
    Cv_mat = [zeros(Float64, number_of_temps) for _ in 1:config_avs]
    Chi_mat = [zeros(Float64, number_of_temps) for _ in 1:config_avs]

    
    # driver code for update the temparature and calling do_eqm() and do_sim() functions
    @showprogress dt=1 desc="Computing...." for v in 1:config_avs

        println("at system $v of total $config_avs systems")

        lattice = LatticeJl.Lattice(L)  # Lattice() method of LatticeJl instantiates a Lattice object with four attributes, L, spins, E and M (see ./lattice.jl)
        lattice.E = find_energy(lattice)
        lattice.M = find_mag(lattice)

        for i in 1:number_of_temps
            T = T_init + i*T_step
            lattice = do_eqm(lattice, T)
            lattice, E_arr, M_arr = do_sim(lattice, T)

            E_avg_mat[v][i] = mean(E_arr)
            M_avg_mat[v][i] = mean(M_arr)
            Cv_mat[v][i] = N*var(E_arr)/T^2
            Chi_mat[v][i] = N*var(M_arr)/T
        end
    end

    println("Starting averaging across various configs")
    for i in 1:number_of_temps
        T = T_init + i*T_step
        T_arr[i] = T
        E_avg_arr[i] = sum(E_avg_mat[v][i] for v in 1:config_avs)/config_avs
        M_avg_arr[i] = sum(M_avg_mat[v][i] for v in 1:config_avs)/config_avs
        Cv_arr[i] = sum(Cv_mat[v][i] for v in 1:config_avs)/config_avs
        Chi_arr[i] = sum(Chi_mat[v][i] for v in 1:config_avs)/config_avs
    end

    
    # # curve fitting using a polynomial of degree 15
    # order = 15
    # poly_Cv = polyfit(T_arr, Cv_arr, order)         # this gives you coefficients of various powers of x from x^0 upto x^(order)
    # poly_Chi = polyfit(T_arr, Chi_arr, order)
    # poly_E = polyfit(T_arr, E_avg_arr, order)
    # poly_M = polyfit(T_arr, M_avg_arr, order)

    # # function to extract the coefficients and compute the interpolating polynomial at point x
    # function fitting(poly, x)
    #     value = 0
    #     for i in 0:order
    #         value += poly[i]*x^i
    #     end
    #     return value
    # end


    # plots
    # plot(T_arr, [fitting(poly_E, x) for x in T_arr], color="navy", label="Interpolation($(order)th order Poly)")
    plot(T_arr, E_avg_arr, lw=0.4, color="navy", alpha=0.6, label="<E>")
    scatter(T_arr, E_avg_arr, color="navy", s=2, alpha=0.5)
    xlabel("T", fontsize=20)
    ylabel("<E> per spin", fontsize=20)
    legend()
    show()


    # plot(T_arr, [fitting(poly_M, x) for x in T_arr], color="maroon", label="Interpolation ($(order)th order Poly)")
    plot(T_arr, M_avg_arr, lw=0.4, color="maroon", alpha=0.6, label="<M>")
    scatter(T_arr, M_avg_arr, color="maroon", s=2, alpha=0.5)
    xlabel("T", fontsize=20)
    ylabel("<M> per spin", fontsize=20)
    legend()
    show()

    # plot(T_arr, [fitting(poly_Cv, x) for x in T_arr], color="grey", label="Interpolation ($(order)th order Poly)")
    plot(T_arr, Cv_arr, lw=0.4, color="black", alpha=0.6, label="Cv")
    scatter(T_arr, Cv_arr, color="black", s=2, alpha=0.5)
    xlabel("T", fontsize=20)
    ylabel(L"$C_v$", fontsize=20)
    legend()
    show()

    # plot(T_arr, [fitting(poly_Chi, x) for x in T_arr], color="green", label="Interpolation ($(order)th order Poly)")
    plot(T_arr, Chi_arr, lw=0.4, color="brown", alpha=0.6, label=L"\chi")
    scatter(T_arr, Chi_arr, color="brown", s=2, alpha=0.5)
    xlabel("T", fontsize=20)
    ylabel(L"\chi", fontsize=20)
    legend()
    show()
end
main()
