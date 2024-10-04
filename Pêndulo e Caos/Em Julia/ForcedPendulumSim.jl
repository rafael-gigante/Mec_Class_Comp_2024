module ForcedPendulumSim
using Plots
using DifferentialEquations, Plots

###############################################
#                  Structs                    #
############################################### 

# Define the struct for the pendulum
mutable struct ForcedPendulum
    OMEGA::Float64   # Frequency of external force
    omega_r::Float64 # omega_0 from the slides (natural frequency of the system)
    beta::Float64    # Damping coefficient, dimensionless
    f_0::Float64     # Amplitude of the driving force, dimensionless
    theta_0::Float64 # Initial angle
    omega_0::Float64 # Initial angular velocity
end

###############################################
#               Main Functions                #
###############################################

# Differential equations for the forced pendulum system
function pendulum_eqs!(du, u, p, t)
    theta, omega = u
    OMEGA, omega_r, beta, f_0 = p
    du[1] = omega  # dθ/dt = ω
    du[2] = (omega_r^2 / OMEGA^2) * (-sin(theta) - 2 * beta * omega * (OMEGA / omega_r^2) + f_0 * sin(t))
end

# Function to run one simulation
function run_simulation(pendulum::ForcedPendulum, t_span::Tuple{Float64, Float64}, dt::Float64)
    # Setting up the problem
    u0 = [pendulum.theta_0, pendulum.omega_0]
    p = [pendulum.OMEGA, pendulum.omega_r, pendulum.beta, pendulum.f_0]
    prob = ODEProblem(pendulum_eqs!, u0, t_span, p)
    
    # Solving the problem
    sol = solve(prob, RK4(), dt=dt)

    # Plotting the results
    plot_fig = plot(sol.t, [u[1] for u in sol.u], label=false, title="Forced Pendulum Simulation", xlabel="Time (s)", ylabel="θ(t) (radians)", color=:black)
    savefig(plot_fig, "images/theta.png")

    plot_fig = plot([u[1] for u in sol.u], [u[2] for u in sol.u], label=false, title="Forced Pendulum Simulation \n Phase Diagram", xlabel="θ(t) (radians)", ylabel="ω(t) (radians/s)", color=:black)
    savefig(plot_fig, "images/phase_diagram.png")

    # Poincaré section  calculations#
    T = 2π  # Period of the driving force
    period_samples = Int(round(T/dt))  # Calculate how many samples per period based on dt

    # Collect every nth sample where n is the number of samples per period
    θ_values = [u[1] for (i, u) in enumerate(sol.u) if i % period_samples == 0]
    ω_values = [u[2] for (i, u) in enumerate(sol.u) if i % period_samples == 0]

    plot_fig = scatter([u[1] for u in sol.u], [u[2] for u in sol.u], label=false, title="Forced Pendulum Simulation \n Poincaré section", xlabel="θ(t) (radians)", ylabel="ω(t) (radians/s)", markersize=3, markercolor=:black)
    savefig(plot_fig, "images/poincare.png")
end

# Function to run the simulation and plot bifurcation diagram
function generate_bifurcation_diagram(pendulum::ForcedPendulum, F_D_range, t_final, dt)
    theta_values = []
    F_D_values = []
    t_span = (0.0, t_final) 

    for F_D in F_D_range
        # Update the driving force amplitude
        pendulum.f_0 = F_D

        # Setup the problem with the current F_D
        u0 = [pendulum.theta_0, pendulum.omega_0]
        p = [pendulum.OMEGA, pendulum.omega_r, pendulum.beta, pendulum.f_0]
        prob = ODEProblem(pendulum_eqs!, u0, t_span, p)

        # Solve the differential equation
        sol = solve(prob, RK4(), dt=dt)

        # Collect steady state or long-term behavior of theta after initial transients
        push!(theta_values, sol.u[end][1])  # Collect the final value of theta
        push!(F_D_values, F_D)
    end

    # Plot the bifurcation diagram
    plot_fig = scatter(F_D_values, theta_values, title="Forced Pendulum Simulation \n Bifurcation Diagram", xlabel="F_D", ylabel="θ_final (Radians)", label=false, markersize=1, markercolor=:black)
    savefig(plot_fig, "images/bifurcation.png")

end

###############################################
#               Usage Example                 #
###############################################

# Example initialization and simulation run
function example1()
    pendulum = ForcedPendulum(2.0, 1.0, 0.1, 0.5, π/4, 0.0)
    run_simulation(pendulum, (0.0, 1000.0), 0.01)
end

# Example function to set parameters and run the diagram generation
function example2()
    pendulum = ForcedPendulum(2.0, 1.0, 0.1, 0.5, π/4, 0.0)
    F_range = 1.00:0.01:2.00
    t_final = 100000.0           
    dt = 0.01                 
    generate_bifurcation_diagram(pendulum, F_range, t_final, dt)
end

end