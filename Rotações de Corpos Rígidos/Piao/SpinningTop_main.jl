include("SpinningTopSim.jl")
using .SpinningTopSim
using Plots

# Parâmetros iniciais do sistema
I = [1.0, 1.0, 2.0]         # Momentum of inertia
M = 100.0                     # Mass
ℓ = 1.0                     # Center of mass distance from the point (0,0)
Ω = 10.0                     # Angular velocity of its physical axis
θ₀ = 0.4                    # θ inicial
dt = 0.01                 # Time step
simulation_time = 4.0      # Simulation duration

# Define o piao
spinningtop = SpinningTopSim.SpinningTop(I, M, ℓ, Ω, θ₀, θ₀- 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
output = SpinningTopSim.Output([], [], [], [])

# Simulation
SpinningTopSim.run_simulation!(spinningtop, output, simulation_time, dt)

# Analyze data
SpinningTopSim.generate_gif(simulation_time, dt, output.θ, output.φ, output.ψ, output.E)