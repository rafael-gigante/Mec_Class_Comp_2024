using DifferentialEquations
using Plots

# Function to define the equations of motion
function top_dynamics!(du, u, p, t)
    # Unpack variables
    ω1, ω2, ω3, φ, θ, ψ = u  # Angular velocities and Euler angles
    I1, I2, I3 = p           # Moments of inertia

    # Euler's equations for the angular velocities
    du[1] = (I2 - I3) / I1 * ω2 * ω3
    du[2] = (I3 - I1) / I2 * ω3 * ω1
    du[3] = (I1 - I2) / I3 * ω1 * ω2

    # Equations for the Euler angles (orientation of the top)
    du[4] = ω1 + ω2 * sin(φ) * tan(θ) + ω3 * cos(φ) * tan(θ)
    du[5] = ω2 * cos(φ) - ω3 * sin(φ)
    du[6] = ω2 * sin(φ) / cos(θ) + ω3 * cos(φ) / cos(θ)
end

# Initial conditions: [ω1, ω2, ω3, φ, θ, ψ]
u0 = [5.0, 0.0, 10.0, 0.0, 0.1, 0.0]

# Parameters: Moments of inertia I1, I2, I3
p = [2.0, 2.0, 1.0]

# Time span for the simulation
tspan = (0.0, 10.0)

# Define the problem and solve it
prob = ODEProblem(top_dynamics!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Plot the angular velocities over time
plot(sol.t, sol[1, :], label="ω₁", xlabel="Time", ylabel="Angular Velocities", title="Spinning Top Dynamics")
plot!(sol.t, sol[2, :], label="ω₂")
plot!(sol.t, sol[3, :], label="ω₃")
