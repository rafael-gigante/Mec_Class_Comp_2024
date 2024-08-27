using DifferentialEquations, Plots

# Circuit parameters
L = 0.5   # Inductance (H)
R = 1.0   # Resistance (Ohms)
C = 0.1   # Capacitance (F)
V0 = 5.0  # Voltage amplitude (V)
ω = 50.0  # Frequency (rad/s)

# Define the ODE system
function circuit!(du, u, p, t)
    V = V0 * sin(ω * t)
    I, VC = u

    du[1] = (V - R * I - VC) / L           # dI/dt
    du[2] = I / C                          # dVC/dt (Voltage across the capacitor)
end

# Initial conditions
u0 = [0.0, 0.0]  # Initial current and capacitor voltage

# Time span for the simulation
tspan = (0.0, 2)

# Solve the system of ODEs
prob = ODEProblem(circuit!, u0, tspan)
sol = solve(prob, Tsit5())

# Plotting the results
plot(sol, xlabel="Time (s)", ylabel="Current (A) and Voltage (V)", layout=(2, 1))
plot!(sol[1], label="Current (I)", subplot=1)
plot!(sol[2], label="Capacitor Voltage (V_C)", subplot=2)
