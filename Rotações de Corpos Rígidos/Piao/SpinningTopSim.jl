module SpinningTopSim
using Plots

###############################################
#                  Structs                    #
############################################### 

# Struct of a SpinningTop
mutable struct SpinningTop
    I::Vector{Float64} # Momentum of inertia
    M::Float64         # Mass
    ℓ::Float64         # Center of mass distance from the point (0,0)
    Ω::Float64         # Angular velocity of its physical axis
    θ₀::Float64        # θ inicial

    θ::Float64         # Nutation angle 
    φ::Float64         # Precession angle
    ψ::Float64         # Rotation of the spinning top physical axis

    θ_dot::Float64     # Velocity of the nutation angle 
    φ_dot::Float64     # Velocity of the precession angle
    ψ_dot::Float64     # Velocity of the rotation of the spinning top physical axis

    E::Float64         # Energy of the system

    θ_dir::Int64         # Direction of θ
end

# Struct of the arrays we will use as output
mutable struct Output
    θ::Vector{Float64}  # Nutation angle 
    φ::Vector{Float64}  # Precession angle
    ψ::Vector{Float64}  # Rotation of the spinning top physical axis
    E::Vector{Float64}  # Energy of the system
end

###############################################
#             Auxiliar Functions              #
############################################### 

# Function f(θ)
function f(spinningtop::SpinningTop)
    a = (2 * spinningtop.M * 9.81 * spinningtop.ℓ) / spinningtop.I[1]
    b = (spinningtop.I[3] * spinningtop.Ω / spinningtop.I[1])^2
    f1 = a * (cos(spinningtop.θ₀) - cos(spinningtop.θ))
    f2 = - b * ((cos(spinningtop.θ₀) - cos(spinningtop.θ))^2 / sin(spinningtop.θ)^2)
    f = sqrt(abs(f1 + f2))
    return f * spinningtop.θ_dir
end

# Function g(θ)
function g(spinningtop::SpinningTop)
    a = (spinningtop.I[3] * spinningtop.Ω / spinningtop.I[1])
    g = a * ((cos(spinningtop.θ₀) - cos(spinningtop.θ)) / sin(spinningtop.θ)^2)
    return g
end

###############################################
#               Main Functions                #
###############################################

# Energy calculation
function calculate_energy(θ::Float64, θ_dot::Float64, φ_dot::Float64, ψ_dot::Float64, I1::Float64, I3::Float64, M::Float64, ℓ::Float64)

    E1 = (I3 / 2) * (ψ_dot + φ_dot * cos(θ))^2
    E2 = (I1 / 2) * (θ_dot^2 + φ_dot^2 * sin(θ)^2)
    E3 = M * 9.81 * ℓ * cos(θ)

    E = E1 + E2 + E3
    return E
end

# What happens in one iteration
function atualizar_grandezas!(spinningtop::SpinningTop, dt::Float64)

    alpha = 0.01

    if abs(spinningtop.θ - spinningtop.θ₀) < alpha
        spinningtop.θ += alpha
        spinningtop.θ_dir *= -1
    end
    # Velocities
    spinningtop.θ_dot = f(spinningtop)
    spinningtop.φ_dot = g(spinningtop)
    spinningtop.ψ_dot = spinningtop.Ω - g(spinningtop) * cos(spinningtop.θ)

    # Angles
    spinningtop.θ = spinningtop.θ + spinningtop.θ_dot * dt
    spinningtop.φ = spinningtop.φ + spinningtop.φ_dot * dt
    spinningtop.ψ = spinningtop.ψ + spinningtop.ψ_dot * dt

    # Energy
    spinningtop.E = calculate_energy(spinningtop.θ, spinningtop.θ_dot, spinningtop.φ_dot, spinningtop.ψ_dot, 
                                     spinningtop.I[1], spinningtop.I[3], spinningtop.M, spinningtop.ℓ)
end

function run_simulation!(spinningtop::SpinningTop, output::Output, simulation_time::Float64, dt::Float64)
    
    # Iterations
    iterations = Int(simulation_time / dt)

    # Initializing energy
    spinningtop.θ_dot = f(spinningtop)
    spinningtop.φ_dot = g(spinningtop)
    spinningtop.ψ_dot = spinningtop.Ω - g(spinningtop) * cos(spinningtop.θ)
    spinningtop.E = calculate_energy(spinningtop.θ, spinningtop.θ_dot, spinningtop.φ_dot, spinningtop.ψ_dot, 
                     spinningtop.I[1], spinningtop.I[3], spinningtop.M, spinningtop.ℓ)

    # Initializing output
    push!(output.θ, spinningtop.θ)
    push!(output.φ, spinningtop.φ)
    push!(output.ψ, spinningtop.ψ)
    push!(output.E, spinningtop.E)

    # Main loop of iterations
    for _ in 1:iterations
        # Generate new values
        atualizar_grandezas!(spinningtop, dt)

        # Append new value to the list
        push!(output.θ, spinningtop.θ)
        push!(output.φ, spinningtop.φ)
        push!(output.ψ, spinningtop.ψ)
        push!(output.E, spinningtop.E)
    end
end

###############################################
#           Visualization Functions           #
############################################### 

# Function to convert spherical coordinates to Cartesian coordinates
function sphericalToCartesian(r, theta, phi)
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z
end

# Plot the piao in 3D and generate an animation for every t in time
function animatePiao(SimulationTime, TimeStep, thetaValues, phiValues)
    r = 1.0  # Assume the piao has length 1.0
    anim = Animation()  # Create an Animation object

    # Set the desired animation length and frame rate
    duration = SimulationTime  # Duration in seconds
    max_fps = 60  # Maximum frames per second
    total_frames = Int(duration * max_fps)  # Total number of frames
    t_values = range(0, stop=duration, length=total_frames)  # Time steps

    # Lists for the trajectory
    trajectory_x = Float64[]
    trajectory_y = Float64[]
    trajectory_z = Float64[]

    # Create the animation
    for t in t_values

        # Extract the angles at the current time
        t_index = Int(round(t / TimeStep))  # Find the closest index for time t
        t_index = min(max(t_index, 1), Int(SimulationTime / TimeStep))  # Ensure t_index stays within bounds
        theta = thetaValues[t_index]
        phi = phiValues[t_index]

        # Convert the spherical coordinates to Cartesian
        x, y, z = sphericalToCartesian(r, theta, phi)

        # Coordinates of the trajectory
        push!(trajectory_x, x)
        push!(trajectory_y, y)
        push!(trajectory_z, z)

        # Plot the piao in 3D
        plot3d([0, x], [0, y], [0, z], label=false, legend=false, 
                title="Evolução do Pião \nt = $(round(t, digits=1)) s",
                xlims=(-1, 1), ylims=(-1, 1), zlims=(0, 1), lw=2) # Piao
        scatter!([x], [y], [z], color=:red, ms=4) # Red dot
        plot3d!(trajectory_x, trajectory_y, trajectory_z, color=:red, lw=1) # Trajectory of the top

        frame(anim)  # Capture the frame
    end
    
    return anim
end

# Plot the piao in 3D and generate an animation for every t in time and a graph of all the angles evolution
function generate_gif(SimulationTime, TimeStep, thetaValues, phiValues, psiValues, EnergyValues)
    # Create time array for plotting
    timeValues = collect(0:TimeStep:SimulationTime)

    # Plot the angles evolution over time
    display(plot(timeValues, [thetaValues phiValues psiValues], 
        label=["Theta" "Phi" "Psi"], xlabel="Time (s)", 
        ylabel="Angle (rad)", title="Evolution of Angles over Time", linewidth=2))
    savefig("Piao/piao_angles.png")

    # Plot the total energy over time
    normalizedEnergyValues = EnergyValues / EnergyValues[1]  # Normalize using the first energy value

    display(plot(timeValues, normalizedEnergyValues, 
        label=false, xlabel="Time (s)", 
        ylabel="E / E_0", title="Energy of the system over time", linewidth=2))
    savefig("Piao/piao_energy.png")

    # Create an animation of the piao
    anim = animatePiao(SimulationTime, TimeStep, thetaValues, phiValues)
    gif(anim, "Piao/piao_evolution.gif")  # Save the animation as a GIF
end

end