using Plots

# Constants
g = [0.0, -9.81]              # Gravity (m/s^2)
L = 1.0                     # Length of the baton (meters)
ω = 10.0                     # Angular velocity (rad/s)
dt = 0.01                   # Time step (seconds)
t_total = 2.0               # Total simulation time (seconds)

# Initial conditions
r_cm = [0.0, 0.0]           # Initial position of the center of mass
v_cm = [7.5, 10.0]           # Initial velocity of the center of mass
θ = 0.0                  # Initial angle of the baton (radians)

# Arrays to store positions for animation
r1_list = []
r2_list = []

x_cm = []
y_cm = []
# Simulation loop
for t in 0:dt:t_total
    # Update center of mass velocity and position
    v_cm .= v_cm .+ g .* dt
    r_cm .= r_cm .+ v_cm .* dt

    push!(x_cm, r_cm[1])
    push!(y_cm, r_cm[2])

    # Update orientation
    θ += ω * dt

    # Calculate positions of the masses
    r1 = r_cm .+ (L / 2) .* [cos(θ), sin(θ)]
    r2 = r_cm .- (L / 2) .* [cos(θ), sin(θ)]

    # Store positions for animation
    push!(r1_list, r1)
    push!(r2_list, r2)
end

# Create the animation
animation = @animate for i in eachindex(r1_list)
    # Extract positions
    r1 = r1_list[i]
    r2 = r2_list[i]

    # Plot the baton
    plot([r1[1], r2[1]], [r1[2], r2[2]], lw=3, label="Baton", xlims=(-1, 20), ylims=(-1, 10))
    scatter!([r1[1], r2[1]], [r1[2], r2[2]], label="Masses", color=:red)
    plot!([x_cm[1:i]], [y_cm[1:i]], color=:black, linestyle =:dot, label=false)
    xlabel!("X Position (m)")
    ylabel!("Y Position (m)")
    title!("Throwing Baton Simulation")
end

# Save or display the animation
gif(animation, "throwing_baton_simulation.gif", fps=30)
