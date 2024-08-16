using Plots

# Define parameters
g = 9.81 # Acceleration due to gravity (m/s²)
L = 1.0  # Length of the baton (m)
θ₀ = 0.0 # Initial angle (radians)
ω₀ = -sqrt(50.0) # Initial angular velocity (rad/s)
r₀ = [0.0, 1.0] # Initial position of the center of mass (m)
v₀ = [5.0, 5.0] # Initial velocity of the center of mass (m/s)
t_final = 2.0 # Simulation time (s)
dt = 0.01 # Time step (s)

# Define time array
t = 0:dt:t_final

# Initialize arrays to store positions
r_cm = [r₀ .+ v₀ * ti .- [0, 0.5 * g * ti^2] for ti in t]
θ = [θ₀ + ω₀ * ti for ti in t]

# Compute positions of the ends of the baton
r_left = [r .- 0.5L * [cos(θi), sin(θi)] for (r, θi) in zip(r_cm, θ)]
r_right = [r .+ 0.5L * [cos(θi), sin(θi)] for (r, θi) in zip(r_cm, θ)]

x_cm = []
y_cm = []
for i in 1:length(t)
    push!(x_cm, r_cm[i][1])
    push!(y_cm, r_cm[i][2])
end

# Plot the motion
anim = @animate for i in 1:length(t)
    plot([r_left[i][1], r_right[i][1]], [r_left[i][2], r_right[i][2]], lw=2, xlims=(-5, 10), ylims=(0, 5), color=:blue, label=false)
    scatter!([r_left[i][1], r_right[i][1]], [r_left[i][2], r_right[i][2]], color=:red, label=false)
    plot!([x_cm[1:i]], [y_cm[1:i]], color=:black, linestyle =:dot, label=false)
end


gif(anim, "throwing_baton.gif", fps=20)
