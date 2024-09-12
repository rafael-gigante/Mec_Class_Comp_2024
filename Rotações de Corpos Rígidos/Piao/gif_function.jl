###############################################
#                  Packages                   #
############################################### 

using Plots  # Import the Plots.jl package

###############################################
#            Variable Declarations            #
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

    # Create the animation
    for t in t_values

        # Extract the angles at the current time
        t_index = Int(round(t / TimeStep))  # Find the closest index for time t
        t_index = min(max(t_index, 1), Int(SimulationTime / TimeStep))  # Ensure t_index stays within bounds
        theta = thetaValues[t_index]
        phi = phiValues[t_index]

        # Convert the spherical coordinates to Cartesian
        x, y, z = sphericalToCartesian(r, theta, phi)

        # Plot the piao in 3D
        plot3d([0, x], [0, y], [0, z], label=false, legend=false, 
                title="Evolução do Pião \nt = $(round(t, digits=1)) s",
                xlims=(-1, 1), ylims=(-1, 1), zlims=(0, 1), lw=2)
        scatter!([x], [y], [z], color=:red, ms=4)
        frame(anim)  # Capture the frame
    end
    
    return anim
end



# Plot the piao in 3D and generate an animation for every t in time and a graph of all the angles evolution
function generate_gif(SimulationTime, TimeStep, thetaValues, phiValues, psiValues)
    # Create time array for plotting
    timeValues = collect(0:TimeStep:SimulationTime-TimeStep)

    # Plot the angles evolution over time
    display(plot(timeValues, [thetaValues phiValues psiValues], 
        label=["Theta" "Phi" "Psi"], xlabel="Time (s)", 
        ylabel="Angle (rad)", title="Evolution of Angles over Time", linewidth=1))

    # Create an animation of the piao
    anim = animatePiao(SimulationTime, TimeStep, thetaValues, phiValues)
    gif(anim, "piao_evolution.gif")  # Save the animation as a GIF
end