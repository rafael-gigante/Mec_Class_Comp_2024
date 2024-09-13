###############################################
#                  Packages and Files                   #
############################################### 

using Plots  # Import the Plots.jl package
include("gif_function.jl") # Import function that creates a gif


###############################################
#            Variable Declarations            #
############################################### 

# Constants
const SimulationTime = 4.0  # Simulation length in seconds
const TimeStep = 1e-2  # Size of the time step (seconds)
const Iterations = Int(SimulationTime / TimeStep)  # Number of iterations

# Initial angles (in radians)
const InitialTheta = 0.7
const InitialPhi = 0.0
const InitialPsi = 0.0

# Lists of angles
thetaValues = []
phiValues = []
psiValues = []

###############################################
#               Function Definitions          #
###############################################

# Function to generate a simple list of angles
function generateAngles()
    # Set initial values
    oldTheta = InitialTheta
    oldPhi  = InitialPhi
    oldPsi  = InitialPsi

    for _ in 1:Iterations
        # Generate new value
        newTheta = oldTheta + (rand() - 0.5) * TimeStep * 10
        newPhi = oldPhi + TimeStep * 5
        newPsi = oldPsi + TimeStep * 10

        # Append new value to the list
        push!(thetaValues, newTheta)
        push!(phiValues, newPhi)
        push!(psiValues, newPsi)

        # Update the old value
        oldTheta = newTheta
        oldPhi  = newPhi
        oldPsi  = newPsi
    end
end

# Function to run the simulation for all three angles and plot the results
function main()
    # Generating angles for theta, phi, and psi
    generateAngles()
    
    generate_gif(SimulationTime, TimeStep, thetaValues, phiValues, psiValues)
end

###############################################
#                 Main Program                #
###############################################

# Run the simulation
main()