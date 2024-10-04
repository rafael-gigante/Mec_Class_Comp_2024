include("ForcedPendulumSim.jl")
import .ForcedPendulumSim as FP

###############################################
#                 Variables                   #
############################################### 

g = 9.81
m = 2
l = 1
F_0 = 20
gamma = sqrt(g/l) / 4

omega_r = sqrt(g/l)     # omega_0 from the slides (natural frequency of the system)
OMEGA = 0.66 * omega_r   # Frequency of external force
beta = gamma/(2*m*l)    # Damping coefficient, dimensionless
f_0 = F_0 / (m*g)       # Amplitude of the driving force, dimensionless
theta_0 = 0.1           # Initial angle (radians)
omega_0 = 0             # Initial angular velocity

pendulum = FP.ForcedPendulum(OMEGA, omega_r, beta, f_0, theta_0, omega_0)

###############################################
#                 Main Code                   #
############################################### 

# Normal simulation 
# Generates 3 graphs: Theta evolution, Phase diagram and Poincaré section
FP.run_simulation(pendulum, (0.0, 600.0), 0.01)

# Study of the bifurcation diagram (Takes much longer)
F_range = 1.060:0.00001:1.085
t_final = 600.0           
dt = 0.01                 
FP.generate_bifurcation_diagram(pendulum, F_range, t_final, dt)