include("BatonSim.jl")
using .BatonSim

# Parâmetros iniciais do sistema
L = 1.0 # Comprimento do bastão 
m₁ = 1.0 # Massa da extremidade 1
m₂ = 1.0 # Massa da extremidade 2
posição = Vetor2D(0.0, 0.0) # Posição inicial do CM do bastão
velocidade = Vetor2D(7.5, 10.0) # Velocidade inicial do CM do bastão
aceleração = Vetor2D(0.0, -9.81) # Aceleração inicial do CM do bastão
φ = 0.0 # Orientação inicial do bastão (0: bastão na horizontal, π/2: bastão na vertical)
ω = 5.0 # Velocidade angular do bastão
α = 0.0 # Aceleração angular do bastão
t_max = 2.0 # Tempo máximo de simulação
dt = 0.01  # Intervalo de tempo


# Define o objeto bastão
bastão = criar_bastão(L, (m₁, m₂), posição, velocidade, aceleração, φ, ω, α)

simular_movimento(bastão, t_max, dt)