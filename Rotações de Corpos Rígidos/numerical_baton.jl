mutable struct VetorPosiçãoCM
    x::Float64
    y::Float64
end

mutable struct VetorVelocidadeCM
    v_x::Float64
    v_y::Float64
end

mutable struct VetorAceleraçãoCM
    a_x::Float64
    a_y::Float64
end

mutable struct Bastão
    L::Float64 # Tamanho do bastão
    m₁::Float64 # Massa da esfera 1
    m₂::Float64 # Massa da esfera 2
    r::VetorPosiçãoCM # Vetor posição do centro de massa
    v::VetorVelocidadeCM # Vetor velocidade do centro de massa
    a::VetorAceleraçãoCM # Vetor aceleração do centro de massa
    φ::Float64 # Ângulo de orientação do bastão
    ω::Float64 # Velocidade angular do bastão
    α::Float64 # Aceleração angular do bastão
end

function criar_bastão(L::Float64, m₁::Float64, m₂::Float64, r::VetorPosiçãoCM, v::VetorVelocidadeCM, a::VetorAceleraçãoCM, φ::Float64, ω::Float64, α::Float64)
    return Bastão(L, m₁, m₂, r, v, a, φ, ω, α)
end

function atualizar_ângulo!(bastão::Bastão, dt::Float64)
    # Atualiza o ângulo de orientação do bastão com base em sua velocidade angular
    bastão.φ += bastão.ω * dt
    bastão.ω += bastão.α * dt
end

function atualizar_posição!(bastão::Bastão, dt::Float64)
    # Atualiza a posição do centro de massa do bastão com base em sua velocidade
    bastão.r.x += bastão.v.v_x * dt
    bastão.r.y += bastão.v.v_y * dt
end

function atualizar_velocidade!(bastão::Bastão, dt::Float64)
    # Atualiza a velocidade do centro de massa do bastão com base em sua aceleração
    bastão.v.v_x += bastão.a.a_x * dt
    bastão.v.v_y += bastão.a.a_y * dt
end

function simular_movimento(bastão::Bastão, t_max::Float64, dt::Float64)
    t = 0.0
    while t < t_max
        atualizar_posição!(bastão, dt)
        atualizar_velocidade!(bastão, dt)
        atualizar_ângulo!(bastão, dt)
        
        println("Time: ", t, "s, Position: (", bastão.r.x , ", ", bastão.r.y, 
                "), Angle: ", bastão.φ, " radians")

        t += dt
    end
end



posição = VetorPosiçãoCM(0.0, 0.0)
velocidade = VetorVelocidadeCM(5.0, 10.0)
aceleração = VetorAceleraçãoCM(0.0, -9.81)
bastão = criar_bastão(1.0, 1.0, 1.0, posição, velocidade, aceleração, 0.0, 5.0, 0.0)

simular_movimento(bastão, 10.0, 0.1)


