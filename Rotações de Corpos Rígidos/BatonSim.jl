module BatonSim

export Bastão, ExtremidadeBastão, Vetor2D, simular_movimento, criar_bastão

mutable struct Vetor2D
    x::Float64
    y::Float64
end

mutable struct ExtremidadeBastão
    posição::Vetor2D # Posição da esfera
    distância_cm::Vetor2D # Vetor com a distância entre o CM e a esfera
    massa::Float64 # Massa da esfera
    sinal::Float64 # Determina se a massa está na esquerda ou na direita (-1 ou 1)
end

mutable struct Bastão
    comprimento::Float64 # Tamanho do bastão
    extremidades::Tuple{ExtremidadeBastão, ExtremidadeBastão}
    posição::Vetor2D # Vetor posição do centro de massa
    velocidade::Vetor2D # Vetor velocidade do centro de massa
    aceleração::Vetor2D # Vetor aceleração do centro de massa
    ângulo::Float64 # Ângulo de orientação do bastão
    velocidade_angular::Float64 # Velocidade angular do bastão
    aceleração_angular::Float64 # Aceleração angular do bastão
end

function criar_bastão(L::Float64, m::Tuple{Float64, Float64}, posição::Vetor2D, velocidade::Vetor2D, aceleração::Vetor2D, φ::Float64, ω::Float64, α::Float64)
    # Calcula o módulo da distância de cada esfera até o centro de massa
    L1 = (m[2] / (m[1] + m[2])) * L
    L2 = (m[1] / (m[1] + m[2])) * L

    # Calcula o vetor distância entre o CM e a esfera
    distância_cm1 = calcular_vetor_distância_esfera(L1, φ, -1.0)
    distância_cm2 = calcular_vetor_distância_esfera(L2, φ, 1.0)

    # Calcula o vetor posição de cada esfera
    posição_ext1 = Vetor2D(posição.x + distância_cm1.x, posição.y + distância_cm1.y)
    posição_ext2 = Vetor2D(posição.x + distância_cm2.x, posição.y + distância_cm2.y)

    # Cria as duas extremidades do bastão
    esfera1 = ExtremidadeBastão(posição_ext1, distância_cm1, m[1], -1.0)
    esfera2 = ExtremidadeBastão(posição_ext2, distância_cm2, m[2], 1.0)

    # Cria e retorna o objeto bastão
    return Bastão(L, (esfera1, esfera2), posição, velocidade, aceleração, φ, ω, α)
end

function atualizar_ângulo_euler!(bastão::Bastão, dt::Float64)
    # Atualiza o ângulo de orientação do bastão com base em sua velocidade angular
    bastão.ângulo += bastão.velocidade_angular * dt
    bastão.velocidade_angular += bastão.aceleração_angular * dt
end

function atualizar_posição_bastão_euler!(bastão::Bastão, dt::Float64)
    # Atualiza a posição do centro de massa do bastão com base em sua velocidade
    bastão.posição.x += bastão.velocidade.x * dt
    bastão.posição.y += bastão.velocidade.y * dt
end

function atualizar_velocidade_euler!(bastão::Bastão, dt::Float64)
    # Atualiza a velocidade do centro de massa do bastão com base em sua aceleração
    bastão.velocidade.x += bastão.aceleração.x * dt
    bastão.velocidade.y += bastão.aceleração.y * dt
end

function atualizar_posição_esferas!(bastão::Bastão)
    # Atualiza as posições das esferas em relação ao centro de massa
    norm1 = sqrt(bastão.extremidades[1].distância_cm.x^2 + bastão.extremidades[1].distância_cm.y^2)
    bastão.extremidades[1].distância_cm = calcular_vetor_distância_esfera(norm1, bastão.ângulo, bastão.extremidades[1].sinal)
    norm2 =sqrt(bastão.extremidades[2].distância_cm.x^2 + bastão.extremidades[2].distância_cm.y^2)
    bastão.extremidades[2].distância_cm = calcular_vetor_distância_esfera(norm2, bastão.ângulo, bastão.extremidades[1].sinal)

    # Atualiza a posição das extremidades do bastão
    bastão.extremidades[1].posição.x = bastão.posição.x + bastão.extremidades[1].distância_cm.x
    bastão.extremidades[1].posição.y = bastão.posição.x + bastão.extremidades[1].distância_cm.y
    bastão.extremidades[2].posição.x = bastão.posição.x + bastão.extremidades[2].distância_cm.x
    bastão.extremidades[2].posição.y = bastão.posição.x + bastão.extremidades[2].distância_cm.y
end

function calcular_vetor_distância_esfera(distância::Float64, ângulo::Float64, sinal_extremidade::Float64)
    dx = sinal_extremidade * distância * cos(ângulo)
    dy = sinal_extremidade * distância * sin(ângulo)
    return Vetor2D(dx, dy)
end

function simular_movimento(bastão::Bastão, t_max::Float64, dt::Float64)
    t = 0.0
    tempo = [t]
    posição_cm = [Vetor2D(bastão.posição.x, bastão.posição.y)]
    posição_esf1 = [Vetor2D(bastão.extremidades[1].posição.x, bastão.extremidades[1].posição.y)]
    posição_esf2 = [Vetor2D(bastão.extremidades[2].posição.x, bastão.extremidades[2].posição.y)]
    while t < t_max
        atualizar_posição_bastão_euler!(bastão, dt)
        atualizar_velocidade_euler!(bastão, dt)
        atualizar_ângulo_euler!(bastão, dt)
        atualizar_posição_esferas!(bastão)
        
        t += dt
        push!(tempo, t)
        aux = Vetor2D(bastão.posição.x, bastão.posição.y)
        push!(posição_cm, aux)
        aux =  Vetor2D(bastão.extremidades[1].posição.x, bastão.extremidades[1].posição.y)
        push!(posição_esf1,aux)
        aux = Vetor2D(bastão.extremidades[2].posição.x, bastão.extremidades[2].posição.y)
        push!(posição_esf2,aux)
    end
    println(tempo)
    println(posição_cm)
end

end 