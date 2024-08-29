module BatonSim

using Plots

export Vetor2D, simular_movimento, criar_bastão

mutable struct Vetor2D
    x::Float64
    y::Float64
end

mutable struct ExtremidadeBastão
    posição::Vetor2D # Posição da esfera
    distância_cm::Vetor2D # Vetor com a distância entre o CM e a esfera
    norm_distância::Float64 # Valor absoluto da distância entre o CM e a esfera
    massa::Float64 # Massa da esfera
    sinal::Float64 # Determina se a massa está na esquerda ou na direita, ou embaixo ou em cima (-1 ou 1)
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
    esfera1 = ExtremidadeBastão(posição_ext1, distância_cm1, L1, m[1], -1.0)
    esfera2 = ExtremidadeBastão(posição_ext2, distância_cm2, L2, m[2], 1.0)

    # Cria e retorna o objeto bastão
    return Bastão(L, (esfera1, esfera2), posição, velocidade, aceleração, φ, ω, α)
end

function calcular_vetor_distância_esfera(distância::Float64, ângulo::Float64, sinal_extremidade::Float64)
    dx = sinal_extremidade * distância * cos(ângulo)
    dy = sinal_extremidade * distância * sin(ângulo)
    return Vetor2D(dx, dy)
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
    norm1 = bastão.extremidades[1].norm_distância
    distância_cm1 = calcular_vetor_distância_esfera(norm1, bastão.ângulo, bastão.extremidades[1].sinal)
    bastão.extremidades[1].distância_cm = distância_cm1
    norm2 =bastão.extremidades[2].norm_distância
    distância_cm2 = calcular_vetor_distância_esfera(norm2, bastão.ângulo, bastão.extremidades[2].sinal)
    bastão.extremidades[2].distância_cm = distância_cm2

    # Atualiza a posição das extremidades do bastão
    posição_ext1 = Vetor2D(bastão.posição.x + distância_cm1.x, bastão.posição.y + distância_cm1.y)
    posição_ext2 = Vetor2D(bastão.posição.x + distância_cm2.x, bastão.posição.y + distância_cm2.y)
    bastão.extremidades[1].posição = posição_ext1
    bastão.extremidades[2].posição = posição_ext2
end

function medir_energia_trans(bastão::Bastão)
    # Massa total
    M = bastão.extremidades[1].massa + bastão.extremidades[2].massa
    v_x = bastão.velocidade.x
    v_y = bastão.velocidade.y
    K_trans = 0.5 * M * (v_x^2 + v_y^2)
    return K_trans
end

function medir_energia_rot(bastão::Bastão)
    # Momento de inércia
    I = (bastão.extremidades[1].massa * bastão.extremidades[1].norm_distância^2) + (bastão.extremidades[2].massa * bastão.extremidades[2].norm_distância^2)
    ω = bastão.velocidade_angular 
    K_rot = 0.5 * I * ω^2
    return  K_rot
end

function medir_energia_pot(bastão::Bastão)
    # Massa total
    M = bastão.extremidades[1].massa + bastão.extremidades[2].massa
    g = bastão.aceleração.y
    h = bastão.posição.y
    U = -1.0 * (M * g * h) 
    return U
end

function simular_movimento(bastão::Bastão, t_max::Float64, dt::Float64)
    t = 0.0
    # Arrays para armazenar as posições
    tempo = [t]
    aux = Vetor2D(bastão.posição.x, bastão.posição.y)
    posição_cm = [aux]
    aux =  Vetor2D(bastão.extremidades[1].posição.x, bastão.extremidades[1].posição.y)
    posição_esf1 = [aux]
    aux = Vetor2D(bastão.extremidades[2].posição.x, bastão.extremidades[2].posição.y)
    posição_esf2 = [aux]

    # Arrays para armazenar as energias
    k_trans = [medir_energia_trans(bastão)]
    k_rot = [medir_energia_rot(bastão)]
    U = [medir_energia_pot(bastão)]
    E = [k_trans[1] + k_rot[1] + U[1]]
    while t < t_max
        atualizar_velocidade_euler!(bastão, dt)
        atualizar_posição_bastão_euler!(bastão, dt)
        atualizar_ângulo_euler!(bastão, dt)
        atualizar_posição_esferas!(bastão)

        if bastão.extremidades[1].posição.y < 0.0 || bastão.extremidades[2].posição.y < 0.0
            println("Bastão atingiu o chão no tempo: ", t)
            break 
        end
        
        t += dt
        push!(tempo, t)
        aux = Vetor2D(bastão.posição.x, bastão.posição.y)
        push!(posição_cm, aux)
        aux =  Vetor2D(bastão.extremidades[1].posição.x, bastão.extremidades[1].posição.y)
        push!(posição_esf1, aux)
        aux = Vetor2D(bastão.extremidades[2].posição.x, bastão.extremidades[2].posição.y)
        push!(posição_esf2, aux)

        push!(k_trans, medir_energia_trans(bastão))
        push!(k_rot, medir_energia_rot(bastão))
        push!(U, medir_energia_pot(bastão))
        push!(E, k_trans[end] + k_rot[end] + U[end])
    end
    animar_trajetória(tempo, posição_cm, posição_esf1, posição_esf2)
    fazer_gráfico_energia(tempo, k_trans, k_rot, U, E)
end

function animar_trajetória(tempo::Vector{Float64}, posição_cm::Vector{Vetor2D}, posição_esf1::Vector{Vetor2D}, posição_esf2::Vector{Vetor2D})
    # Extraindo as componentes x e y de cada objeto
    x_cm = [ponto.x for ponto in posição_cm]
    y_cm = [ponto.y for ponto in posição_cm]
    x_esf1 = [ponto.x for ponto in posição_esf1]
    y_esf1 = [ponto.y for ponto in posição_esf1]
    x_esf2 = [ponto.x for ponto in posição_esf2]
    y_esf2 = [ponto.y for ponto in posição_esf2]

    # Cria a animação do lançamento do bastão
    animação = @animate for i in eachindex(tempo)
        plot(xlabel="Posição X", ylabel="Posição Y", minorgrid=true, title="Movimento do bastão", framestyle=:box, legend=:topright)
        plot!([x_cm[1:i]], [y_cm[1:i]], color=:blue, linestyle =:dot, label="Trajetória do CM")
        plot!([x_esf1[i], x_esf2[i]], [y_esf1[i], y_esf2[i]], lw=3, color=:blue, label=false)
        scatter!([x_esf1[i]], [y_esf1[i]], label=false, color=:red)
        plot!([x_esf1[1:i]], [y_esf1[1:i]], color=:red, linestyle =:dot, label="Trajetória da Esfera 1")
        scatter!([x_esf2[i]], [y_esf2[i]], label=false, color=:green)
        plot!([x_esf2[1:i]], [y_esf2[1:i]], color=:green, linestyle =:dot, label="Trajetória da Esfera 2")
        xlims!(x_cm[1]-1, maximum(x_cm)+2)
        ylims!(y_cm[1]-1, maximum(y_cm)+2)
        annotate!(0.55, 6.0, text("t = $(round(tempo[i], digits=2))", :black, 12))
    end
    gif(animação, "Rotações de Corpos Rígidos/simulação_bastão.gif", fps=40)
end 

function fazer_gráfico_energia(tempo::Vector{Float64}, k_trans::Vector{Float64}, k_rot::Vector{Float64}, U::Vector{Float64}, E::Vector{Float64})
    pl = plot(xlabel="Tempo", ylabel="Energia", minorgrid=true, title="Evolução da energia - Lançamento de bastão", framestyle=:box, legend=:topright)
    plot!(tempo, k_trans, label="Cinética")
    plot!(tempo, k_rot, label="Rotacional")
    plot!(tempo, U, label="Potencial")
    plot!(tempo, E, label="Total")

    savefig(pl, "Rotações de Corpos Rígidos/energia.png")
end 

end 