module SpinningTopSim

mutable struct SpinningTop
    I::Vector{Float64} # Momentum of inertia
    M::Float64         # Mass
    𝑙::Float64         # Center of mass distance from the point (0,0)
    Ω::Float64         # Angular velocity of its physical axis
    θ₀::Float64
    θ::Float64         # Nutation angle 
    φ::Float64         # Precession angle
    ψ::Float64         # Rotation of the spinning top physical axis
end

function f(spinningtop::SpinningTop)
    a = (2 * spinningtop.M * 9.81 * spinningtop.𝑙) / spinningtop.I[1]
    b = (spinningtop.I[3] * spinningtop.Ω / spinningtop.I[1])^2
    f1 = a * (cos(spinningtop.θ₀) - cos(spinningtop.θ))
    f2 = - b * ((cos(spinningtop.θ₀) - cos(spinningtop.θ))^2 / sin(spinningtop.θ)^2)
    f = sqrt(f1 + f2)
    return f
end

end