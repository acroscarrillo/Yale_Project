const ħ = 6.62607015e-34
const ec = 1.602176634e-19 
const h = 2*π * ħ
const ϕ_0 = ħ/(2*ec)
#########
# Snail #
#########

# struct and constructor
struct Snail
    α 
    n
    ϕ_E
    ϕ_array
    i_min
    c2
    c3
    c4
    c5
    c6
end

function Snail(α, n, ϕ_E)
    ϕ_array = Vector(range(-3*π,3*π,length=10001))
    i_min =  argmin(U.(α,ϕ_E,n,ϕ_array))
    c2 = U_der2(α,ϕ_E,n,ϕ_array[i_min])  
    c3 = U_der3(α,ϕ_E,n,ϕ_array[i_min])
    c4 = U_der4(α,ϕ_E,n,ϕ_array[i_min])
    c5 = U_der5(α,ϕ_E,n,ϕ_array[i_min])
    c6 = U_der6(α,ϕ_E,n,ϕ_array[i_min])

    return Snail(α, n, ϕ_E, ϕ_array, i_min, c2, c3, c4, c5, c6)
end

# potentials for calculating Taylor expansion
U(α,ϕ_E,n,ϕ) = -α * cos(ϕ) - n * cos((ϕ_E-ϕ)/n)
U_der(α,ϕ_E,n,ϕ) = α * sin(ϕ) - sin((ϕ_E-ϕ)/n)
U_der2(α,ϕ_E,n,ϕ) = α*cos(ϕ) + (1/n)*cos((ϕ_E-ϕ)/n)
U_der3(α,ϕ_E,n,ϕ) = -α*sin(ϕ) + (1/n^2)*sin((ϕ_E-ϕ)/n)
U_der4(α,ϕ_E,n,ϕ) = -α*cos(ϕ) - (1/n^3)*cos((ϕ_E-ϕ)/n)
U_der5(α,ϕ_E,n,ϕ) = α * sin(ϕ) - (1/n^4)*sin((ϕ_E - ϕ)/n)
U_der6(α,ϕ_E,n,ϕ) = α * cos(ϕ) + (1/n^5)*cos((ϕ_E - ϕ)/n)


###########
# Circuit #
###########
struct Circuit
    L
    C
    n_modes
    ω
    Z
    ϕ_zpf
    q_zpf
    E_L
    E_C
    φ_zpf
    n_zpf
end

function Circuit(L, C, n_modes=1)
    ω = 1 / √(L*C)
    Z = 1 / ω
    ϕ_zpf = √(ħ * Z / 2)
    q_zpf = √(ħ / (2 * Z))
    E_L = (ϕ_0^2)/(L * ħ)
    E_C = (ec^2)/(2 * C * ħ)
    φ_zpf = (2 * E_C / E_L)^0.25
    n_zpf = (E_L / (32 * E_C))^0.25
    return Circuit(L, C, n_modes, ω, Z, ϕ_zpf, q_zpf, E_L, E_C, φ_zpf, n_zpf)
end

function Base.show(io::IO, circuit::Circuit)
    print(io,"'Circuit' type with properties: \n")
    print(io,"f0 = "*string(round(circuit.ω / (2*π) * 1e-9, sigdigits=3))," GHz \n")
    print(io,"L_tot = "*string(round(circuit.L * 1e-9, sigdigits=3))," nH \n")
    print(io,"C_tot = "*string(round(circuit.C * 1e12, sigdigits=3))," pF \n")
    print(io,"Z = "*string(round(circuit.Z, sigdigits=3))," Ohm \n")
    print(io,"E_L = "*string(round(circuit.E_L / (2*π)  * 1e-9, sigdigits=3))," GHz \n")
    print(io,"E_L/E_C = "*string(round(circuit.E_L / circuit.E_C, sigdigits=3))," \n")
    print(io,"φ_zpf = "*string(round(circuit.φ_zpf, sigdigits=3)))
end

    
################
# SnailCircuit #
################
struct SnailCircuit
    L_j
    snail
    M
    E_J
    L_s
    p
    c2
    c3
    c4
    c5
    c6
    g3
    g4
    g5
    g6
    g4_eff
end

function SnailCircuit(L_lin, C, L_j, snail=Snail(0.1, 3, π), M=1)
    cir = Circuit(L_lin + M*L_j/snail.c2, C) # is this just to call the pretty print?? lol
    
    E_J = (ϕ_0^2)/(L_j * ħ) #TODO - is this a problem? shouldn't it be L_j?
    L_s = L_j / snail.c2
    p = M * L_s / cir.L
    c2 = p * snail.c2 / M
    c3 = p^3 * snail.c3 / M^2
    c4 = p^4 / M^3 * (snail.c4 - 3*snail.c3^2 / snail.c2 * (1-p))
    
    temp = snail.c5 - 6*snail.c3*snail.c4 / snail.c2 * (1-p) + 5 * snail.c3^3 / (snail.c2^2) * (1-p)^2
    c5 = temp * p^5 / M^4 
    
    temp = snail.c6 - (11 * snail.c3 * snail.c5 + 6 * snail.c4^2) / snail.c2 * (1-p)
    temp += (51 * snail.c3^2 * snail.c4) / snail.c2^2 * (1-p)^2
    temp -= 35 * snail.c3^4 / snail.c2^3 * (1-p)^3
    c6 = temp * p^6 / M^5 
    
    
    g3 = c3 * E_J * cir.φ_zpf^3 / factorial(3) 
    g4 = c4 * E_J * cir.φ_zpf^4 / factorial(4) 
    g5 = c5 * E_J * cir.φ_zpf^5 / factorial(5) 
    g6 = c6 * E_J * cir.φ_zpf^6 / factorial(6) 
    g4_eff =  g4 - 5*(g3^2)/cir.ω
    return SnailCircuit(L_j, snail, M, E_J, L_s, p, c2, c3, c4, c5, c6, g3, g4, g5, g6, g4_eff)
end

function Base.show(io::IO, snailcircuit::SnailCircuit)
    print(io,"'SnailCircuit' type with properties: \n")
    print(io,"M = "*string(snailcircuit.M)," \n")
    print(io,"p = "*string(round(snailcircuit.p, sigdigits=3))," \n")
    print(io,"L_j = "*string(round(snailcircuit.L_j, sigdigits=3))," nH \n")
    print(io,"g3 = "*string(round(snailcircuit.g3 * 1e-6/(2*π), sigdigits=3))," MHz \n")
    print(io,"g4 = "*string(round(snailcircuit.g4 * 1e-6/(2*π), sigdigits=3))," MHz \n")
    print(io,"g5 = "*string(round(snailcircuit.g5 * 1e-3/(2*π), sigdigits=3))," kHz \n")
    print(io,"g6 = "*string(round(snailcircuit.g6 * 1e-3/(2*π), sigdigits=3))," kHz \n")
    print(io,"g4_eff = "*string(round(snailcircuit.g4_eff * 1e-6/(2*π), sigdigits=3))," MHz \n")
    print(io,"ξ_qq = "*string(round(12*snailcircuit.g4_eff * 1e-6/(2*π), sigdigits=3))," MHz \n")
end
    
function get_KC_summary(snailcircuit, nbar=2.0)
    print(snailcircuit)
    print("Kerr Cat for nbar = $nbar")
    print("K = ",round(6 * snailcircuit.g4 * 1e-6 / (2*π), sigdigits=3), " MHz \n")
    print("K_eff = ", round(6 * snailcircuit.g4_eff * 1e-6 / (2*π), sigdigits=3), " MHz \n")
    print("etc")
end


xi_for_cat(snailcircuit, nbar=2.0) = 2*nbar*snailcircuit.g4_eff / snailcircuit.g3

gap_cat(snailcircuit, nbar=2.0) = 4 * nbar * snailcircuit.g4_eff * 6.0

####################
# helper functions #
####################
calc_C(ω, L) = 1 / (ω^2 * L)
calc_L(ω, C) = 1 / (ω^2 * C)
