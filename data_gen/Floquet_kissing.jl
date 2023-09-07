include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using CSV 
using DifferentialEquations
using ProgressBars

function f(u, p, t)
    # H, H_0, H_d, ω_d = p[1], p[2], p[3], p[4] 
    H = p[1] + p[2]*cos(p[3]*t)
    return -im * H * u
end

function f!(du,u, p, t)
    # H, H_0, H_d, ω_d = p[1], p[2], p[3], p[4] 
    H = p[1] + p[2]*cos(p[3]*t)
    du .= -im * H * u
end

function quasienergies(N, ω_0, g_n, Ω_d, ω_d)
    p =  H_0(N,ω_0,g_n), H_d(N,Ω_d), ω_d
    T = 2*pi/ω_d
    u_0 = ComplexF64.(Matrix(I(N)))
    tspan = (0.0, 2*T)
    prob = ODEProblem(f!, u_0, tspan, p)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
    U_2T = sol.u[end]
    η_n, _ = eigen(U_2T)
    ϵ_n = im*log.(η_n)/T
    return real.(ϵ_n)
end

N = 30 #from the calculations in the eff model, N>=25 is sufficient 
ω_0 = 1.
g_n = [10^(-5),10^(-6)]
ω_d = 2*ω_0

Ω_d_array = Vector(0:0.005:5) #graph above could be truncated at 10

# Define data form
data_array = zeros( N*length(Ω_d_array), 3 ) # data form: Δϵ_n | Ω_d | N

counter = 1
for (i,Ω_d) in ProgressBar(enumerate(Ω_d_array))
    lamb = quasienergies(N, ω_0, g_n, Ω_d, ω_d)
    for n=1:N
        data_array[counter,:] .= lamb[n]-minimum(lamb), Ω_d, N
        counter += 1
    end
end

# Put data in convenient DataFrame object
df_floquet = DataFrame(data_array, ["Δϵ_n","Ω_d","N"]) 
CSV.write("data/floquet_kissing.csv", df_floquet)