include("../../src/main_charge_basis.jl") # import charge basis stuff

N = 501
E_J = 1 #in units of E_J
E_C = 1/5907
m = 3
α = 0.1
ϕ_ext = π*2


E,V = eigen(H_joseph(N,E_J,α,ϕ_ext,m))

plot(norm.(V[:,1]))
plot(norm.(V[:,1:7]),legend=true)

left_E = []
middle_E = []
right_E = []
for (j,e) in enumerate(E)
    if isodd(j)
        push!(left_E,e)
    else
        push!(right_E,e)
    end
end
reverse!(left_E)
E_2_plot = vcat(left_E,right_E)
plot(E_2_plot)