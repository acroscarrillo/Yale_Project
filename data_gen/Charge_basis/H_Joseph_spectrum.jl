include("../../src/main_charge_basis.jl") # import charge basis stuff

N = 101
E_J = 1 #in units of E_J
m = 3

α = 0.1
ϕ_ext_array = Vector( 0:0.1:0.001)

spectrum_array = zeros(N,length(ϕ_ext_array))
for (n,ϕ) in enumerate(ϕ_ext_array)
    spectrum_array[n,:] = eigen(H_joseph(N,α,E_J,ϕ,m)).values
end

heatmap(spectrum_array)


α_array = Vector( 0:0.01:1)
ϕ_ext = 1

spectrum_array = zeros(N,length(α_array))
for (n,α) in enumerate(α_array)
    spectrum_array[n,:] = eigen(H_joseph(N,α,E_J,ϕ_ext,m)).values
end

heatmap(spectrum_array)
