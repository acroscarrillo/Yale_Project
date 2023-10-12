include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using Plots 
using LaTeXStrings
using ProgressBars
using QuantumOptics

# Define parameter space
N = 100
Δ = 1e-3 #to break the parity degeneracy (anti-epilepsy)
K = 1
n_levels = 5

basis = FockBasis(N-1)

ϵ_1_array = [0,3.33,6.66,10]
ϵ_2_array = Vector( range(0, 10, length=200) )

ϵ_2_2_save = ϵ_2_array[1:20:200]

# Define plot

l = @layout [ttitle; a b c d e ; f g h i j ; k l m n o ; p q r s t]

global frame_counter = 1

anim = @gif for ϵ_2 in ProgressBar(ϵ_2_array)
    plots = []
    for ϵ_1 in ϵ_1_array
        E, v = eigen(-H_eff(N,Δ,K,ϵ_1,ϵ_2))
        if ϵ_1==0
            sorted_levels = 1:n_levels
        else
            X_exp = [v[:,n]'*(a(N)+a(N)')*v[:,n] for n=1:n_levels]
            sorted_levels = sortperm(X_exp,rev=true)
        end
        for n in sorted_levels
            ttl = L"\epsilon_1="*string(ϵ_1)*", "*L"n="*string(n-1)*", "*L"E="*string(round(E[n],digits = 1))
            mat_temp = wigner( Ket(basis, v[:,n]), -7:0.1:7, -7:0.1:7)' # transpose to put p on y-axis
            push!(plots,heatmap(-7:0.1:7, -7:0.1:7,mat_temp,title=ttl,titlefontsize=6,colorbar=false, c = :berlin,xtickfontsize=6,ytickfontsize=6,clim=(-0.3,0.3),xticks=([-6,0,6],["-6","0","6"]),yticks=([-6,0,6],["-6","0","6"])))
        end
    end
    title = plot(title = "H_eff eigenstates at N=$N, "*L"\epsilon_2="*string(round(ϵ_2,digits = 1)), grid = false, showaxis = false, bottom_margin = -50Plots.px)
    plot(title,plots..., layout = l)

    if ϵ_2 in ϵ_2_2_save
        savefig("figs/animations/states_animation_frames/"*string(ϵ_2)*".pdf")
    end
end