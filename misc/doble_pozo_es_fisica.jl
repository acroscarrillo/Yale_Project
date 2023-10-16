include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames
using CSV
using LaTeXStrings # latex support for sexy figs
using Plots
using QuantumOptics


#############
# load data #
#############

df_heatmap = DataFrame(CSV.File("data/h_eff_heatmap_crossing.csv"))
df_crossings = DataFrame(CSV.File("data/h_eff_crossings_epsilon_1.csv"))

df_crossings = filter(row ->  row.n > -1, df_crossings) # discard theory line

N = Int(unique(df_heatmap.N)[1])
Δ = unique(df_heatmap.Δ)[1]
K = unique(df_heatmap.K)[1]

ϵ_1_array = unique(df_heatmap.ϵ_1)
ϵ_2_array = unique(df_heatmap.ϵ_2)

cross_loc = [(1.8,4),(1.95,4), (2,4),(2.05,4)]
cross_labels = ["A","B","C"]
matrix_color = reshape(df_heatmap."min(Δnn)",(length(ϵ_2_array),length(ϵ_1_array)))

basis = FockBasis(N-1)
#############
# plot data #
#############

l = @layout [ a b c d ; e f g h]
plots_array = []




# STATES
y_lab = L"p"
x_lab = L"x"

crossing_states = [(2,3),(2,3),(2,3),(2,3)]
for j=1:2
    for i=1:length(cross_loc)
        ϵ_1,ϵ_2 = cross_loc[i]

        E, v = eigen(-H_eff(N,Δ,K,ϵ_1,ϵ_2))

        n = crossing_states[i][j]

        ttl = L"n="*string(n-1)*", "*L"E="*string(round(E[n],digits = 1))*", "*L"\epsilon_1="*string(ϵ_1)*", "*L"\epsilon_2="*string(ϵ_2)
        mat_temp = wigner( Ket(basis, v[:,n]), -7:0.1:7, -3:0.1:3)' # transpose to put p on y-axis

        push!(plots_array,heatmap(-7:0.1:7, -3:0.1:3,mat_temp,title=ttl,titlefontsize=6,colorbar=false, c = :berlin,xtickfontsize=6,ytickfontsize=6,clim=(-0.3,0.3),xticks=([-6,0,6],["-6","0","6"]),yticks=([-3,0,3],["-3","0","3"]),xlabel=x_lab, ylabel=y_lab,guidefont=font(8)))
        end
end
    
    

    
# PLOT IT!
ttl = L"H_{eff}"*" resonances at "*L"N="*string(N)*", "*L"\Delta="*string(Δ)*", "*L"K="*string(K)*".\n Color bar in log scale. All levels considered."

global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -75Plots.px, titlefontsize=10)

plot(plots_array..., layout = l)