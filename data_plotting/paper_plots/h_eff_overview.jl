include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

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

cross_loc = [(2,4), (4.5,9), (7.95,7)]
cross_labels = ["A","B","C"]
matrix_color = reshape(df_heatmap."min(Δnn)",(length(ϵ_2_array),length(ϵ_1_array)))

basis = FockBasis(N-1)
#############
# plot data #
#############

l = @layout [ttitle; a ; b c d ; e f g; h i j]
plots_array = []


# HEATMAP + CROSSES
y_lab = L"\epsilon_2/K"
x_lab = L"\epsilon_1/K"

heatmap_plot = heatmap(ϵ_1_array,ϵ_2_array, log.(matrix_color),colormap = :lajolla,colorbar=true,ylabel=y_lab,xlabel=x_lab,xtickfontsize=6,ytickfontsize=6,guidefont=font(8))

plot!([2,4.5,7.95], seriestype="vline",color="black")
annotate!([2,4.5,7.95].+0.1, 1,text.(cross_labels, :black, :left, 6))

# plot!([4,9,7], seriestype="hline",color="black",linestyle=:dash)



push!(plots_array,scatter!(cross_loc,marker=:xcross,color=:blue,markersize=7,legend=false))


# CROSSINGS
y_lab = L"(E-E_0)/K"
x_lab = L"\epsilon_1/K"

levels = unique(df_crossings.n)
counter = 1
for (ϵ_1_fixed,ϵ_2) in cross_loc
    ttl = L"\epsilon_2="*string(ϵ_2)

    df_temp = filter(row ->  row.ϵ_2 == ϵ_2, df_crossings)

    df_temp_2 = filter(row ->  row.n == 1, df_temp)
    plot_2_push = scatter(df_temp_2.ΔE_n,df_temp_2.ϵ_1) # to trigger the plot, note the !
    plot!([ϵ_1_fixed], seriestype="vline",color="black")
    annotate!(ϵ_1_fixed+0.5, 10,text(cross_labels[counter], :black, :left, 6))
    counter += 1
    for lvls in levels[2:end]
        df_temp_2 = filter(row ->  row.n == lvls, df_temp)

        plot_2_push = scatter!(df_temp_2.ϵ_1, df_temp_2.ΔE_n,markerstrokewidth=0,markersize=.4,titlefontsize=8,xtickfontsize=6,ytickfontsize=6,legend=false,label= string(lvls),xlabel=x_lab, ylabel=y_lab,title=ttl,guidefont=font(8))

    end
    push!(plots_array,plot_2_push)
end


# STATES
y_lab = L"p"
x_lab = L"x"

crossing_states = [(2,3),(3,4),(4,5)]
for j=1:2
    for i=1:length(cross_loc)
        ϵ_1,ϵ_2 = cross_loc[i]

        E, v = eigen(-H_eff(N,Δ,K,ϵ_1,ϵ_2))

        n = crossing_states[i][j]

        ttl = L"n="*string(n-1)*", "*L"E="*string(round(E[n],digits = 1))*", "*L"\epsilon_1="*string(ϵ_1)*", "*L"\epsilon_2="*string(ϵ_2)
        mat_temp = wigner( Ket(basis, v[:,n]), -7:0.1:7, -7:0.1:7)' # transpose to put p on y-axis

        push!(plots_array,heatmap(-7:0.1:7, -7:0.1:7,mat_temp,title=ttl,titlefontsize=6,colorbar=false, c = :berlin,xtickfontsize=6,ytickfontsize=6,clim=(-0.3,0.3),xticks=([-6,0,6],["-6","0","6"]),yticks=([-6,0,6],["-6","0","6"]),xlabel=x_lab, ylabel=y_lab,guidefont=font(8)))
        end
end
    
    

    
# PLOT IT!
ttl = L"H_{eff}"*" resonances at "*L"N="*string(N)*", "*L"\Delta="*string(Δ)*", "*L"K="*string(K)*".\n Color bar in log scale. All levels considered."

global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -75Plots.px, titlefontsize=10)

plot(global_title,plots_array..., layout = l)

savefig("figs/important_figs/H_eff/H_eff_overview.png")
savefig("figs/important_figs/H_eff/H_eff_overview.pdf")
