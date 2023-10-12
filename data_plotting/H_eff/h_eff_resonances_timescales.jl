using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs

# Define plot structure
l = @layout [ttl; experiment; cut_1 cut_2 ] #a{.5h}
plots = []
hlines_labels = [L"B'",L"A'",L"B",L"A"]
h_lines_pos = [5.9,7.6,9.7,11.25]
cross_loc = [(2.45,5.9),(5.5,7.6),(3.15,9.7),(9.3,9.7),(7.3,5.9),(6.75,11.25),(13.4,11.25) ]

# Theory
df_theo = DataFrame(CSV.File("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv"))

ϵ_1_array = unique(df_theo.ϵ_1)
ϵ_2_array = unique(df_theo.ϵ_2)

combined_mat = zeros((length(ϵ_2_array),length(ϵ_1_array)))
for n=1:7 # ignore last for the sake of plotting
    nn = n+1
    df_2_plot = filter(row ->   row.cross_n == n, df_theo)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array)))
    combined_mat += temp_mat.^-1
end

ttl = "Theory, ILAS. "*L"("*"Inverse logarithmic anti-crossings spacing "*L"1/\log|\delta_{n,n+1}|)"
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title=ttl,clim=(-10,3),titlefontsize=10,xtickfontsize=7,ytickfontsize=7,guidefont=font(10),height=400,colorbar=false,yticks=([2.5,5,7.5,10,12.5]))#,yticks=([0,10],["0","10"])
plot!(h_lines_pos, seriestype="hline",color="black",label=false,style=:dash)
annotate!(1,h_lines_pos .+ .8,text.(hlines_labels, :black, :left, 10))
scatter!(cross_loc,marker=:xcross,color=:limegreen,markersize=6,legend=false)

push!(plots,heatmap_combined)

# ϵ_2 cuts of this plot
# A cut:

# B cut:


# Plot it all
ttl = "Experiment vs. theory of "*L"H_{eff}."
global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -150Plots.px, titlefontsize=13)

plot(global_title,plots..., layout = l)

# savefig("figs/important_figs/H_eff/H_eff_vs_experiment_resonances_oscillations.png")
# savefig("figs/important_figs/H_eff/H_eff_vs_experiment_resonances_oscillations.pdf")