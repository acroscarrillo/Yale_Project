using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs


df = DataFrame(CSV.File("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv"))

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

l = @layout [ttl; a{.5h}; a b c; d e f; g h i ] 

plots = []

################
# Reso Heatmap #
################
combined_mat = zeros((length(ϵ_2_array),length(ϵ_1_array)))
for n=2:7 # ignore last for the sake of plotting
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array)))
    combined_mat += temp_mat.^-1
end
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="All plots bellow summed",clim=(-5,5),titlefontsize=8,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),height=400,c = :oleron100)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

plot!([5.9,7.6,9.7], seriestype="hline",color="black",label=false,style=:dash)

push!(plots,heatmap_combined)

##########
# Levels #
##########
for n=2:7 # ignore last for the sake of plotting
    nn = n+1
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array))).^-1
    heatmap_tmp = heatmap(ϵ_1_array,ϵ_2_array,temp_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="levels $n and $nn",clim=(-5,5),titlefontsize=7,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),xticks=([0,12],["0","12"]),yticks=([0,12],["0","12"]),c = :oleron100, colorbar=false)
    push!(plots,heatmap_tmp)
end

############
# Plot it! #
############
ttl = "Logarithmic inverse of anti-crossings spacing, "*L"1/\log|\delta_{n,n+1}|"*" of "*L"H_{eff}."
global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -55Plots.px, titlefontsize=8)

plot(global_title,plots..., layout = l,size=(720,400),dpi=600)


# savefig("figs/important_figs/H_eff/H_eff_anticrossings_oscillations.png")
# savefig("figs/important_figs/H_eff/H_eff_anticrossings_oscillations.pdf")