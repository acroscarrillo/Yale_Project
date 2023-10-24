using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs


df = DataFrame(CSV.File("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv"))

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

################
# Reso Heatmap #
################
combined_mat = zeros((length(ϵ_2_array),length(ϵ_1_array)))
for n=2:7 # ignore last for the sake of plotting
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array)))
    combined_mat += temp_mat.^-1
end
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="Numerical theory, ILAS. "*L"1/\log | \delta_{n,n+1}|",clim=(-5,5),titlefontsize=8,xtickfontsize=8,ytickfontsize=8,guidefont=font(10),height=400,c = :oleron10,dpi=600)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

savefig("figs/important_figs/H_eff/heatmap_plaga.png")