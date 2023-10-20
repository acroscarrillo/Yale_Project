using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs
using NPZ

# Define plot structure
hlines_labels = [L"B'",L"A'",L"B",L"A"]
h_lines_pos = [5.9,7.6,9.7,11.25]
cross_loc = [(2.45,5.9),(5.5,7.6),(3.15,9.7),(9.3,9.7),(7.3,5.9),(6.75,11.25),(13.4,11.25) ]

# Experiment
exp_data = npzread("data/experimental/DecayData.npz")

K_exp = 0.505 # MHz or 505KHz

heatmap_mat_exp = exp_data["arr_2"] # μs, but irrelant for now
ϵ_1_array_exp = exp_data["arr_0"][1,:] # MHz
ϵ_2_array_exp = exp_data["arr_1"][:,1] # MHz

ϵ_1_array_exp = ϵ_1_array_exp./K_exp # Units of Kerr
ϵ_2_array_exp = ϵ_2_array_exp./K_exp # Units of Kerr

# I sliced the speckels off
heatmap_exp = heatmap(ϵ_1_array_exp[1:end-20],ϵ_2_array_exp[1:17],heatmap_mat_exp[1:17,1:end-20],xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="Experiment, life-time.",titlefontsize=12,xtickfontsize=8,ytickfontsize=8,guidefont=font(12),height=400,colorbar=true,dpi=600, size=(720,400))
plot!(h_lines_pos, seriestype="hline",color="white",label=false,style=:dash)
annotate!(1,h_lines_pos .+ .8,text.(hlines_labels, :black, :left, 10))
scatter!(cross_loc,marker=:xcross,color=:white,markersize=6,legend=false)


savefig("figs/important_figs/H_eff/H_eff_experiment_resonances_oscillations.png")
savefig("figs/important_figs/H_eff/H_eff_experiment_resonances_oscillations.pdf")


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

ttl = "Numerical theory, ILAS. "*L"("*"Inverse logarithmic anti-crossings spacing "*L"1/\log|\delta_{n,n+1}|)"
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title=ttl,clim=(-10,3),titlefontsize=11,xtickfontsize=8,ytickfontsize=8,guidefont=font(12),height=400,colorbar=true,yticks=([2.5,5,7.5,10,12.5]),dpi=600, size=(720,400))#,yticks=([0,10],["0","10"])
h_lines_pos = range(2.5,12.5,length=10)
v_lines_pos = range(0.5,13,length=11)
plot!(h_lines_pos, seriestype="hline",color="black",label=false,style=:dash,alpha=0.3)
plot!(v_lines_pos, seriestype="vline",color="black",label=false,style=:dash,alpha=0.3)
# annotate!(1,h_lines_pos .+ .8,text.(hlines_labels, :black, :left, 12))
# scatter!(cross_loc,marker=:xcross,color=:limegreen,markersize=6,legend=false)


savefig("figs/important_figs/H_eff/H_eff_theory_resonances_oscillations.png")
savefig("figs/important_figs/H_eff/H_eff_theory_resonances_oscillations.pdf")