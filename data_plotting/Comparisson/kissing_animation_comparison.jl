include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LaTeXStrings # latex support for sexy figs
using DataFrames # this is like pandas
using CSV 
using ProgressBars
using Plots

# load data file
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))

# focus on this first
df_h_eff = filter(row ->  row.ΔE_n < 700, df_h_eff)
df_floq = filter(row ->  row.ΔE_n < 700, df_floq)

# Define parameter space

ϵ_1_array_h_eff = unique(df_h_eff.ϵ_1)
ϵ_1_array_floq = unique(df_floq.ϵ_1)

ϵ_2_array_h_eff = unique(df_h_eff.ϵ_2)
ϵ_2_array_floq = unique(df_floq.ϵ_2)

ylbl = L"(E_n-E_0)/K \quad vs. \quad  (\epsilon_n-\epsilon_n_0)/K"

# Kissings
anim =  @animate for (i,ϵ__1) in ProgressBar(enumerate(ϵ_1_array_floq))
    df_temp_h_eff = filter(row ->  row.ϵ_1 == ϵ_1_array_h_eff[i], df_h_eff)
    df_temp_floq = filter(row ->  row.ϵ_1 == ϵ_1_array_floq[i], df_floq)

    # display(df_temp_floq.ΔE_n)

    scatter(df_temp_h_eff.ϵ_2, df_temp_h_eff.ΔE_n,ms=2,markerstrokewidth=0,title="ϵ_1=$ϵ__1",label="H_eff",ylabel=ylbl,xlabel=L"\epsilon_2/K")

    scatter!(df_temp_floq.ϵ_2, df_temp_floq.ΔE_n,ms=2,markerstrokewidth=0,label="Floq")
end
gif(anim, "anim_fps15.gif", fps = 10)

# Crossings
lemnis_surf(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)
linea_fija(ϵ_1,ϵ_2) = 4*ϵ_2
linea_primera(ϵ_1,ϵ_2) = 4*ϵ_1*sqrt(ϵ_2)
anim =  @animate for (i,ϵ__2) in ProgressBar(enumerate(ϵ_2_array_h_eff))
    y_max = 1.15*lemnis_surf(ϵ_1_array_h_eff[end],ϵ__2)
    df_temp_h_eff = filter(row ->  row.ϵ_2 == ϵ_2_array_h_eff[i], df_h_eff)
    scatter(df_temp_h_eff.ϵ_1, df_temp_h_eff.ΔE_n,ms=1,markerstrokewidth=0,title="ϵ_2=$ϵ__2",label="H_eff",ylabel=ylbl,xlabel=L"\epsilon_1/K")
    
    df_temp_floq = filter(row ->  row.ϵ_2 == ϵ_2_array_floq[i], df_floq)
    scatter!(df_temp_floq.ϵ_1, df_temp_floq.ΔE_n,ms=1,markerstrokewidth=0,label="Floq",ylim=(0,y_max))
end
gif(anim, "anim_fps15.gif", fps = 10)
