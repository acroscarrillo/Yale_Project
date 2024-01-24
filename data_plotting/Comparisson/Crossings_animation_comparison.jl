include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LaTeXStrings # latex support for sexy figs
using DataFrames # this is like pandas
using CSV 
using Plots

# load data file
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))

# Define parameter space (for reasons they are not numerically identical)
ϵ_1_h_eff_array = unique(df_h_eff.ϵ_1)
ϵ_2_h_eff_array = unique(df_h_eff.ϵ_2)

ϵ_1_floq_array = unique(df_floq.ϵ_1)
ϵ_2_floq_array = unique(df_floq.ϵ_2)

# filter out large levels
df_h_eff = filter(row ->  row.ΔE_n < 600, df_h_eff)
df_floq = filter(row ->  row.ΔE_n < 600, df_floq)

# filter out large n_photons
n_photons_cutoff = 20
df_h_eff = filter(row ->  row.n_photons < n_photons_cutoff, df_h_eff)
df_floq = filter(row ->  row.n_photons < n_photons_cutoff, df_floq)

# axis labels
ylbl = L"\epsilon_2/K"
xlbl = L"\epsilon_1/K"
ttl = "n_photons_cutoff=$n_photons_cutoff"

# ϵ_2_2_save = Vector(range(0.5,10,length=20)) #save these screenshots

#plotting stuff
lemnis_surf(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)
# linea_fija(ϵ_1,ϵ_2) = 4*ϵ_2
# linea_primera(ϵ_1,ϵ_2) = 4*ϵ_1*sqrt(ϵ_2)
@gif for i=ProgressBar(1:length(ϵ_2_h_eff_array))
    ϵ_2_h_eff = ϵ_2_h_eff_array[i]
    ϵ_2_floq = ϵ_2_floq_array[i]

    max_ΔE = 1.25*maximum( lemnis_surf.(ϵ_1_array,ϵ_2_h_eff) ) # lemniscate max height

    df_h_eff_temp = filter(row ->  row.ϵ_2 == ϵ_2_h_eff, df_h_eff)
    df_floq_temp = filter(row ->  row.ϵ_2 == ϵ_2_floq, df_floq)

    scatter(df_h_eff_temp.ϵ_1, df_h_eff_temp.ΔE_n,ylim=(-5,max_ΔE),markersize=1,markerstrokewidth=0,label="H_eff",legend=:bottomright,xlabel=xlbl,ylabel=ylbl)
    scatter!(df_floq_temp.ϵ_1, df_floq_temp.ΔE_n,ylim=(-5,max_ΔE),markersize=1,markerstrokewidth=0,label="Floq",title=ttl)
   
    # plot!(ϵ_1_array, linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 2*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 3*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 4*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)

    # if ϵ_2 in ϵ_2_2_save
        # savefig("figs/animations/anticrossings_frames/epsilon_2_"*string(ϵ_2)*".pdf")
    # end
end