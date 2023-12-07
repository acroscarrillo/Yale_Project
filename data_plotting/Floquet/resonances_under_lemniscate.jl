include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LaTeXStrings # latex support for sexy figs
using DataFrames # this is like pandas
using CSV 
using ProgressBars
using Plots

# load data file
df_floquet =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))

# focus on the first
df_floquet = filter(row ->  row.ΔE_n < 1000, df_floquet)

# Define parameter space

ϵ_1_array_floq = unique(df_floquet.ϵ_1)
ϵ_2_array_floq = unique(df_floquet.ϵ_2)

ylbl = L"\epsilon_2/K"
xlbl = L"\epsilon_1/K"
# Crossings
lemnis_surf(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)
reso_htmp_array = zeros(length(ϵ_2_array_floq),length(ϵ_1_array_floq))
for (i,ϵ_2) in ProgressBar(enumerate(ϵ_2_array_floq))
    df_temp_1 = filter(row ->  row.ϵ_2 == ϵ_2, df_floquet)
    for (j,ϵ_1) in enumerate(ϵ_1_array_floq)

        y_max = 1.15 * lemnis_surf(ϵ_1,ϵ_2)
        df_temp_2 = filter(row ->  row.ΔE_n < y_max && row.ϵ_1 == ϵ_1, df_temp_1)
        ΔE_lemnis = df_temp_2.ΔE_n
        
        if length(ΔE_lemnis) > 1
            gap_array = zeros(length(ΔE_lemnis)-1)
            for n=1:(length(ΔE_lemnis)-1)
                gap_array[n] = ΔE_lemnis[n+1]-ΔE_lemnis[n]
            end
            reso_htmp_array[i,j] = minimum(gap_array)
        end

    end
end

heatmap(ϵ_1_array_floq,ϵ_2_array_floq,log.(reso_htmp_array),xlab=xlbl,ylab=ylbl)