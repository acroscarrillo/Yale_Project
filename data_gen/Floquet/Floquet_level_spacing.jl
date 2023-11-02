include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars



df = DataFrame(CSV.File("data/floquet_resonances.csv"))

K = unique(df.K)
ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

r_array = zeros(length(ϵ_1_array),length(ϵ_2_array))

pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
for (j,ϵ_1) in enumerate(ϵ_1_array)
    for (k,ϵ_2) in enumerate(ϵ_2_array)
        df_temp = filter(row -> row.ϵ_1 == ϵ_1 && row.ϵ_2 == ϵ_2, df)
        Δnn_array = df_temp.Δnn
        r_array[j,k] = avg_spacing_ratio(Δnn_array)

        update(pbar)
    end
end

heatmap(ϵ_1_array,ϵ_2_array,r_array')