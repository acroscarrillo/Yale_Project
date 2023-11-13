include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LaTeXStrings # latex support for sexy figs
using DataFrames # this is like pandas
using CSV 
using Plots

# load data file
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))


# Define parameter space

ϵ_1_array_h_eff = unique(df.ϵ_1)
ϵ_2_array_h_eff = unique(df.ϵ_2)

# ϵ_2_2_save = Vector(range(0.5,10,length=20)) #save these screenshots

#plotting stuff
# linea_negra(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)
# linea_fija(ϵ_1,ϵ_2) = 4*ϵ_2
# linea_primera(ϵ_1,ϵ_2) = 4*ϵ_1*sqrt(ϵ_2)
@gif for ϵ_2 in ProgressBar(ϵ_2_array)
    ΔE_plots = zeros(length(ϵ_1_array), length(n_levels))
    for n in n_levels
        n = Int(n)
        ΔE_plots[:,n] = filter(row ->  row.ϵ_2 == ϵ_2 && row.n == n, df).ΔE_n
    end
    scatter(ϵ_1_array, ΔE_plots,markersize=1,markerstrokewidth=0,title="H_eff, N="*string(N[1])*", n_levels="*string(maximum(n_levels))*", "*L"\epsilon_2="*string(ϵ_2),xlabel=L"\epsilon_1/K",ylabel=L"(E_n-E_0)/K",legend=false,ylim=(-5, 150))
    plot!(ϵ_1_array, linea_negra.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    plot!(ϵ_1_array, linea_primera.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 2*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 3*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)
    # plot!(ϵ_1_array, 4*linea_fija.(ϵ_1_array,ϵ_2),linestyle=:dash,color=:black)

    if ϵ_2 in ϵ_2_2_save
        # savefig("figs/animations/anticrossings_frames/epsilon_2_"*string(ϵ_2)*".pdf")
    end
end