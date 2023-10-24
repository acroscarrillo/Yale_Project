using DataFrames # this is like pandas
using CSV 
using Plots 
using LaTeXStrings

df = DataFrame(CSV.File("data/h_eff_crossings_epsilon_1.csv"))

# Define parameter space
N = unique(df.N)
Δ = unique(df.Δ)
K = unique(df.K)
n_levels = unique(df.n)

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

ϵ_2_2_save = Vector(range(0.5,10,length=20)) #save these screenshots

#plotting stuff
linea_negra(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)
linea_fija(ϵ_1,ϵ_2) = 4*ϵ_2
linea_primera(ϵ_1,ϵ_2) = 4*ϵ_1*sqrt(ϵ_2)
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

# same plot but t is now ϵ_1 instead of ϵ_2
# @gif for ϵ_1 in ϵ_1_array
#     ΔE_plots = zeros(length(ϵ_2_array), length(n_levels))
#     for n in n_levels
#         n = Int(n)
#         ΔE_plots[:,n] = filter(row ->  row.ϵ_1 == ϵ_1 && row.n == n, df).ΔE_n
#     end
#     scatter(ϵ_2_array, ΔE_plots,markersize=1,markerstrokewidth=0,title="H_eff, N="*string(N[1])*", n_levels="*string(maximum(n_levels)),xlabel=L"\epsilon_2/K",ylabel=L"(E_n-E_0)/K",legend=false)#,ylim=(-5, 150))
# end