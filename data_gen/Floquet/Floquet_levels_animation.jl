df_h_eff = DataFrame(CSV.File("data/h_eff_crossings_epsilon_1.csv"))

# Define parameter space
N_h_eff = unique(df_h_eff.N)
Δ_h_eff = unique(df_h_eff.Δ)
K_h_eff = unique(df_h_eff.K)
n_levels_h_eff = unique(df_h_eff.n)

ϵ_1_array_h_eff = unique(df_h_eff.ϵ_1)
ϵ_2_array_h_eff = unique(df_h_eff.ϵ_2)


df_floquet = DataFrame(CSV.File("data/floquet_resonances.csv"))

lemniscate_surface(ϵ_2,ϵ_1_array) = 1 .* ( (ϵ_2)^2 .+ 2 .* ϵ_1_array .* √(ϵ_2) ) #(0.15/240)
fudge_fac = 3
linea_madre(ϵ_2,ϵ_1_array) = 4 .* ϵ_1_array .* √(ϵ_2)
linea_padre(ϵ_2,ϵ_1_array) = 4 * (ϵ_2)

anim = @animate for i in ProgressBar(range(1,length(unique(df_floquet.ϵ_2))))
    ϵ_2 = unique(df_floquet.ϵ_2)[i]
    df_temp = filter(row -> row.ϵ_2 == ϵ_2, df_floquet)
    ϵ_1_array_2_plot = (df_temp.ϵ_1 .* ω_0_exp)./ (K*ω_0_exp)
    ϵ_2_2_plot = ω_0_exp*ϵ_2/(K*ω_0_exp)
    y_max = lemniscate_surface(ϵ_2_2_plot,ϵ_1_array[end])
    scatter(ϵ_1_array_2_plot, df_temp.ϵ_n .* ω_0_exp ./ (K*ω_0_exp),ms=1,ylims=(0,y_max),legend=false,title=L"\epsilon_2/K="*string((ω_0_exp*ϵ_2)/(K*ω_0_exp)),markerstrokewidth=0,xlab=L"\epsilon_1/K",ylab=L"(\epsilon_n-\epsilon_0)/K",c=:black,dpi=300)
    plot!(ϵ_1_array_2_plot,lemniscate_surface(ϵ_2_2_plot,ϵ_1_array_2_plot),c=:red,lw=1.5)
    plot!(ϵ_1_array_2_plot,linea_madre(ϵ_2_2_plot,ϵ_1_array_2_plot),c=:red,lw=1.5)
    plot!(ϵ_1_array_2_plot,fudge_fac.*linea_madre(ϵ_2_2_plot,ϵ_1_array_2_plot),c=:red,lw=1.5)
    plot!(ϵ_1_array_2_plot,fudge_fac.*linea_padre.(ϵ_2_2_plot,ϵ_1_array_2_plot),c=:red,lw=1.5)


end
gif(anim, "anim_fps15.gif", fps = 10)


#plotting stuff
linea_negra(ϵ_1,ϵ_2) = ϵ_2^2 + 2*ϵ_1*sqrt(ϵ_2)

@gif for ϵ_2 in ProgressBar(ϵ_2_array_h_eff)
    ΔE_plots = zeros(length(ϵ_1_array_h_eff), length(n_levels_h_eff))
    for n in n_levels_h_eff
        n = Int(n)
        ΔE_plots[:,n] = filter(row ->  row.ϵ_2 == ϵ_2 && row.n == n, df_h_eff).ΔE_n
    end
    scatter(ϵ_1_array_h_eff, ΔE_plots,markersize=1,markerstrokewidth=0,title="H_eff, N="*string(N_h_eff[1])*", n_levels="*string(maximum(n_levels_h_eff))*", "*L"\epsilon_2/K="*string(ϵ_2),xlabel=L"\epsilon_1/K",ylabel=L"(E_n-E_0)/K",legend=false,ylim=(-5, 150),dpi=300)
    lemniscate_surf = linea_negra.(ϵ_1_array_h_eff,ϵ_2)
    plot!(ϵ_1_array_h_eff, lemniscate_surf,linestyle=:dash,color=:black,dpi=300,ylims=(0,lemniscate_surf[end]*1.15))
    plot!(ϵ_1_array_h_eff, lemniscate_surf.*1.15,linestyle=:dash,color=:black,dpi=300,ylims=(0,lemniscate_surf[end]*1.15))
end
