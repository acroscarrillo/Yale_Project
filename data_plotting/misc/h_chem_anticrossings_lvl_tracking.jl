using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs


# Quick plot to see levels
df = DataFrame(CSV.File("data/h_chem_anticrossings_lvl_tracking.csv"))

k_1_array = unique(df.k_1)
k_2_array = unique(df.k_2)
N_cross = maximum(unique(df.cross_n))

l = @layout [ttl; a{.5h}; a b c; d e f; g h i ] 

plots = []

combined_mat = zeros((length(k_2_array),length(k_1_array)))
for n=1:N_cross
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(k_2_array),length(k_1_array)))
    combined_mat += temp_mat.^-1
end
heatmap_combined = heatmap(k_1_array,k_2_array,combined_mat,xlab=L"k_1/k_4",ylab=L"k_2/k_4",title="All plots bellow summed (plus some more)",c = :viridis,clim=(-10,3),titlefontsize=8,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),height=400)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])
plot!([9,13,16.5,19.6,22.5], seriestype="hline",color="black",label=false,style=:dash)
push!(plots,heatmap_combined)

for n=2:7 # ignore last for the sake of plotting
    nn = n+1
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(k_2_array),length(k_1_array))).^-1
    heatmap_tmp = heatmap(k_1_array,k_2_array,temp_mat,xlab=L"k_1/k_4",ylab=L"k_2/k_4",title="levels $n and $nn",c = :viridis,clim=(-10,3),titlefontsize=7,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),xticks=([0,14],["0","14"]),yticks=([0,25],["0","25"]),colorbar=false)
    push!(plots,heatmap_tmp)
end


ttl = "Logarithmic inverse of anti-crossings spacing, "*L"1/\log|\delta_{n,n+1}|"*" of "*L"H_{chem}=p^2 + k_1 x - k_2 x^2 + k_4 x^4."
global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -55Plots.px, titlefontsize=8)

plot(global_title,plots..., layout = l)


savefig("figs/H_chem_anticrossings_oscillations.png")
savefig("figs/H_chem_anticrossings_oscillations.pdf")

fumada(A,x) = exp(-1/(1-A*exp(-x^2))) #forma de los anti-crossings