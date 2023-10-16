using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs
using ShiftedArrays: lag

function df_first_ϵ_1_derivative(df,n,ϵ_2)
    df_temp = filter(row ->  row.ϵ_2 == ϵ_2 && row.cross_n == n, df)
    y_n = lag(df_temp.Δnn, 0)
    y_nn = lag(df_temp.Δnn, 1)
    h = df_temp.ϵ_1[2]-df_temp.ϵ_1[1]
    return (y_nn + y_n)./h^2
end

function df_second_ϵ_1_derivative(df,n,ϵ_2)
    df_temp = filter(row ->  row.ϵ_2 == ϵ_2 && row.cross_n == n, df)
    y_n = lag(df_temp.Δnn, -1)
    y_nn = lag(df_temp.Δnn, 0)
    y_nnn = lag(df_temp.Δnn, 1)
    h = df_temp.ϵ_1[2]-df_temp.ϵ_1[1]
    return (y_nnn - 2*y_nn + y_n)./h^2
end


function df_first_ϵ_2_derivative(df,n,ϵ_1)
    df_temp = filter(row ->  row.ϵ_1 == ϵ_1 && row.cross_n == n, df)
    y_n = lag(df_temp.Δnn, 0)
    y_nn = lag(df_temp.Δnn, 1)
    h = df.ϵ_2[2]-df.ϵ_2[1]
    return (y_nn + y_n)./h
end

function df_second_ϵ_2_derivative(df,n,ϵ_1)
    df_temp = filter(row ->  row.ϵ_1 == ϵ_1 && row.cross_n == n, df)
    y_n = lag(df_temp.Δnn, -1)
    y_nn = lag(df_temp.Δnn, 0)
    y_nnn = lag(df_temp.Δnn, 1)
    h = df_temp.ϵ_2[2]-df_temp.ϵ_2[1]
    return (y_nnn - 2*y_nn + y_n)./h^2
end





# Quick plot to see levels
df = DataFrame(CSV.File("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv"))

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

l = @layout [ttl; a{.5h}; a b c; d e f; g h i ] 

plots = []

second_derivative_mat = zeros((length(ϵ_2_array),length(ϵ_1_array[2:end-1])))
for (j,ϵ_2) in ProgressBar(enumerate(ϵ_2_array))
    second_derivative_mat[j,:] .= df_second_ϵ_1_derivative(df,3,ϵ_2)[2:end-1]
end

first_derivative_mat = zeros((length(ϵ_2_array),length(ϵ_1_array[2:end])))
for (j,ϵ_2) in ProgressBar(enumerate(ϵ_2_array))
    first_derivative_mat[j,:] .= df_first_ϵ_1_derivative(df,4,ϵ_2)[2:end]
end

heatmap(ϵ_1_array[2:end],ϵ_2_array, first_derivative_mat)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])



heatmap(ϵ_1_array[2:end-1],ϵ_2_array,second_derivative_mat)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

heatmap(ϵ_1_array[2:end-1],ϵ_2_array,log.(second_derivative_mat.^-2))#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

heatmap(ϵ_1_array[2:end-1],ϵ_2_array,log.(second_derivative_mat.^2))#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])



# second_derivative_mat_2 = zeros((length(ϵ_2_array[2:end-1]),length(ϵ_1_array)))
# for (j,ϵ_1) in ProgressBar(enumerate(ϵ_1_array))
#     second_derivative_mat_2[j,:] .= df_second_ϵ_2_derivative(df,2,ϵ_1)[2:end-1]
# end

# first_derivative_mat_2 = zeros((length(ϵ_2_array[2:end]),length(ϵ_1_array)))
# for (j,ϵ_1) in ProgressBar(enumerate(ϵ_1_array))
#     first_derivative_mat_2[j,:] .= df_first_ϵ_2_derivative(df,2,ϵ_1)[2:end]
# end

# heatmap(ϵ_1_array,ϵ_2_array[2:end], first_derivative_mat)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

# heatmap(ϵ_1_array,ϵ_2_array[2:end-1], second_derivative_mat)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])



for n=2:7 # ignore last for the sake of plotting
    nn = n+1
    df_2_plot = filter(row ->   row.cross_n == n, df)
    # second_deriv = df_second_ϵ_1_derivative(df,n,ϵ_2)
    temp_mat = reshape(df_2_plot.Δnn,(length(ϵ_2_array),length(ϵ_1_array)))
    combined_mat += temp_mat.^-1
end
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="All plots bellow summed",c = :viridis,clim=(-10,3),titlefontsize=8,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),height=400)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])
plot!([5.9,7.6,9.7], seriestype="hline",color="black",label=false,style=:dash)
push!(plots,heatmap_combined)

for n=2:7 # ignore last for the sake of plotting
    nn = n+1
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array))).^-1
    heatmap_tmp = heatmap(ϵ_1_array,ϵ_2_array,temp_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="levels $n and $nn",c = :viridis,clim=(-10,3),titlefontsize=7,xtickfontsize=6,ytickfontsize=6,guidefont=font(7),xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"]),colorbar=false)
    push!(plots,heatmap_tmp)
end


ttl = "Logarithmic inverse of anti-crossings spacing, "*L"1/\log|\delta_{n,n+1}|"*" of "*L"H_{eff}."
global_title = plot(title = ttl, grid = false, showaxis = false, bottom_margin = -55Plots.px, titlefontsize=8)

plot(global_title,plots..., layout = l)


savefig("figs/important_figs/H_eff/H_eff_anticrossings_oscillations.png")
savefig("figs/important_figs/H_eff/H_eff_anticrossings_oscillations.pdf")

fumada(A,x) = exp(-1/(1-A*exp(-x^2))) #forma de los anti-crossings