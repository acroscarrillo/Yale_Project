using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs


function find_barrier_height(ϵ_1,ϵ_2)
    H_p0_cut(x) = H_cl(x,0,ϵ_1,ϵ_2)
    H_p0_cut_derivative(x) = gradient(H_p0_cut, x)[1]
    all_zeros = find_zeros(H_p0_cut_derivative,x_array[1],x_array[end])
    if length(all_zeros)==3
        x_barrier = sort(all_zeros)[2]
        E_barrier = H_p0_cut(x_barrier)
        return E_barrier, x_barrier
    else 
        return nothing, nothing
    end
end


function left_well_area_MC(ϵ_1,ϵ_2,E_barrier,x_barrier,x_array,p_array)
    n=0
    for x in x_array
        for p in p_array
            if (H_cl(x,p,ϵ_1,ϵ_2)-E_barrier <= 0) && (x <= x_barrier) # left well
                n += 1
            end
        end
    end
    total_points = length(x_array) * length(p_array) 
    total_area = (x_array[end]-x_array[1])*(p_array[end]-p_array[1])
    return total_area*(n/total_points)
end


function right_well_area_MC(ϵ_1,ϵ_2,E_barrier,x_barrier,x_array,p_array)
    n=0
    for x in x_array
        for p in p_array
            if (H_cl(x,p,ϵ_1,ϵ_2)-E_barrier <= 0) && (x >= x_barrier) # left well
                n += 1
            end
        end
    end
    total_points = length(x_array) * length(p_array) 
    total_area = (x_array[end]-x_array[1])*(p_array[end]-p_array[1])
    return total_area*(n/total_points)
end


function contour_line(mat,ϵ_1_array,ϵ_2_array,val, tol)
    indxs =  findall(val+tol .>= mat .>= val-tol)
    x_2_plot = []
    y_2_plot = []
    for ind in indxs
        i,j = ind[1], ind[2]
        push!(x_2_plot,ϵ_1_array[i])
        push!(y_2_plot,ϵ_2_array[j])
    end
    return x_2_plot, y_2_plot
end






df = DataFrame(CSV.File("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv"))

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

################
# Reso Heatmap #
################
combined_mat = zeros((length(ϵ_2_array),length(ϵ_1_array)))
for n=2:7 # ignore last for the sake of plotting
    df_2_plot = filter(row ->   row.cross_n == n, df)
    temp_mat = reshape(log.(df_2_plot.Δnn),(length(ϵ_2_array),length(ϵ_1_array)))
    combined_mat += temp_mat.^-1
end
heatmap_combined = heatmap(ϵ_1_array,ϵ_2_array,combined_mat,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title="",clim=(-5,5),titlefontsize=8,xtickfontsize=8,ytickfontsize=8,guidefont=font(10),height=400,c = :oleron100,dpi=600)#,xticks=([0,12],["0","12"]),yticks=([0,10],["0","10"])

# cross_loc = [(2.45,5.9),(5.5,7.6),(3.15,9.7),(3.9,3.8),(9.3,9.7),(10.8,7.2),(7.1,5.5),(6.75,11.25),(13.4,11.25) ]
# plot!(cross_loc,marker=:o,color=:red,markersize=6,legend=false)

################
# Area Heatmap #
################
# x_array = Vector(range(-10,10,length=1000)) #it's a bit choppy but oh well
# p_array = Vector(range(-6,6,length=1000))

# left_well_area_matrix = zeros((length(ϵ_2_array),length(ϵ_1_array)))
# right_well_area_matrix = zeros((length(ϵ_2_array),length(ϵ_1_array)))
# for (j,ϵ_1) in ProgressBar(enumerate(ϵ_1_array))
#     for (k,ϵ_2) in enumerate(ϵ_2_array)

#         E_barrier, x_barrier = find_barrier_height(ϵ_1,ϵ_2)

#         if E_barrier !== nothing
#             left_well_area_matrix[j,k] = left_well_area_MC(ϵ_1,ϵ_2,E_barrier,x_barrier,x_array,p_array)
#             right_well_area_matrix[j,k] = right_well_area_MC(ϵ_1,ϵ_2,E_barrier,x_barrier,x_array,p_array)
#         else
#             left_well_area_matrix[j,k] = -0.5
#             right_well_area_matrix[j,k] = -0.5
#         end

#     end
# end

# contour(ϵ_1_array,ϵ_2_array,left_well_area_matrix')
# contour!(ϵ_1_array,ϵ_2_array,right_well_area_matrix')

###################
# Analytics plots #
###################
# plot!(ϵ_1_array,(2*ϵ_1_array).^(2/3),ylim=(2.5,12.5),linestyle=:dash,legend=false,lw=2,c=:black)

plot!(ϵ_1_array,(ϵ_1_array/1).^2,ylim=(2.5,12.5),legend=false,linewidth=3,ls=:solid,alpha=1,c=:black)
plot!(ϵ_1_array,(ϵ_1_array/2).^2,ylim=(2.5,12.5),legend=false,linewidth=3,ls=:solid,alpha=1,c=:black)
plot!(ϵ_1_array,(ϵ_1_array/3).^2,ylim=(2.5,12.5),legend=false,linewidth=3,ls=:solid,alpha=1,c=:black)
plot!(ϵ_1_array,(ϵ_1_array/4).^2,ylim=(2.5,12.5),legend=false,linewidth=3,ls=:solid,alpha=1,c=:black)


###################
# Equi-Area plots #
###################
plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*1.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*1.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*2.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*2.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*3.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*3.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*4.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*4.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*5.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*5.5, 0.05),ylim=(2.5,12),ms=2,alpha=1, ls=:solid,c=:black,lw=2)



# plot!(ϵ_1_array[1:160],(-4.7*(ϵ_1_array[1:160] .- 7.5)).^(1/2) .+ 5)

savefig("figs/important_figs/H_eff/heatmap_and_analytics.png")
