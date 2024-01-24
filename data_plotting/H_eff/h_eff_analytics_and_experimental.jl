include("../../src/src.jl") # import charge basis stuff

using DataFrames # this is like pandas
using CSV 
using ProgressBars
using LaTeXStrings # latex support for sexy figs
using NPZ
using EasyFit
using Roots
using Flux
using JLD


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

ϵ_1_array = Vector( range(0,13.69,length=300) )
ϵ_2_array = Vector( range(2.5,12.5,length=300) )

# Define crosses and stuff
hlines_labels = [L"B'",L"A'",L"B",L"A"]
h_lines_pos = [5.9,7.6,9.7,11.25]
cross_loc = [(2.45,5.9),(5.5,7.6),(3.15,9.7),(9.3,9.7),(7.3,5.9),(6.75,11.25),(13.4,11.25) ]

# Experiment
exp_data = npzread("data/experimental/DecayData.npz")

K_exp = 0.505 # MHz or 505KHz

heatmap_mat_exp = exp_data["arr_2"] # μs, but irrelant for now
ϵ_1_array_exp = exp_data["arr_0"][1,:] # MHz
ϵ_2_array_exp = exp_data["arr_1"][:,1] # MHz

# remove the upward trend
for (j,ϵ2) in enumerate(ϵ_2_array_exp)
    heatmap_mat_exp[j,:] .= movingaverage(heatmap_mat_exp[j,:],1).x .- movingaverage(heatmap_mat_exp[j,:],50).x
end

ϵ_1_array_exp = ϵ_1_array_exp./K_exp # Units of Kerr
ϵ_2_array_exp = ϵ_2_array_exp./K_exp # Units of Kerr


# without trend formatting
# heatmap(ϵ_1_array_exp,ϵ_2_array_exp,heatmap_mat_exp,xlab=L"\epsilon_x/K",ylab=L"\epsilon_2/K",title="",titlefontsize=12,xtickfontsize=12,ytickfontsize=12,guidefont=font(16),height=400,colorbar=true,dpi=600, size=(720,400),cbar_title="\nLifetime "*L"T_{coh}(\mu s)",widen=false,tickdirection=:out,right_margin = 8Plots.mm,colorbar_tickfontsize=16,left_margin = 3Plots.mm,bottom_margin = 3Plots.mm)

# with trend formatting
heatmap(ϵ_1_array_exp[1:end-20],ϵ_2_array_exp,heatmap_mat_exp[1:end,1:end-20],xlab=L"\epsilon_x/K",ylab=L"\epsilon_2/K",title="",titlefontsize=12,xtickfontsize=14,ytickfontsize=14,guidefont=font(16),height=400,colorbar=false,dpi=600, size=(720,400),clim=(-125,80),cbar_title="Lifetime",widen=false,tickdirection=:out,right_margin = 8Plots.mm,colorbar_tickfontsize=16,left_margin = 3Plots.mm,bottom_margin = 3Plots.mm)

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

# save("data/well_area_matrices.jld", "left", left_well_area_matrix, "right", right_well_area_matrix)

left_well_area_matrix = load("data/well_area_matrices.jld", "left")
right_well_area_matrix = load("data/well_area_matrices.jld", "right")

# contour(ϵ_1_array,ϵ_2_array,left_well_area_matrix')
# contour!(ϵ_1_array,ϵ_2_array,right_well_area_matrix')

###################
# Analytics plots #
###################
# plot(ϵ_1_array,(2*ϵ_1_array).^(2/3),ylim=(2.5,12),linestyle=:solid,legend=false,lw=2,xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",dpi=600)

plot!(ϵ_1_array,(ϵ_1_array/1).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,dpi=600,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/2).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/3).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/4).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)


plot!(ϵ_1_array,(ϵ_1_array/1).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,dpi=600,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/2).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/3).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)
plot!(ϵ_1_array,(ϵ_1_array/4).^2,ylim=(2.5,12),legend=false,linewidth=2,alpha=0.5,linestyle=:solid,c=:limegreen)


###################
# Equi-Area plots #
###################
# plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0, 0.05),ylim=(2.5,12),lw=1,alpha=0.25,c=:white)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*0.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*1.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*1.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*2.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*2.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*3.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*3.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*4.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*4.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

# plot!(contour_line(left_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*5.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)
plot!(contour_line(right_well_area_matrix,ϵ_1_array,ϵ_2_array,2*π*5.5, 0.05),ylim=(2.5,12),lw=1.5,alpha=0.5,c=:white)

# plot!(ϵ_1_array[1:160],(-4.7*(ϵ_1_array[1:160] .- 7.5)).^(1/2) .+ 5,c=:white)

savefig("figs/important_figs/H_eff/analytics_and_experimental.png")
savefig("figs/important_figs/H_eff/analytics_and_experimental.pdf")
savefig("figs/important_figs/H_eff/analytics_and_experimental.svg")


