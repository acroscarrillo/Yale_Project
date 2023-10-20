include("../../src/src.jl")

using Plots
using ProgressBars

#########################
# Lemniscate parameters #
#########################
x_array = Vector(range(-6,6,length=100)) #it's a bit choppy but oh well
p_array = Vector(range(-6,6,length=100))

ϵ_1, ϵ_2 = 8, 8

l = @layout [potential lemni ] 


###########################################################
# Plot potential, its derivative and its barrier position #
###########################################################
H_p0_cut(x) = H_cl(x,0,ϵ_1,ϵ_2)
potential = plot(x_array, H_p0_cut.(x_array),label=L"V(x)",ylab=L"E",xlab=L"x",title="Double well at "*L"\epsilon_1="*"$ϵ_1, "*L"\epsilon_2="*"$ϵ_2",titlefontsize=10)

H_p0_cut_derivative(x) = gradient(H_p0_cut, x)[1]
plot!(x_array, H_p0_cut_derivative.(x_array),label=L"dV(x)/dx")

x_barrier = find_zeros(H_p0_cut_derivative,x_array[1],x_array[end])[2] #should be the middle value
E_barrier = H_p0_cut(x_barrier)
plot!([E_barrier],seriestype="hline",color="black",label=false,style=:dash)
plot!([x_barrier],seriestype="vline",color="black",label=false,style=:dash)
plot!([0],seriestype="vline",color="black",label=false,style=:solid)
plot!([0],seriestype="hline",color="black",label=false,style=:solid)


####################
# Plot lemniscates #
####################
lem_x_points_array = []
lem_p_points_array = []
for x in ProgressBar(x_array)
    lemniscate(p) = H_cl(x,p,ϵ_1,ϵ_2) - E_barrier
    p_roots = find_zeros(lemniscate,x_array[1],x_array[end])
    for p in p_roots
        push!(lem_x_points_array,x)
        push!(lem_p_points_array,p)
    end
end

H0_x_points_array = []
H0_p_points_array = []
for x in ProgressBar(x_array)
    lemniscate(p) = H_cl(x,p,ϵ_1,ϵ_2) - 0
    p_roots = find_zeros(lemniscate,x_array[1],x_array[end])
    for p in p_roots
        push!(H0_x_points_array,x)
        push!(H0_p_points_array,p)
    end
end

lemni = scatter(lem_x_points_array,lem_p_points_array,markersize=1,xlab=L"x",ylab=L"p",legend=false,title="Lemniscate at barrier top and "*L"H(x,p)=0",titlefontsize=10,color=:blue,markerstrokewidth=0)
scatter!(H0_x_points_array,H0_p_points_array,markersize=1,color=:red,markerstrokewidth=0)

# send it!
plot(potential, lemni, layout = l,size=(750,400),dpi=650)

# savefig("figs/important_figs/H_eff/lemniscate_plots.png")