include("../../src/src.jl")

using Plots
using ProgressBars

cross_loc = [(2.45,5.9),(5.5,7.6),(3.15,9.7),(3.9,3.8),(9.3,9.7),(10.8,7.2),(7.1,5.5),(6.75,11.25),(13.4,11.25) ]

# scatter!(cross_loc,marker=:o,color=:red,markersize=6,legend=false)

#########################
# Lemniscate parameters #
#########################
x_array = Vector(range(-8,8,length=100)) #it's a bit choppy but oh well
p_array = Vector(range(-6,6,length=100))

cross_n = 9
ϵ_1, ϵ_2 = cross_loc[cross_n]
# ϵ_1, ϵ_2 = 13.4, (2*13.4)^(2/3)

l = @layout [potential lemni ] 


###########################################################
# Plot potential, its derivative and its barrier position #
###########################################################
H_p0_cut(x) = H_cl(x,0,ϵ_1,ϵ_2)
potential = plot(x_array, H_p0_cut.(x_array),label=L"V(x)",ylab=L"E",xlab=L"x",title="Double well, "*L"H_{cl}(x,p=0)"*", at "*L"\epsilon_1="*"$ϵ_1, "*L"\epsilon_2="*"$ϵ_2",titlefontsize=10,dpi=600)
plot!([0],seriestype="vline",color="black",label=false,style=:solid)
plot!([0],seriestype="hline",color="black",label=false,style=:solid)

H_p0_cut_derivative(x) = gradient(H_p0_cut, x)[1]
# plot!(x_array, H_p0_cut_derivative.(x_array),label=L"dV(x)/dx")
plot!(x_array, H_cl.(x_array,0,0,ϵ_2),label=L"V(x)"*" at "*L"\epsilon_1=0")


x2 = find_zeros(H_p0_cut_derivative,x_array[1],x_array[end])[2] #should be the middle value
E2 = H_p0_cut(x2)
plot!([E2],seriestype="hline",color="red",label=false,style=:dash)
plot!([x2],seriestype="vline",color="red",label=false,style=:dash)
x,E =  round(x2,sigdigits=3), round(E2,sigdigits=3)
scatter!([x2],[E2],marker=:xcross,color="red",label=L"E="*"$E, "*L"x="*"$x")

x1 = find_zeros(H_p0_cut_derivative,x_array[1],x_array[end])[1] #should be the middle value
E1 = H_p0_cut(x1)
plot!([E1],seriestype="hline",color="green",label=false,style=:dash)
plot!([x1],seriestype="vline",color="green",label=false,style=:dash)
x,E =  round(x1,sigdigits=3), round(E1,sigdigits=3)
scatter!([x1],[E1],marker=:xcross,color="green",label=L"E="*"$E, "*L"x="*"$x")

x3 = find_zeros(H_p0_cut_derivative,x_array[1],x_array[end])[3] #should be the middle value
E3 = H_p0_cut(x3)
plot!([E3],seriestype="hline",color="blue",label=false,style=:dash)
plot!([x3],seriestype="vline",color="blue",label=false,style=:dash)
x,E =  round(x3,sigdigits=3), round(E3,sigdigits=3)
scatter!([x3],[E3],marker=:xcross,color="blue",label=L"E="*"$E, "*L"x="*"$x")


# Arrows: Horizontal 
plot!([x1,0],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = L"lenght(L_1)="*string(round(abs(x1-x2),sigdigits=3)), legend=true,arrow=true)
plot!([0,x1],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = false, legend=true,arrow=true)
annotate!(-abs(x1)/2,E3*0.95,L"L_1")

plot!([-sqrt(2*ϵ_2),0],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = L"lenght(L_1)="*string(round(abs(x1-x2),sigdigits=3)), legend=true,arrow=true)
plot!([0,-sqrt(2*ϵ_2)],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = L"lenght(L_1)="*string(round(abs(x1-x2),sigdigits=3)), legend=true,arrow=true)




plot!([0,x3],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = L"lenght(L_2)="*string(round(abs(x3-x2),sigdigits=3)), legend=true,arrow=true)
plot!([x3,0],[E3*1.1,E3*1.1],linewidth=2,c=:hotpink,label = false, legend=true,arrow=true)
annotate!(abs(x3)/2,E3*0.95,L"L_2")

# Arrows: Vertical 
plot!([x2,x2],[E1,E2],linewidth=2,c=:hotpink,label = L"lenght(L_3)="*string(round(abs(E1-E2),sigdigits=3)), legend=true,arrow=true)
plot!([x2,x2],[E2,E1],linewidth=2,c=:hotpink,label = false, legend=true,arrow=true)
annotate!(x2*1.5,-abs(E1-E2)/2+E2,L"L_3")

plot!([x2,x2],[E1,E3],linewidth=2,c=:hotpink,label = L"lenght(L_4)="*string(round(abs(E1-E3),sigdigits=3)), legend=true,arrow=true)
plot!([x2,x2],[E3,E1],linewidth=2,c=:hotpink,label = false, legend=true,arrow=true)
annotate!(x2*1.5,-abs(E3-E2)/2+E2,L"L_4")

# savefig("figs/important_figs/H_eff/Classic_potential_plot.png")
