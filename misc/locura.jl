using LinearAlgebra # to diagonalise stuff
using Plots

N = 70
ω_0 = 1
g3 = 1*(0.00075)*1^-1 .*ω_0
g4 = 1*(1.27*10^(-7))*1^-1 .*ω_0

H_test_1 = ω_0*a(N)'*a(N) + (g3/3)*(a(N)'+a(N))^3 + (g4/4)*(a(N)'+a(N))^4

lamb_test_1, _ = eigen(H_test_1)

display(lamb_test_1[3])
display(lamb_test_1[2])
display(lamb_test_1[1])

alpha_test_1 = lamb_test_1[3]+lamb_test_1[1]-2*lamb_test_1[2]
display("-α/2ω_0 = "*string(-alpha_test_1/2))

lamb_test_1 = lamb_test_1 .- Vector(0:N-1)*(lamb_test_1[2]-lamb_test_1[1])
    

###############
K_test_2 = 10*(g3)^2/(3*ω_0) - 3*(g4)/2
display("K/ω_0 = "*string(K_test_2))

###############

H_test_2 = (ω_0-2*K_test_2)*a(N)'*a(N) - K_test_2*a(N)'*a(N)'*a(N)*a(N)
lamb_test_2, _ = eigen(H_test_2)
# display(lamb_test_2) 
lamb_test_2 = lamb_test_2 .- Vector(0:N-1)*(lamb_test_2[2]-lamb_test_2[1])
# display(lamb_test_2)

# scatter(lamb_test_1[1:20] .-lamb_test_2[1:20])
scatter(lamb_test_2[1:20]./(-alpha_test_1/2),linestyle=:dash)
scatter!(lamb_test_1[1:20]./K_test_2,linestyle=:dash,ylim=(-100,5))


# plot(lamb_test_1-lamb_test_2, ylim=(-5e-4,5e-4),xlim=(0,25))






function H_cos(x,p,E_C,E_J,α,m,ϕ_ext)
    H_kin = 4*E_C*(p)^2
    return H_kin + potential(x,E_J,α,m,ϕ_ext)
end

function potential(x,E_J,α,m,ϕ_ext)
    N = size(x)[1]
    H_cos1 = - α*E_J*cos( x ) 
    H_cos2 = - m*E_J*cos( I(N)*2*π*ϕ_ext/m - x/m )
    return H_cos1 + H_cos2
end

function potential_classic(x,E_J,α,m,ϕ_ext)
    H_cos1 = - α*E_J*cos( x ) 
    H_cos2 = - m*E_J*cos( 2*π*ϕ_ext/m - x/m )
    return H_cos1 + H_cos2
end




function H_charge_basis(N,m,α,ϕ_ext,E_C,E_J)
    N_op = Diagonal(Vector(-N:1:N-1))

    step = Matrix(Tridiagonal(repeat([0],2*N-1),repeat([0],2*N),repeat([1],2*N-1)))
    step[1,end] = 0
    step[end,1] = 1
    display(step)
    cos_m = (exp(-im*ϕ_ext*2*π)*step)^(1/m)
    cos_m = (cos_m + cos_m')/2

    cos_1 = Tridiagonal(repeat([1],2*N-1),repeat([0],2*N),repeat([1],2*N-1))/2

    return 4*E_C*N_op^2 - E_J*(α*cos_1 + m*cos_m)
end




function taylor(f,order_n,x_0)
    coeffs = zeros(order_n)
    coeffs[1] = f(x_0)/factorial(1)
    g_n = x -> f(x)
    for n=2:order_n
        g_n = x -> gradient(g_n,x)[1]
        coeffs[n] = g_n(x_0)/factorial(n)
    end
    return coeffs
end

f_0(x) = potential_classic(x,5907,0.1,3,0.33)
f_1(x) = gradient(f_0,x)[1]
f_2(x) = gradient(f_1,x)[1]
f_3(x) = gradient(f_2,x)[1]
f_4(x) = gradient(f_3,x)[1]
f_5(x) = gradient(f_4,x)[1]
f_6(x) = gradient(f_5,x)[1]

x_0 = 1.7794902664439491
coeff_0 = f_0(x_0)/factorial(0)
coeff_1 = f_1(x_0)/factorial(1)
coeff_2 = f_2(x_0)/factorial(2)
coeff_3 = f_3(x_0)/factorial(3)
coeff_4 = f_4(x_0)/factorial(4)
coeff_5 = f_5(x_0)/factorial(5)
coeff_6 = f_5(x_0)/factorial(6)

taylor(x) = coeff_0 + coeff_1*(x-x_0) + coeff_2*(x-x_0)^2+ coeff_3*(x-x_0)^3+ coeff_4*(x-x_0)^4+ coeff_5*(x-x_0)^5 + coeff_6*(x-x_0)^6