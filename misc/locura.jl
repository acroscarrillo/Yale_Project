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
