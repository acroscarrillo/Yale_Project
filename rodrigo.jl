using LinearAlgebra
using BenchmarkTools


function C_pwr(N,k)
    temp = im*zeros((N,N))
    omega = exp(im * 2 *k *π/N)
    for i=1:N
        temp[i,i] = omega^i
    end
    return temp
end

function S_pwr(N,k)
    temp = zeros((N,N))
    k = k-1
    for i=1:(N-1)
        temp[i,mod(i+k,N)+1] = 1
    end
    temp[N,mod(k,N)+1]=1
    return temp
end


function Tkl(k,l,N)
    C_pwr_k = C_pwr(N,k)
    S_pwr_l = S_pwr(N,l)
    return exp(-im * π / N * k * l) * C_pwr_k * S_pwr_l
end


function tnm_M(n1,n2,m1,m2,N)
    Cnm = Tkl(n1,n2,N)*Tkl(m1,m2,N) - Tkl(m1,m2,N)*Tkl(n1,n2,N)
    return tr(Cnm*Cnm')
end

function f(N)
    temp = zeros((N,N,N,N))
    for n1=1:N
        C_pwr_n1 = C_pwr(N,n1)
        for n2=1:N
            S_pwr_n2 = S_pwr(N,n2)
            Tkl_n1n2 = exp(-im * π / N * n1 * n2) * C_pwr_n1 * S_pwr_n2
            for m1=1:N
                C_pwr_m1 = C_pwr(N,m1)
                for m2=1:N
                    S_pwr_m2 = S_pwr(N,m2)
                    Tkl_m1m2 = exp(-im * π / N * m1 * m2) * C_pwr_m1 * S_pwr_m2
                    Cnm = Tkl_n1n2*Tkl_m1m2 - Tkl_m1m2*Tkl_n1n2
                    temp[n1,n2,m1,m2] = tr(Cnm*Cnm')
                end
            end
        end
    end
end

function S_to_the(N,k)
    temp = zeros((N,N))
    for i=1:(N-1)
        temp[i*K,i+1] = 1
        temp[N-k+1,k]
    end
    temp[N,1]=1
    return temp^k
end

