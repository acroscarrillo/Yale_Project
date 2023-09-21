# utility code goes here
function coherent_state(N,α)
    A = a(N)
    vac = zeros(N)
    vac[1] = 1
    coh_state =  exp(α*A-(α*A)')*vac
    return coh_state/norm(coh_state)
end

function IPR(vec)
    return sum( [norm(vec[n])^4 for n=1:length(vec)] )
end