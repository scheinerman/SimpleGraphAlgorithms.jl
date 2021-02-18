export hom, hom_check

"""
`hom_check(G,H,d)` checks if the dictionary `d` is a graph homomorphism
from `G` to `H`.
"""
function hom_check(G::SimpleGraph, H::SimpleGraph, d::Dict)::Bool
    for e in G.E
        v, w = e
        x = d[v]
        y = d[w]
        if !H[x, y]
            return false
        end
    end
    return true
end

"""
`hom(G,H)` finds a graph homomorphism from `G` to `H`.
"""
function hom(G::SimpleGraph{S}, H::SimpleGraph{T}) where {S,T}


    m = Model(get_solver())

    @variable(m, P[G.V, H.V], Bin)

    for v in G.V
        @constraint(m, sum(P[v, x] for x in H.V) == 1)
    end

    for e in G.E
        u, v = e
        for x in H.V
            for y in H.V
                if !H[x, y]
                    @constraint(m, P[u, x] + P[v, y] <= 1)
                end
            end
        end
    end


    optimize!(m)
    status = Int(termination_status(m))

    if status != 1
        error("No homomorphism")
    end

    PP = Int.(value.(P))

    result = Dict{S,T}()

    for v in G.V
        for x in H.V
            if PP[v, x] > 0
                result[v] = x
            end
        end
    end

    return result

end
