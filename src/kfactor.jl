export kfactor

"""
`kfactor(G,k=1)` returns a `k`-factor of `G`. This is a set of edges
with the property that every vertex of `G` is incident with exactly `k`
edges in the set. An error is thrown if no such set exists.
"""
function kfactor(G::SimpleGraph, k::Int = 1)
    @assert k > 0 "The parameter k must be positive"

    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    MOD = Model(get_solver())

    @variable(MOD, x[EE], Bin)

    for v in VV
        star = [get_edge(G, v, w) for w in G[v]]
        @constraint(MOD, sum(x[e] for e in star) == k)
    end

    optimize!(MOD)

    status = Int(termination_status(MOD))

    if status != 1
        error("This graph does not have a $k-factor")
    end

    X = value.(x)
    result = Set(e for e in EE if X[e] > 0.5)
    return result
end
