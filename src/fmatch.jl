export fractional_matching

"""
`fractional_matching(G)` returns a maximum fractional matching in the form of 
a dictionary mapping edges of `G` to rational values. These values are 
in the set `{0, 1/2, 1}`.
"""
function fractional_matching(G::SimpleGraph{T}) where {T}

    if cache_check(G, :fractional_matching)
        return cache_recall(G, :fractional_matching)
    end

    m = Model(get_solver())
    @variable(m, w[G.E], Int)

    for e in G.E
        @constraint(m, 0 <= w[e] <= 2)
    end

    for v in G.V
        star = [e for e in G.E if e[1] == v || e[2] == v]
        @constraint(m, sum(w[e] for e in star) <= 2)
    end

    @objective(m, Max, sum(w[e] for e in G.E))
    optimize!(m)
    wts = Int.(value.(w))

    d = Dict{Tuple{T,T},Rational{Int}}()
    for e in G.E
        d[e] = wts[e] // 2
    end

    cache_save(G, :fractional_matching, d)
    return d
end
