
export k_edge_color, edge_chromatic_number

"""
`k_edge_color(G,k)` returns a proper `k`-edge coloring of `G` or
throws an error if one does not exist.
"""
function k_edge_color(G::SimpleGraph, k::Int)
    err_msg = "This graph is not " * string(k) * "-edge colorable"
    Delta = maximum(deg(G))
    if k<Delta
        error(err_msg)
    end

    err_msg = "This graph is not " * string(k) * " colorable"

    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    VT = vertex_type(G)
    ET = Tuple{VT,VT}

    result = Dict{ET,Int}()

    MOD = Model(solver=GurobiSolver())
    @variable(MOD, x[EE,1:k], Bin)

    # Every edge must have exactly one color
    for ed in EE
        @constraint(MOD, sum(x[ed,i] for i=1:k) == 1)
    end

    # Two edges incident with the same vertex must have different colors
    for i=1:m-1
        for j=i+1:m
            e1 = EE[i]
            e2 = EE[j]
            meet = intersect(Set(e1), Set(e2))
            if length(meet)>0
                for t=1:k
                    @constraint(MOD, x[e1,t]+x[e2,t] <= 1)
                end
            end
        end
    end

    status = solve(MOD,suppress_warnings=true)
    if status != :Optimal
        error(err_msg)
    end

    X = getvalue(x)
    for ee in EE
        for t=1:k
            if X[ee,t] > 0
                result[ee] = t
            end
        end
    end
    return result
end


"""
`edge_chromatic_number(G)` returns the edge chromatic number
of the graph `G`.
"""
function edge_chromatic_number(G::SimpleGraph)::Int
    if cache_check(G,:edge_chromatic_number)
        return cache_recall(G,:edge_chromatic_number)
    end

    d = maximum(deg(G))
    try
        f = k_edge_color(G,d)
        cache_save(G, :edge_chromatic_number, d)
        return d
    end

    cache_save(G, :edge_chromatic_number, d+1)
    return d+1
end
