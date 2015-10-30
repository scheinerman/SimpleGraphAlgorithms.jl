module SimpleGraphAlgorithms
using SimpleGraphs
using MathProgBase
using JuMP

export max_indep_set, max_clique, max_matching, min_dom_set

"""
`max_indep_set(G)` returns a maximum size independent set of a
`SimpleGraph`.
"""
function max_indep_set(G::SimpleGraph)
    n = NV(G)
    A = incidence(G,false)'
    c = ones(n)

    X = mixintprog(-c,A,'<',1,:Int,0,1)

    indices = find(round(Int,X.sol))
    VV = vlist(G)
    return Set(VV[indices])
end


"""
`max_clique(G)` returns a maximum size clique of a `SimpleGraph`.
"""
max_clique(G::SimpleGraph) = max_indep_set(G')


"""
`max_matching(G)` returns a maximum matching of a `SimpleGraph`.
"""
function max_matching(G::SimpleGraph)
    m = NE(G)
    A = incidence(G,false)
    c = ones(m)

    X = mixintprog(-c,A,'<',1,:Int,0,1)

    indices = find(round(Int,X.sol))
    EE = elist(G)
    return Set(EE[indices])
end

"""
`min_dom_set(G)` returns a smallest dominating set of a
`SimpleGraph`. That is, a smallest set `S` with the property that
every vertex of `G` either is in `S` or is adjacent to a vertex of
`S`.
"""
function min_dom_set(G::SimpleGraph)
    n = NV(G)
    A = adjacency(G)+eye(Int,n)
    c = ones(n)

    X = mixintprog(c,A,'>',1,:Int,0,1)

    indices = find(round(Int,X.sol))
    VV = vlist(G)
    return Set(VV[indices])
end


include("iso.jl")


end # of module SimpleGraphAlgorithms
