module SimpleGraphAlgorithms
using SimpleGraphs
using MathProgBase
using JuMP

export max_indep_set, max_clique, max_matching, min_dom_set
export min_vertex_cover, min_edge_cover

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
`min_vertex_cover(G)` returns a smallest vertex cover `S` of `G`.
This is a smallest possible set of vertices such that every edge 
of `G` has one or both end points in `S`.
"""
min_vertex_cover(G) = setdiff(G.V, max_indep_set(G))

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
`min_edge_cover(G)` finds a smallest subset `F` of the edges such that
every vertex of `G` is the end point of at least one edge in `F`. An
error is raised if `G` has a vertex of degree 0.
"""
function min_edge_cover(G::SimpleGraph)
    if minimum(deg(G))==0
        error("Graph has an isolated vertex; no minimum edge cover exists.")
    end
    m = NE(G)
    M = incidence(G,false)
    c = ones(m)

    X = mixintprog(c,M,'>',1,:Int,0,1)

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
include("kcolor.jl")
include("chrome_poly.jl")

end # of module SimpleGraphAlgorithms
