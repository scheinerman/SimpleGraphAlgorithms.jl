module SimpleGraphAlgorithms
using SimpleGraphs
using MathProgBase
using JuMP

# using Cbc
using Gurobi

# _SOLVER() = CbcSolver()
_SOLVER() = GurobiSolver(OutputFlag=0)

export max_indep_set, max_clique, max_matching, min_dom_set
export min_vertex_cover, min_edge_cover

import SimpleGraphs.cache_save

"""
`max_indep_set(G)` returns a maximum size independent set of a
`SimpleGraph`.
"""
function max_indep_set(G::SimpleGraph)
    if NV(G)==0
      T = vertex_type(G)
      return Set{T}()
    end
    if cache_check(G,:max_indep_set)
      return cache_recall(G,:max_indep_set)
    end
    n = NV(G)
    A = incidence(G,false)'
    c = ones(n)

    X = mixintprog(-c,A,'<',1,:Int,0,1,_SOLVER())

    indices = find(round.(Int,X.sol))
    VV = vlist(G)
    result = Set(VV[indices])
    cache_save(G,:max_indep_set,result)
    return result
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
function max_clique(G::SimpleGraph)
  if cache_check(G,:max_clique)
    return cache_recall(G,:max_clique)
  end
  result = max_indep_set(G')
  cache_save(G,:max_clique,result)
  return result
end

"""
`max_matching(G)` returns a maximum matching of a `SimpleGraph`.
"""
function max_matching(G::SimpleGraph)
    if cache_check(G,:max_matching)
      return cache_recall(G,:max_matching)
    end
    m = NE(G)
    if m == 0
      T = vertex_type(G)
      S = Tuple{T,T}
      return Set{S}()
    end
    A = incidence(G,false)
    c = ones(m)

    X = mixintprog(-c,A,'<',1,:Int,0,1,_SOLVER())

    indices = find(round.(Int,X.sol))
    EE = elist(G)
    result = Set(EE[indices])
    cache_save(G,:max_matching,result)
    return result
end

"""
`min_edge_cover(G)` finds a smallest subset `F` of the edges such that
every vertex of `G` is the end point of at least one edge in `F`. An
error is raised if `G` has a vertex of degree 0.
"""
function min_edge_cover(G::SimpleGraph)
    if cache_check(G,:min_edge_cover)
      return cache_recall(G,:min_edge_cover)
    end

    if NV(G) == 0
      T = vertex_type(G)
      S = Tuple{T,T}
      return Set{S}()
    end

    if minimum(deg(G))==0
        error("Graph has an isolated vertex; no minimum edge cover exists.")
    end
    m = NE(G)
    M = incidence(G,false)
    c = ones(m)

    X = mixintprog(c,M,'>',1,:Int,0,1,_SOLVER())

    indices = find(round.(Int,X.sol))
    EE = elist(G)
    result = Set(EE[indices])
    cache_save(G,:min_edge_cover,result)
    return result
end

"""
`min_dom_set(G)` returns a smallest dominating set of a
`SimpleGraph`. That is, a smallest set `S` with the property that
every vertex of `G` either is in `S` or is adjacent to a vertex of
`S`.
"""
function min_dom_set(G::SimpleGraph)
    if cache_check(G,:min_dom_set)
      return cache_recall(G,:min_dom_set)
    end
    n = NV(G)
    if n==0
      T = vertex_type(G)
      return Set{T}()
    end

    A = adjacency(G)+eye(Int,n)
    c = ones(n)

    X = mixintprog(c,A,'>',1,:Int,0,1,_SOLVER())

    indices = find(round.(Int,X.sol))
    VV = vlist(G)
    result = Set(VV[indices])
    cache_save(G,:min_dom_set,result)
    return result
end

include("iso.jl")
include("kcolor.jl")
include("chrome_poly.jl")
include("mad.jl")

end # of module SimpleGraphAlgorithms
