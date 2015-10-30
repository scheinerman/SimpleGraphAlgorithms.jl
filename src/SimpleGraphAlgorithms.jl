using SimpleGraphs
using MathProgBase

module SimpleGraphAlgorithms

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


end # of module SimpleGraphAlgorithms
