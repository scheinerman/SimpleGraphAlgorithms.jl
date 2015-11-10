export iso, iso_matrix, iso_check

"""
`iso_matrix(G,H)` returns a permutation matrix `P` such that
`A*P==P*B` where `A` is the adjacency matrix of `G` and `B` is the
adjacency matrix of `H`.  If the graphs are not isomorphic, an
error is raised.
"""

function iso_matrix(G::SimpleGraph, H::SimpleGraph)
    n = NV(G)

    err_msg = "The graphs are not isomorphic"

    # quickly rule out some easily detected nonisomorphic graphs
    if NV(H) != n || NE(G) != NE(H) || deg(G) != deg(H)
        error(err_msg)
    end

    m = Model()
    @defVar(m,P[1:n,1:n],Bin)
    A = adjacency(G)
    B = adjacency(H)

    for i=1:n
        for k=1:n
            @addConstraint(m,
                           sum{A[i,j]*P[j,k], j=1:n} ==
                           sum{P[i,j]*B[j,k], j=1:n}
                           )
        end
        @addConstraint(m, sum{P[i,j], j=1:n}==1)
        @addConstraint(m, sum{P[j,i], j=1:n}==1)
    end

    status = solve(m, suppress_warnings=true)

    if status != :Optimal
        error(err_msg)
    end

    return round(Int,getValue(P))
end

"""
`iso(G,H)` finds an isomorphism from `G` to `H` (or throws an error if
the graphs are not isomorphic). Returns a dictionary mapping the
vertices of `G` to `H`.
"""
function iso(G::SimpleGraph, H::SimpleGraph)
    err_msg = "The graphs are not isomorphic"

    # quickly rule out some easily detected nonisomorphic graphs
    if NV(G) != NV(H) || NE(G) != NE(H) || deg(G) != deg(H)
        error(err_msg)
    end

    VG = vlist(G)
    VH = vlist(H)
    n = NV(G)

    MOD = Model()

    @defVar(MOD, x[VG,VH],Bin)

    for v in VG
        @addConstraint(MOD, sum{x[v,VH[k]], k=1:n}==1)
    end

    for w in VH
        @addConstraint(MOD, sum{x[VG[k],w], k=1:n}==1)
    end

    for v in VG
        for w in VH
            @addConstraint(MOD,
                           sum{ has(G,v,VG[k])*x[VG[k],w], k=1:n } ==
                           sum{ x[v,VH[k]]*has(H,VH[k],w), k=1:n }
                           )
        end
    end

    status = solve(MOD, suppress_warnings=true)

    if status != :Optimal
        error(err_msg)
    end

    X = getValue(x)

    result = Dict{vertex_type(G), vertex_type(H)}()

    for v in VG
        for w in VH
            if X[v,w] > 0
                result[v] = w
            end
        end
    end
    return result
end


"""
`iso_check(G,H,d)` checks if `d` is an isomorphism from `G` to `H`.
"""
function iso_check(G::SimpleGraph, H::SimpleGraph, d::Dict)
    n = NV(G)

    # standard quick check
    if n != NV(H) || NE(G) != NE(H) || deg(G) != deg(H) || length(d) != n
        return false
    end

    # check keys and values in d match VG and VH, respectively
    for k in keys(d)
        if !has(G,k)
            return false
        end
    end

    for v in values(d)
        if !has(H,v)
            return false
        end
    end

    # Check edges match. Only need to go one way since we determined
    # that the two graphs have the same number of edges.
    EE = elist(G)
    for e in EE
        v = e[1]
        w = e[2]
        if !has(H,d[v],d[w])
            return false
        end
    end

    return true
end
