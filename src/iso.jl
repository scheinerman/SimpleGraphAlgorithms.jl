export iso

"""
`iso(G,H)` returns a permutation matrix `P` such that `A*P==P*B` where
`A` is the adjacency matrix of `G` and `B` is the adjacency matrix of `H`.
An empty matrix is returned if no such `P` can be found.
"""
function iso(G::SimpleGraph, H::SimpleGraph)
    m = Model()
    n = NV(G)

    P_fail = Array{Int,2}()

    if NV(H) != n || NE(G) != NE(H)
        warn("Two graphs have different number of vertices or edges")
        return P_fail
    end

    @defVar(m,P[1:n,1:n],Bin)

    A = adjacency(G)
    B = adjacency(H)

    for i=1:n
        for k=1:n
            @addConstraint(m,
                           sum{A[i,j]*P[j,k], j=1:n}
                           == sum{P[i,j]*B[j,k], j=1:n}
                           )
        end
        @addConstraint(m, sum{P[i,j], j=1:n}==1)
        @addConstraint(m, sum{P[j,i], j=1:n}==1)
    end

    status = solve(m)
    # println(status)
    if status != :Optimal
        return P_fail
    end
    return round(Int,getValue(P))
end
