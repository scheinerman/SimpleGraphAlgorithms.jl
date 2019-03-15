# This will house vertex and edge connectivity functions.
# WARNING: These implementations are TERRIBLE. This is just to get something
# going.

export edge_connectivity

"""
`edge_connectivity(G,s,t)` determines the minimum size of an edge cut
separating `s` and `t`.
"""
function edge_connectivity(G::SimpleGraph, s, t, verbose::Bool=false)::Int
    has(G,s) || error("$s is not a vertex of this graph")
    has(G,t) || error("$t is not a vertex of this graph")
    s!=t || error("source and sink cannot be the same")

    if length(find_path(G,s,t)) == 0
        return 0
    end

    VV = vlist(G)
    n = NV(G)

    MOD = Model(with_optimizer(_SOLVER.Optimizer; _OPTS...))

    @variable(MOD, x[VV,VV], Bin)  # X[i,j]=1 if unit flow from i to j

    # no flow on a self loop
    for v in VV
        @constraint(MOD,x[v,v]==0)
    end

    # no flow on a non edge
    for u in VV
        for v in VV
            if u!=v && !has(G,u,v)
                @constraint(MOD, x[u,v]==0)
            end
        end
    end

    # don't have flows on both antiparallel edges
    for u in VV
        for v in VV
            if u!=v
                @constraint(MOD,x[u,v]+x[v,u] <= 1)
            end
        end
    end

    # flow in = flow out  everywhere except s,t
    for v in VV
        if v==s || v==t
            continue
        end

        Nv = G[v]
        @constraint(MOD, sum(x[v,w] for w in Nv) == sum(x[w,v] for w in Nv))
    end

    Ns = G[s]
    Nt = G[t]

    @constraint(MOD, sum(x[s,w] for w in Ns) == sum(x[w,t] for w in Nt))

    @objective(MOD, Max, sum(x[s,w] for w in Ns))

    optimize!(MOD)
    status = Int(termination_status(MOD))

    X = value.(x)


    if verbose
        for u in VV
            for v in VV
                if X[u,v] > 0.4
                    println("($u,$v)")
                end
            end
        end
    end

    return Int(objective_value(MOD))

end
