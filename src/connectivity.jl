export edge_connectivity, min_edge_cut
export connectivity, min_cut


# min_cut(G,s,t) -- min cut separating s and t
# min_cut(G) -- min cut separating the graph
# flag == true means we want s/t cut

"""
`min_cut(G)` returns a minimum size set of vertices that disconnects `G`.

`min_cut(G,s,t)` returns a minimum size cut that separates `s` and `t`.
"""
function min_cut(G::UG{T}, s::T, t::T, flag::Bool = true)::Set{T} where {T}
    n = NV(G)
    m = NE(G)
    n * (n - 1) != 2m || error("Graph must not be complete")

    if flag
        has(G, s) || error("$s is not a vertex of this graph")
        has(G, t) || error("$t is not a vertex of this graph")
        s != t || error("two vertices must be different")
        !has(G, s, t) || error("vertices $s and $t cannot be adjacent")
    else
        if cache_check(G, :min_cut)
            return cache_recall(G, :min_cut)
        end

        if !is_connected(G)
            X = Set{T}()
            cache_save(G, :min_cut, X)
            return X
        end
    end

    VV = vlist(G)
    EE = elist(G)

    MOD = Model(get_solver())

    @variable(MOD, a[VV], Bin)  # in part 1
    @variable(MOD, b[VV], Bin)  # in part 2
    @variable(MOD, c[VV], Bin)  # in cut set

    for v in VV
        @constraint(MOD, a[v] + b[v] + c[v] == 1)
    end

    if flag
        @constraint(MOD, a[s] == 1)
        @constraint(MOD, b[t] == 1)
    end

    @constraint(MOD, sum(a[v] for v in VV) >= 1)
    @constraint(MOD, sum(b[v] for v in VV) >= 1)

    for e in EE
        v, w = e
        @constraint(MOD, a[v] + b[w] <= 1)
        @constraint(MOD, a[w] + b[v] <= 1)
    end

    @objective(MOD, Min, sum(c[v] for v in VV))

    optimize!(MOD)
    status = Int(termination_status(MOD))

    C = value.(c)
    X = Set(v for v in VV if C[v] > 0.5)

    if !flag
        cache_save(G, :min_cut, X)
    end
    return X
end

function min_cut(G::UG)::Set
    n = NV(G)
    n > 1 || error("Graph must have at least two vertices")
    s = first(G.V)
    return min_cut(G, s, s, false)
end

"""
`connectivity(G)` returns the (vertex) connectivity of `G`, i.e.,
the minimum size of a vertex cut set. If `G` is a complete graph with
`n` vertices, the connectivity is `n-1` (or `0` for the empty graph).
"""
function connectivity(G::UG)::Int
    n = NV(G)
    m = NE(G)
    if n * (n - 1) == 2m   # graph is complete
        return max(n - 1, 0)  # what is the connectivity of K_0?
    end

    return length(min_cut(G))
end


# min_edge_cut(G,s,t) -- min edge cut separating s and t
# min_edge_cut(G) -- min edge cut separating the graph
# flag == true means we want s/t cut

"""
`min_edge_cut(G)` returns a minimum size set of edges whose removal
disconnects `G`. The graph must have at least two vertices.

`min_edge_cut(G,s,t)` is a minimum size set of edges whose removal
separates vertices `s` and `t`.
"""
function min_edge_cut(G::UG{T}, s::T, t::T, flag::Bool = true)::Set{Tuple{T,T}} where {T}
    n = NV(G)
    n > 1 || error("Graph must have at least two vertices")

    if flag
        has(G, s) || error("$s is not a vertex of this graph")
        has(G, t) || error("$t is not a vertex of this graph")
        s != t || error("vertices must not be the same")
    else
        if cache_check(G, :min_edge_cut)
            return cache_recall(G, :min_edge_cut)
        end

        if !is_connected(G)
            X = Set{Tuple{T,T}}()
            cache_save(G, :min_edge_cut, X)
            return X
        end
    end

    VV = vlist(G)
    EE = elist(G)

    MOD = Model(get_solver())

    @variable(MOD, a[VV], Bin)  # in part 1
    @variable(MOD, b[VV], Bin)  # in part 2
    @variable(MOD, c[EE], Bin)  # span between the parts

    for v in VV
        @constraint(MOD, a[v] + b[v] == 1)
    end

    @constraint(MOD, sum(a[v] for v in VV) >= 1)
    @constraint(MOD, sum(b[v] for v in VV) >= 1)

    if flag
        @constraint(MOD, a[s] == 1)
        @constraint(MOD, b[t] == 1)
    end


    for e in EE
        u, v = e
        @constraint(MOD, a[u] + b[v] - 1 <= c[e])
        @constraint(MOD, a[v] + b[u] - 1 <= c[e])
    end

    @objective(MOD, Min, sum(c[e] for e in EE))

    optimize!(MOD)
    status = Int(termination_status(MOD))

    C = value.(c)

    X = Set(e for e in EE if C[e] > 0.5)

    if !flag
        cache_save(G, :min_edge_cut, X)
    end
    return X
end

function min_edge_cut(G::UG)::Set
    n = NV(G)
    n > 1 || error("Graph must have at least two vertices")
    s = first(G.V)
    return min_edge_cut(G, s, s, false)
end

"""
`edge_connectivity(G)` returns the size of a minimum edge cut of `G`.

`edge_connectivity(G,s,t)` determines the minimum size of an edge cut
separating `s` and `t`.
"""
edge_connectivity(G::UG)::Int = length(min_edge_cut(G))
edge_connectivity(G::UG, s, v) = length(min_edge_cut(G, s, t))
