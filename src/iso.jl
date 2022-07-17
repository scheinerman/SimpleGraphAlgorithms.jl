export iso, iso2, iso_matrix, is_iso, info_map, uhash, degdeg, fast_iso_test


const iso_err_msg = "The graphs are not isomorphic"


"""
`iso_matrix(G,H)` returns a permutation matrix `P` such that
`A*P==P*B` where `A` is the adjacency matrix of `G` and `B` is the
adjacency matrix of `H`.  If the graphs are not isomorphic, an
error is raised.
"""
function iso_matrix(G::UG, H::UG)
    n = NV(G)

    if !fast_iso_test_basic(G, H)
        error(iso_err_msg)
    end

    m = Model(get_solver())

    @variable(m, P[1:n, 1:n], Bin)
    A = adjacency(G)
    B = adjacency(H)

    for i = 1:n
        for k = 1:n
            @constraint(
                m,
                sum(A[i, j] * P[j, k] for j = 1:n) == sum(P[i, j] * B[j, k] for j = 1:n)
            )
        end
        @constraint(m, sum(P[i, j] for j = 1:n) == 1)
        @constraint(m, sum(P[j, i] for j = 1:n) == 1)
    end

    optimize!(m)
    status = Int(termination_status(m))

    if status != 1
        error(iso_err_msg)
    end

    return Int.(value.(P))
end



"""
`is_iso(G,H,d)` checks if `d` is an isomorphism from `G` to `H`.
"""
function is_iso(G::UG, H::UG, d::Dict)::Bool
    n = NV(G)

    # standard quick check
    if n != NV(H) || NE(G) != NE(H) || deg(G) != deg(H) || length(d) != n
        return false
    end

    # check keys and values in d match VG and VH, respectively
    for k in keys(d)
        if !has(G, k)
            return false
        end
    end

    for v in values(d)
        if !has(H, v)
            return false
        end
    end

    # Check edges match. Only need to go one way since we determined
    # that the two graphs have the same number of edges.
    EE = elist(G)
    for e in EE
        v = e[1]
        w = e[2]
        if !has(H, d[v], d[w])
            return false
        end
    end

    return true
end

"""
`is_iso(G,H)` returns `true` if the two graphs are isomorphic and `false`
otherwise.
"""
function is_iso(G::UG, H::UG)::Bool
    try
        f = iso(G, H)
    catch
        return false
    end
    return true
end

"""
`fast_iso_test_basic(G,H)` is a very quick test to rule out graphs being
isomorphic by considering the number of vertices, number of edges, and
degree sequences. Returns `false` if the graphs fail this very basic
test of isomorphism. A `true` result does *not* imply the graphs are
isomorphic.
"""
function fast_iso_test_basic(G::UG, H::UG)::Bool
    if NV(G) != NV(H) || NE(G) != NE(H) || deg(G) != deg(H)
        return false
    end
    return true
end

"""
`fast_iso_test(G,H)` applies various heuristics to see if
the graphs *might* be isomorphic. A `false` return certifies the
graphs are **not** isomorphic; a `true` result indicates the
graphs might be (indeed, likely are) isomorphic.
"""
function fast_iso_test(G::UG, H::UG)::Bool
    return fast_iso_test_basic(G, H) && uhash(G) == uhash(H)
end

"""
`degdeg(G)` creates an n-by-n matrix where the nonzero entries in a
row are the degrees of the neighbors of that vertex. The rows are
lexicographically sorted.
"""
function degdeg(G::UG)
    if cache_check(G, :degdeg)
        return cache_recall(G, :degdeg)
    end

    n = NV(G)
    result = zeros(Int, n, n)
    VV = vlist(G)

    for k = 1:n
        v = VV[k]
        dv = [deg(G, w) for w in G[v]]
        dv = sort(dv)
        dv = [dv; zeros(Int, n - deg(G, v))]
        result[k, :] = dv
    end

    result = sortslices(result, dims = 1)
    cache_save(G, :degdeg, result)
    return result
end


"""
`info_map(G)`:
We create a dictionary mapping the vertices of `G` to 128-bit integer
values in such a way that twin vertices *will* have the same value
but, we hope, nontwin vertices will have different values. (By *twin*
we mean a pair of vertices such that there is an automorphism of the
graph mapping one to the other.)
"""
function info_map(G::UG)
    if cache_check(G, :info_map)
        return cache_recall(G, :info_map)
    end
    n = NV(G)
    VV = vlist(G)
    H = G'

    # Section 1: Neighborhood degrees
    M1 = zeros(Int, n, n - 1)
    for k = 1:n
        v = VV[k]
        rv = [deg(G, w) for w in G[v]]
        rv = [sort(rv); zeros(Int, n - 1 - deg(G, v))]
        M1[k, :] = rv
    end

    # Section 2: Complement degrees
    M2 = zeros(Int, n, n - 1)
    for k = 1:n
        v = VV[k]
        rv = [deg(H, w) for w in H[v]]
        rv = [sort(rv); zeros(Int, n - 1 - deg(H, v))]
        M2[k, :] = rv
    end

    # Section 3: Distances
    M3 = zeros(Int, n, n - 1)
    for k = 1:n
        v = VV[k]
        dv = dist(G, v)
        M3[k, :] = sort(collect(values(dv)))[1:n-1]
    end

    # Section 4: Distances in complement
    M4 = zeros(Int, n, n - 1)
    for k = 1:n
        v = VV[k]
        dv = dist(H, v)
        M4[k, :] = sort(collect(values(dv)))[1:n-1]
    end

    M = [M1 M2 M3 M4]

    # println(sortrows(M))  # debug

    T = eltype(G)
    d = Dict{T,Int128}()
    for k = 1:n
        v = VV[k]
        d[v] = Int128(hash(M[k, :]))
    end
    cache_save(G, :info_map, d)
    return d
end


"""
`iso(G,H)` finds an isomorphism from `G` to `H` (or throws an error if
the graphs are not isomorphic). Returns a dictionary mapping the
vertices of `G` to `H`.

See also: `iso2`.
"""
function iso(G::UG, H::UG)

    # quickly rule out some easily detected nonisomorphic graphs
    if !fast_iso_test(G, H) || uhash(G) != uhash(H)
        error(iso_err_msg)
    end

    VG = vlist(G)
    VH = vlist(H)
    n = NV(G)

    TG = eltype(G)
    TH = eltype(H)

    dG = info_map(G)
    dH = info_map(H)

    # Create mappings from vertex key values back to the vertices themselves

    xG = Dict{Int128,Set{TG}}()
    xH = Dict{Int128,Set{TH}}()

    valsG = sort(collect(values(dG)))
    valsH = sort(collect(values(dH)))


    for x in valsG
        xG[x] = Set{TG}()
        xH[x] = Set{TH}()
    end

    for v in VG
        x = dG[v]
        push!(xG[x], v)
    end

    for v in VH
        x = dH[v]
        push!(xH[x], v)
    end

    MOD = Model(get_solver())

    @variable(MOD, x[VG, VH], Bin)

    for v in VG
        @constraint(MOD, sum(x[v, VH[k]] for k = 1:n) == 1)
    end

    for w in VH
        @constraint(MOD, sum(x[VG[k], w] for k = 1:n) == 1)
    end

    for v in VG
        for w in VH
            @constraint(
                MOD,
                sum(has(G, v, VG[k]) * x[VG[k], w] for k = 1:n) ==
                sum(x[v, VH[k]] * has(H, VH[k], w) for k = 1:n)
            )
        end
    end


    # Add contraints based on vertex numbers

    for val in unique(valsG)
        SG = collect(xG[val])
        SH = collect(xH[val])
        s = length(SG)
        @constraint(MOD, sum(x[u, v] for u in SG for v in SH) == s)
    end

    optimize!(MOD)
    status = Int(termination_status(MOD))

    if status != 1
        error(iso_err_msg)
    end

    X = value.(x)

    result = Dict{eltype(G),eltype(H)}()

    for v in VG
        for w in VH
            if X[v, w] > 0
                result[v] = w
            end
        end
    end
    return result
end


"""
`iso2(G,H)` is another version of `iso(G,H)` that does much less preprocessing.
It may be faster for highly symmetric graphs (e.g., vertex transitive graphs).
"""
function iso2(G::UG{S}, H::UG{T}) where {S,T}
    if NV(G) != NV(H)
        error(iso_err_msg)
    end

    m = Model(get_solver())


    @variable(m, P[G.V, H.V], Bin)

    for v in G.V
        for x in H.V
            @constraint(
                m,
                sum(G[v, w] * P[w, x] for w in G.V) == sum(P[v, y] * H[y, x] for y in H.V)
            )
        end
    end

    for v in G.V
        @constraint(m, sum(P[v, x] for x in H.V) == 1)
    end

    for x in H.V
        @constraint(m, sum(P[v, x] for v in G.V) == 1)
    end

    optimize!(m)
    status = Int(termination_status(m))

    if status != 1
        error(iso_err_msg)
    end

    PP = Int.(value.(P))

    result = Dict{S,T}()

    for v in G.V
        for x in H.V
            if PP[v, x] > 0
                result[v] = x
            end
        end
    end


    return result

end






using LinearAlgebra


function matrix_moments(A, max_exp::Int = 10)
    result = zeros(Int, max_exp)
    B = copy(A)

    for t = 1:max_exp
        B = A * B
        result[t] = tr(B)
    end
    return result
end




"""
`uhash(G)` creates an `UInt64` hash of the graph such that isomorphic
graphs will produce the same value. We hope that nonisomorphic graphs
will create different values, but, alas, that need not be the case.
"""
function uhash(G::UG)
    if cache_check(G, :uhash)
        return cache_recall(G, :uhash)
    end

    v1 = sort(collect(values(info_map(G))))

    A = adjacency(G)
    B = deepcopy(A)
    depth = min(10, NV(G))

    v2 = matrix_moments(adjacency(G))
    v3 = matrix_moments(laplace(G))


    vv = [v1; v2; v3]


    x = hash(vv)
    cache_save(G, :uhash, x)
    return x
end
