using SimplePolynomials, SimplePartitions
import Base.length, Base.setindex!, Base.getindex
export chromatic_poly, reset_cpm, size_cpm

const CPM_pair = Tuple{UG,SimplePolynomial}
const CPM_dict = Dict{UInt64,Vector{CPM_pair}}

# This is a device to record previously computed chromatic polynomials
mutable struct ChromePolyMemo
    D::CPM_dict
    function ChromePolyMemo()
        new(CPM_dict())
    end
end

_CPM = ChromePolyMemo()

"""
`reset_cpm()` clears the datastructure of previously computed
chromatic polynomials.
"""
function reset_cpm()
    global _CPM = ChromePolyMemo()
    return true
end

"""
`size_cpm()` reports the number of graphs held in the datastructure of
previously computed chromatic polynomials.
"""
function size_cpm()
    return length(_CPM)
end

function show(io::IO, CPM::ChromePolyMemo)
    print(io, "ChromePolyMemo with $(length(CPM)) graphs")
end

function length(CPM::ChromePolyMemo)
    result = 0
    for ds in keys(CPM.D)
        result += length(CPM.D[ds])
    end
    return result
end

"""
`remember!(CPM,G,P)` save the polynomial `P` associated with `G` in
this data structure (if it wasn't in there already).

Short form: `CPM[G] = P`.
"""
function remember!(CPM::ChromePolyMemo, G::UG, P::SimplePolynomial)
    index::Int128 = uhash(G)  # get the signature for this graph

    # have we seen this graph's uhash before?
    # if no, create a new entry in the Memo
    if !haskey(CPM.D, index)
        CPM.D[index] = [(G, P)]
        return
    end

    # Otherwise, look through the values associated with this uhash to
    # see if we've seen this graph (or one isomorphic to it) before.

    # Recall those graphs with matching uhash
    for (H, Q) in CPM.D[index]
        if is_iso(G, H) # if isomorphic, nothing to add
            return
        end
    end

    # otherwise, add this polynomial to the end of the list
    push!(CPM.D[index], (G, P))
end

setindex!(CPM::ChromePolyMemo, P::SimplePolynomial, G::UG) = remember!(CPM, G, P)


"""
`recall(CPM,G)` look up in CPM to see if we have already computed the
chromatic polynomial of this graph (or something isomorphic to it).

Short form: `P = CPM[G]`.
"""
function recall(CPM::ChromePolyMemo, G::UG)
    ds = uhash(G)
    SG = CPM.D[ds]  # This may throw an error if not found. That's good.

    for (H, Q) in SG
        if is_iso(G, H)
            return Q
        end
    end

    # Not found, so it's an error
    error("Graph not found")
end

getindex(CPM::ChromePolyMemo, G::UG) = recall(CPM, G)

"""
`chromatic_poly(G)` computes the chromatic polynomial of the graph `G`.
If `G` is a directe graph, this returns the chromatic polynomial of
`G`'s underlying simple graph.

This function builds a datastructure to prevent computing the
chromatic polynomial of the same graph twice. To do this, it uses
frequent isomorphism checks.

See `size_cpm` and `reset_cpm`.
"""
function chromatic_poly(G::UG) #, CPM::ChromePolyMemo = _CPM)
    if cache_check(G, :chrome_poly)
        return G.cache[:chrome_poly]
    end

    n::Int = NV(G)
    m::Int = NE(G)

    # Special case: no vertices
    if n == 0
        return SimplePolynomial([1])
    end

    # Special case: no edges
    if m == 0
        x = getx()
        return x^n
    end

    # Special case: complete graph
    if 2m == n * (n - 1)
        x = getx()
        return prod(x - k for k = 0:n-1)
        # return fromroots(0:n-1)
    end

    # Special case: disconnected graph
    comps = parts(components(G))
    if length(comps) > 1
        result = SimplePolynomial([1])
        for A in comps
            H = induce(G, A)
            result *= chromatic_poly(H)
        end
        return result
    end

    # Special case: Tree
    if m == n - 1
        return SimplePolynomial([0, 1]) * SimplePolynomial([-1, 1])^(n - 1)
    end

    # See if we've computed this chromatic polynomial of this graph before
    try
        P = _CPM[G]
        return P
    catch
    end

    # Resort to inclusion/exclusion

    # Find a vertex u of minimum degree
    min_d = minimum(deg(G))
    smalls = filter(v -> deg(G, v) == min_d, G.V)
    u = first(smalls)

    # And choose any neighbor
    v = first(G[u])


    # Special case to speed up handling leaves
    if min_d == 1
        GG = deepcopy(G)
        delete!(GG, u)
        p1 = chromatic_poly(GG)
        P = p1 * SimplePolynomial([-1, 1])
        _CPM[G] = P
        return P
    end

    # p1 is chrome_poly of G-e
    GG = deepcopy(G)
    delete!(GG, u, v)
    p1 = chromatic_poly(GG)

    # p2 is chrome_poly of G/e
    GG = deepcopy(G)
    contract!(GG, u, v)
    p2 = chromatic_poly(GG)

    P = p1 - p2

    _CPM[G] = P # save in case we see this graph again
    G.cache[:chrome_poly] = P  # and save in graph's cache too
    return P
end

chromatic_poly(G::DG) = chromatic_poly(simplify(G))
