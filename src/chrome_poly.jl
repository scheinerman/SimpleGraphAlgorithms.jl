using Polynomials
import Base.length, Base.setindex!, Base.getindex, Base.show
export ChromePolyMemo, chrome_poly

typealias CPM_data Dict{Vector{Int}, Dict{SimpleGraph,Poly{Int}}}

# This is a device to record previously computed chromatic polynomials
type ChromePolyMemo
    D::CPM_data
    function ChromePolyMemo()
        new( CPM_data() )
    end
end

function show(io::IO, CPM::ChromePolyMemo)
    print(io,"ChromePolyMemo with $(length(CPM)) graphs")
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
function remember!(CPM::ChromePolyMemo, G::SimpleGraph, P::Poly{Int})
    ds = deg(G)

    # have we seen this graph's degree sequence before?
    # if no, create a new entry in the Memo
    if !haskey(CPM.D,ds)
        dG = Dict{SimpleGraph, Poly{Int}}()
        dG[G] = P
        CPM.D[ds] = dG
        return
    end

    # Otherwise, look through the values associated with this degree
    # sequence to see if we've seen this graph (or one isomorphic to
    # it) before.

    dG = CPM.D[ds]  # Recall those graphs with matching deg seq
    for H in keys(dG)
        try iso2(G,H)  # if isomorphic, we're done!
            return
        end
    end
    # otherwise, add this polynomial to the dG table
    dG[G] = P
end

setindex!(CPM,P,G) = remember!(CPM,G,P)


"""
`recall(CPM,G)` look up in CPM to see if we have already computed the
chromatic polynomial of this graph (or something isomorphic to it.

Short form: `P = CPM[G]`.
"""
function recall(CPM::ChromePolyMemo, G::SimpleGraph)
    ds = deg(G)
    dG = CPM.D[ds]  # This may throw an error if not found. That's good.

    for H in keys(dG)
        try iso(G,H) # we found a copy of this graph!
            return dG[H]
        end
    end

    # Not found, so it's an error
    error("Graph not found")
end

getindex(CPM,G) = recall(CPM,G)

"""
`chrome_poly(G)` computes the chromatic polynomial of the graph `G`.

This  function  builds  a   datastructure  to  prevent  computing  the
chromatic polynomial  of the  same graph  twice. To  do this,  it uses
frequent isomorphism checks. It is possible to save the datastructure 
for future invocations of `chrome_poly` as follows.

+ First, create a `ChromePolyMemo` object like this: `CPM=ChromePolyMemo()`.
+ Then, call `chrome_poly` with `CPM` as a second argument:
  `P = chrome_poly(G,CPM)`. Note that `CPM` will be changed by this; 
  it will be a database of the chromatic polynomials it has already computed.
+ Later, to compute the chromatic polynomial of another graph, just use
  `chrome_poly(H,CPM)`.
"""

function chrome_poly(G::SimpleGraph, CPM::ChromePolyMemo = ChromePolyMemo())
    n::Int = NV(G)    
    m::Int = NE(G)

    # Special case: no vertices
    if n==0
        return Poly([1])
    end

    # Special case: no edges
    if m==0
        return Poly([zeros(Int,n);1])
    end

    # Special case: complete graph
    if 2m == n*(n-1)
        return poly(0:n-1)
    end

    # Special case: disconnected graph
    comps = components(G)
    if length(comps) > 1
        result = Poly([1])
        for A in comps
            H = induce(G,A)
            result *= chrome_poly(H)
        end
        return result
    end 

    # Special case: Tree
    if m == n-1
        return Poly([0,1]) * Poly([-1,1])^(n-1)
    end

    # See if we've computed this chromatic polynomial of this graph before
    try P = CPM[G]
        return P
    end

    # Resort to inclusion/exclusion

    # Find a vertex u of minimum degree
    min_d = minimum(deg(G))
    smalls = filter(v -> deg(G,v) == min_d, G.V)
    u = first(smalls)

    # And choose any neighbor
    v = first(G[u])

    
    # Special case to speed up handling leaves
    if min_d==1
        delete!(G,u)
        p1 = chrome_poly(G,CPM)
        add!(G,u,v)
        P = p1 * Poly([-1,1])
        CPM[G] = P
        return P
    end

    # p1 is chrome_poly of G-e
    delete!(G,u,v)
    p1 = chrome_poly(G,CPM)
    add!(G,u,v)

    # p2 is chrome_poly of G/e
    GG = deepcopy(G)
    contract!(GG,u,v)
    p2 = chrome_poly(GG,CPM)

    P = p1 - p2

    CPM[G] = P # save in case we see this graph again

    return P
end
    
