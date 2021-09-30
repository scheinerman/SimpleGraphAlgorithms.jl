module SimpleGraphAlgorithms
using SimpleGraphs
using JuMP
using Cbc, ChooseOptimizer

_SOLVER = Cbc       # this is used by JuMP
_OPTS = Dict()        # to pass options to JuMP optimizers

export use_Cbc, use_Gurobi, use_optimizer


"""
`use_optimizer(choice::Module, options::Dict)` chooses the
optimization package and package options to be used by functions
in `SimpleGraphAlgorithms` (and `SimplePosetAlgorithms`).

See `use_Cbc` and `use_Gurobi`.
"""
function use_optimizer(choice::Module, options::Dict = Dict())
    global _SOLVER = choice
    global _OPTS = deepcopy(options)

    nothing
end

"""
`use_Gurobi(verbose=false)` sets the optimization software
to be the `Gurobi` solver. Requires the user to invoke
`using Gurobi` first.

See also `use_Cbc`.
"""
function use_Gurobi(verbose::Bool = false)
    set_solver(Main.Gurobi)
    set_solver_verbose(verbose)
    nothing
end

"""
`use_Cbc(verbose=false)` sets the optimization software
to be the `Cbc` solver. This is set by default when `SimpleGraphAlgorithms`
is loaded.

See also `use_Gurobi`.
"""
function use_Cbc(verbose::Bool = false)
    set_solver(Cbc)
    set_solver_verbose(verbose)
    nothing
end

# Use Cbc in nonverbose mode until told otherwise
use_Cbc()


export max_indep_set, max_clique, max_matching, min_dom_set
export min_vertex_cover, min_edge_cover

import SimpleGraphs.cache_save

"""
`max_indep_set(G)` returns a maximum size independent set of a
`SimpleGraph`.
"""
function max_indep_set(G::SimpleGraph)
    if cache_check(G, :max_indep_set)
        return cache_recall(G, :max_indep_set)
    end
    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    MOD = Model(get_solver())
    @variable(MOD, x[VV], Bin)
    for e in EE
        u, v = e
        @constraint(MOD, x[u] + x[v] <= 1)
    end
    @objective(MOD, Max, sum(x[v] for v in VV))
    optimize!(MOD)

    X = value.(x)
    A = Set([v for v in VV if X[v] > 0.1])

    cache_save(G, :max_indep_set, A)
    return A
end

"""
`min_vertex_cover(G)` returns a smallest vertex cover `S` of `G`.
This is a smallest possible set of vertices such that every edge
of `G` has one or both end points in `S`.

`min_vertex_cover(G,k)` returns the smallest set of vertices `S`
such that at least `k` edges are indicent with a vertex in `S`.
"""
min_vertex_cover(G) = setdiff(G.V, max_indep_set(G))

function min_vertex_cover(G::SimpleGraph, k::Int)
    MOD = Model(get_solver())
    Vs = vlist(G)
    Es = elist(G)

    @variable(MOD, x[v in Vs], Bin)
    @variable(MOD, y[f in Es], Bin)

    # edge constraints
    for f in Es
        u, v = f
        @constraint(MOD, x[u] + x[v] >= y[f])
    end

    # assure proper size
    @constraint(MOD, sum(y[f] for f in Es) >= k)

    # we want to minimize the sum of the x[v]'s
    @objective(MOD, Min, sum(x[v] for v in Vs))

    optimize!(MOD)
    status = Int(termination_status(MOD))
    if status != 1
        error("Cannot find such a minimum vertex cover")
    end

    XX = value.(x)
    A = Set([v for v in Vs if XX[v] > 0.1])
    return A
end

"""
`max_clique(G)` returns a maximum size clique of a `SimpleGraph`.
"""
function max_clique(G::SimpleGraph)
    if cache_check(G, :max_clique)
        return cache_recall(G, :max_clique)
    end
    result = max_indep_set(G')
    cache_save(G, :max_clique, result)
    return result
end

"""
`max_matching(G)` returns a maximum matching of a `SimpleGraph`.
"""
function max_matching(G::SimpleGraph)
    if cache_check(G, :max_matching)
        return cache_recall(G, :max_matching)
    end
    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    MOD = Model(get_solver())
    @variable(MOD, x[EE], Bin)

    for v in VV
        star = [e for e in EE if e[1] == v || e[2] == v]
        @constraint(MOD, sum(x[e] for e in star) <= 1)
    end

    @objective(MOD, Max, sum(x[e] for e in EE))

    optimize!(MOD)
    X = value.(x)
    M = Set([e for e in EE if X[e] > 0.1])

    cache_save(G, :max_matching, M)
    return M
end


"""
`min_edge_cover(G)` finds a smallest subset `F` of the edges such that
every vertex of `G` is the end point of at least one edge in `F`. An
error is raised if `G` has a vertex of degree 0.
"""
function min_edge_cover(G::SimpleGraph)
    if cache_check(G, :min_edge_cover)
        return cache_recall(G, :min_edge_cover)
    end

    if NV(G) == 0
        T = eltype(G)
        S = Tuple{T,T}
        return Set{S}()
    end

    if minimum(deg(G)) == 0
        error("Graph has an isolated vertex; no minimum edge cover exists.")
    end

    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    MOD = Model(get_solver())
    @variable(MOD, x[EE], Bin)

    for v in VV
        star = [e for e in EE if e[1] == v || e[2] == v]
        @constraint(MOD, sum(x[e] for e in star) >= 1)
    end

    @objective(MOD, Min, sum(x[e] for e in EE))

    optimize!(MOD)
    X = value.(x)
    A = Set([e for e in EE if X[e] > 0.1])

    cache_save(G, :min_edge_cover, A)
    return A
end

"""
`min_dom_set(G)` returns a smallest dominating set of a
`SimpleGraph`. That is, a smallest set `S` with the property that
every vertex of `G` either is in `S` or is adjacent to a vertex of
`S`.
"""
function min_dom_set(G::SimpleGraph)
    if cache_check(G, :min_dom_set)
        return cache_recall(G, :min_dom_set)
    end
    n = NV(G)
    if n == 0
        T = eltype(G)
        return Set{T}()
    end

    VV = vlist(G)

    MOD = Model(get_solver())
    @variable(MOD, x[VV], Bin)

    for v in VV
        Nv = G[v]
        push!(Nv, v)
        @constraint(MOD, sum(x[w] for w in Nv) >= 1)
    end

    @objective(MOD, Min, sum(x[v] for v in VV))
    optimize!(MOD)

    X = value.(x)

    S = Set(v for v in VV if X[v] > 0.1)

    cache_save(G, :min_dom_set, S)
    return S
end

include("iso.jl")
include("hom.jl")
include("vertex-coloring.jl")
include("chrome_poly.jl")
include("edge-coloring.jl")
include("mad.jl")
include("kfactor.jl")
include("connectivity.jl")
include("fmatch.jl")
include("frac_iso.jl")

end # of module SimpleGraphAlgorithms
