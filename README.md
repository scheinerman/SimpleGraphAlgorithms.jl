# SimpleGraphAlgorithms

[![Build Status](https://travis-ci.org/scheinerman/SimpleGraphAlgorithms.jl.svg?branch=master)](https://travis-ci.org/scheinerman/SimpleGraphAlgorithms.jl)


[![codecov.io](http://codecov.io/github/scheinerman/SimpleGraphAlgorithms.jl/coverage.svg?branch=master)](http://codecov.io/github/scheinerman/SimpleGraphAlgorithms.jl?branch=master)


This module provides additional functions for the `SimpleGraphs`
module that rely on integer programming. In addition to requiring the
`SimpleGraphs` module, it also requires `JuMP` and `MathProgBase`
which, in turn, requires that some solvers be loaded. I've used `Cbc`.

**New**: Now requires the `Polynomials` module.

**Note**: Because these functions rely on solving integer linear
  programs, they can be rather slow for large graphs.

## Functions

* `max_indep_set(G)` returns a maximum size independent set of a
  `SimpleGraph`.

* `max_clique(G)` returns a maximum size clique of a `SimpleGraph`.

* `max_matching(G)` returns a maximum size matching of a
  `SimpleGraph`.

* `min_dom_set(G)` returns a smallest dominating set of a
  `SimpleGraph`. That is, a smallest set `S` with the property that
  every vertex of `G` either is in `S` or is adjacent to a vertex of
  `S`.

* `min_vertex_cover(G)` returns a smallest vertex cover of `G`. This
  is a set of vertices `S` such that every edge of `G` has at least
  one end point in `S`.

* `min_edge_cover(G)` returns a smallest edge cover of `G`. This is
  a set of edges `F` such that every vertex of `G` is the end point
  of at least one edge in `F`. **Note**: If `G` has an isolated
  vertex, then no edge cover is possible and error is generated.

* `iso(G,H)` finds an isomorphism between graphs `G` and
  `H`. Specifically, it finds a `Dict` mapping the vertices of `G` to
  the vertices of `H` that gives the isomorphism. If the graphs are
  not isomorphic, an error is raised.

* `iso2(G,H)` has the same functionality as `iso`, but applies various
  preprocessing to speed up the optimization. If the graphs are vertex
  transitive, this probably won't help. But if they have small
  automorphism groups, this will likely speed things up
  considerably. It will also likely detect when the given graphs are
  not isomorphic faster than `iso` will.

* `is_iso(G,H)` checks if the two graphs are isomorphic.

* `is_iso(G,H,d)` checks if the dictionary `d` is an isomorphism
  from `G` to `H`.

* `iso_matrix(G,H)` finds an isomorphism between graphs `G` and
  `H`. Specifically, it finds a permutation matrix `P` such that
  `A*P==P*B` where `A` and `B` are the adjacency matrices of the
  graphs `G` and `H`, respectively. If the graphs are not isomorphic,
  an error is raised.

* `info_map(G)` creates a mapping from the vertices of `G` to 128-bit
  integers. If there is an automorphism between a pair of vertices,
  then they will map to the same value, and the converse is *likely*
  to be true. This is used by `iso2` as part of the preprocessing
  phase.

* `uhash(G)` creates a hash value for the graph `G` with the property
   that isomorphic graphs have the same hash value.

* `color(G,k)` returns a `k`-coloring of `G` (or throws an error if no
  such coloring exists).

* `chromatic_number(G)` returns the least `k` such that `G` is `k`-colorable.

* `chrome_poly(G)` computes the chromatic polynomial of `G`. (See the
  `help` message for more information.)

* `edge_color(G,k)` returns a `k`-edge-coloring of `G`.

* `edge_chromatic_number(G)` returns the edge chromatic number of `G`.


* `ad(G)` returns the average degree of `G`.

* `mad(G)` returns the maximum average degree of `G`.

* `mad_core(G)` returns a subgraph `H` of `G` for which `ad(H)==mad(G)`.

## Examples

```julia
julia> using SimpleGraphs; using SimpleGraphAlgorithms

julia> G = Paley(17)
Paley(17) graph (n=17, m=68)

julia> max_indep_set(G)
Set([7,4,1])

julia> max_clique(G)
Set([3,5,1])

julia> min_dom_set(G)
Set([0,10,3])

julia> max_matching(G)
Set([(2,3),(11,13),(15,16),(0,1),(10,14),(6,7),(4,5),(8,9)])

julia> color(G,6)
Dict{Int64,Int64} with 17 entries:
  2  => 1
  16 => 5
  11 => 6
  0  => 6
  7  => 2
  9  => 1
  10 => 5
  8  => 3
  6  => 6
  4  => 2
  3  => 4
  5  => 3
  13 => 5
  14 => 4
  15 => 3
  12 => 1
  1  => 2
```

Petersen's graph can be described as either the 5,2-Kneser graph or as
the complement of the line graph of K(5).

```julia
julia> G = Kneser(5,2)
Kneser(10,2) graph (n=10, m=15)

julia> H = complement(line_graph(Complete(5)))
SimpleGraph{Tuple{Int64,Int64}} (n=10, m=15)

julia> iso(G,H)
Dict{Set{Int64},Tuple{Int64,Int64}} with 10 entries:
  Set([4,1]) => (4,5)
  Set([4,5]) => (1,5)
  Set([2,1]) => (3,4)
  Set([3,5]) => (1,2)
  Set([2,5]) => (1,3)
  Set([3,1]) => (2,4)
  Set([4,3]) => (2,5)
  Set([4,2]) => (3,5)
  Set([2,3]) => (2,3)
  Set([5,1]) => (1,4)

julia> iso_matrix(G,H)
10x10 Array{Int64,2}:
 0  0  0  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  1
 1  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0
 0  0  0  0  0  1  0  0  0  0
 0  1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0
 0  0  0  0  0  0  1  0  0  0
 0  0  0  0  1  0  0  0  0  0
 0  0  1  0  0  0  0  0  0  0
```

## Setting Solver and its Options

By default, the `Cbc` solver is used for integer programming
and the optimizer does no output.

The function `use_Cbc()` sets the solver to be the `Cbc` solver.
Called as `use_Cbc(true)` causes the solver to be verbose in
it working.

The `Gurobi` solver may used instead. Since this module is not
dependent on `Gurobi`, do this:
```
julia> using Gurobi
julia> use_Gurobi()
```
Alternatively, `use_Gurobi(true)` for extensive output as the
solver does its work.

To switch back to the `Cbc` solver, do `use_Cbc()`.
