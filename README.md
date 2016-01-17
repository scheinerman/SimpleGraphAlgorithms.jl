# SimpleGraphAlgorithms

This module provides additional functions for the `SimpleGraphs`
module that rely on integer programming. In addition to requiring the
`SimpleGraphs` module, it also requires `JuMP` and `MathProgBase`
which, in turn, requires that some solvers be loaded. I've used `GLPK`
and `GLPKMathProgInterface`.

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

* `iso_check(G,H,d)` checks if the dictionary `d` is an isomorphism
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

* `kcolor(G,k)` returns a `k`-coloring of `G` (or throws an error if no
  such coloring exists).
  
* `chrome_poly(G)` computes the chromatic polynomial
of `G`. (See the `help` message for more information.)

## Examples

```julia
julia> using SimpleGraphs; using SimpleGraphAlgorithms

julia> G = Paley(17)
SimpleGraphs.SimpleGraph{Int64} (17 vertices)

julia> max_indep_set(G)
Set([7,4,1])

julia> max_clique(G)
Set([3,5,1])

julia> min_dom_set(G)
Set([0,10,3])

julia> max_matching(G)
Set([(2,3),(11,13),(15,16),(0,1),(10,14),(6,7),(4,5),(8,9)])

julia> kcolor(G,6)
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
SimpleGraphs.SimpleGraph{Set{Int64}} (10 vertices)

julia> H = complement(line_graph(Complete(5)))
SimpleGraphs.SimpleGraph{Tuple{Int64,Int64}} (10 vertices)

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
