# SimpleGraphAlgorithms

This module provides additional functions for the `SimpleGraphs`
module that rely on integer programming. In addition to requiring the
`SimpleGraphs` module, it also requires `JuMP` and `MathProgBase`
which, in turn, requires that some solvers be loaded. I've used `GLPK`
and `GLPKMathProgInterface`.

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
  `H`. Specifically, it finds a permutation matrix `P` such that
  `A*P==P*B` where `A` and `B` are the adjacency matrices of the
  graphs `G` and `H`, respectively. If the graphs are not isomorphic,
  an empty matrix is returned.


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

julia> iso(G,G')
17x17 Array{Int64,2}:
 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0
 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0
 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1
 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0
 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0
 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> max_matching(G)
Set([(2,3),(11,13),(15,16),(0,1),(10,14),(6,7),(4,5),(8,9)])
```



## To do list

* Chromatic number
