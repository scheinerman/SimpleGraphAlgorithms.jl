# SimpleGraphAlgorithms

This module provides additional functions for the `SimpleGraphs`
module that rely on integer programming. In addition to requiring the
`SimpleGraphs` module, it also requires `MathProgBase` which, in turn,
requires that some solvers be loaded. I've used `GLPK` and
`GLPKMathProgInterface`.

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


## To do list

* Chromatic number
