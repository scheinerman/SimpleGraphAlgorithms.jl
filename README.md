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


## To do list

* Maximum matching
* Minimum dominating set
* Chromatic number
