# SimpleGraphAlgorithms

[![Build Status](https://travis-ci.com/scheinerman/SimpleGraphAlgorithms.jl.svg?branch=master)](https://travis-ci.com/scheinerman/SimpleGraphAlgorithms.jl)


This module provides additional functions for the `SimpleGraphs`
module that rely on integer programming. In addition to requiring the
`SimpleGraphs` module, it also requires `JuMP` and `MathProgBase`
which, in turn, requires that some solvers be loaded. I've used `Cbc`.

**New**: As of version 0.4.2, 
this module uses the `SimplePolynomials` module instead of `Polynomials`.

**Note**: Because these functions rely on solving integer linear
  programs, they can be rather slow for large graphs.

**Note**: Some name changes as of version 0.4.4:
* `color` is now `vertex_color` to avoid conflict with `Plots`.
* `chome_poly` is now `chromatic_poly`.

## Functions

#### Cliques and independent sets

* `max_indep_set(G)` returns a maximum size independent set of a
  `SimpleGraph`.

* `max_clique(G)` returns a maximum size clique of a `SimpleGraph`.

* `max_matching(G)` returns a maximum size matching of a
  `SimpleGraph`.

* `kfactor(G,k)` returns a `k`-factor of `G`. This is a set of edges
  with the property that every vertex of `G` is incident with exactly `k`
  edges of the set. An error is thrown if no such set exists.
  (When `k==1` this returns a perfect matching.)

#### Covering and domination

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


#### Isomorphism

* `iso(G,H)` finds an isomorphism between graphs `G` and
  `H`. Specifically, it finds a `Dict` mapping the vertices of `G` to
  the vertices of `H` that gives the isomorphism. If the graphs are
  not isomorphic, an error is raised.

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
  to be true. 

* `uhash(G)` creates a hash value for the graph `G` with the property
   that isomorphic graphs have the same hash value.

#### Coloring

* `vertex_color(G,k)` returns a `k`-coloring of `G` (or throws an error if no
  such coloring exists).

* `chromatic_number(G)` returns the least `k` such that `G` is `k`-colorable.

* `chromatic_poly(G)` computes the chromatic polynomial of `G`. (See the
  `help` message for more information.)

* `edge_color(G,k)` returns a `k`-edge-coloring of `G`.

* `edge_chromatic_number(G)` returns the edge chromatic number of `G`.


#### Connectivity

* `min_cut(G)` returns a minimum size (vertex) cut set. `min_cut(G,s,t)`
return a smallest set of vertices that separate `s` and `t`.

* `connectivity(G)` or `connectivity(G,s,t)` returns the size of such a cut set.

* `min_edge_cut(G)` returns a minimum size edge cut set.
`min_edge_cut(G,s,t)` returns a minimum set of edges that separate vertices
`s` and `t`.

* `edge_connectivity(G)` or `edge_connectivity(G,s,t)` returns the size of
such an edge cut set.


#### Maximum average degree

* `ad(G)` returns the average degree of `G`.

* `mad(G)` returns the maximum average degree of `G`.

* `mad_core(G)` returns a subgraph `H` of `G` for which `ad(H)==mad(G)`.

## Examples

```julia
julia> using SimpleGraphs; using SimpleGraphAlgorithms; using ChooseOptimizer; using ShowSet

julia> set_solver_verbose(false)
[ Info: Setting verbose option for Cbc to false

julia> G = Paley(17)
Paley (n=17, m=68)

julia> max_indep_set(G)
{1,4,7}

julia> max_clique(G)
{3,4,5}

julia> min_dom_set(G)
{3,6,9}

julia> max_matching(G)
{(1, 16),(2, 4),(3, 12),(5, 9),(6, 15),(7, 8),(10, 11),(13, 14)}

julia> vertex_color(G,6)
Dict{Int64,Int64} with 17 entries:
  2  => 3
  16 => 1
  11 => 4
  0  => 4
  7  => 6
  9  => 2
  10 => 1
  8  => 3
  6  => 4
  4  => 6
  3  => 5
  5  => 3
  13 => 1
  14 => 5
  15 => 2
  12 => 2
  1  => 6
```

Petersen's graph can be described as either the 5,2-Kneser graph or as
the complement of the line graph of K(5).

```julia
julia> G = Kneser(5,2);

julia> H = complement(line_graph(Complete(5)));

julia> iso(G,H)
Dict{Set{Int64},Tuple{Int64,Int64}} with 10 entries:
  {1,4} => (1, 5)
  {2,4} => (1, 4)
  {2,5} => (3, 4)
  {1,3} => (2, 5)
  {3,4} => (1, 2)
  {1,2} => (4, 5)
  {3,5} => (2, 3)
  {4,5} => (1, 3)
  {2,3} => (2, 4)
  {1,5} => (3, 5)

julia> iso_matrix(G,H)
10×10 Array{Int64,2}:
 0  0  0  0  0  0  0  1  0  0
 0  0  0  0  0  0  0  0  1  0
 0  0  0  1  0  0  0  0  0  0
 0  0  0  0  0  1  0  0  0  0
 0  0  0  0  1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  1
 1  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  0
 0  0  1  0  0  0  0  0  0  0
```

## Setting Solver and its Options

By default, the `Cbc` solver is used for integer programming
and the optimizer does no output.

The function `use_Cbc()` sets the solver to be the `Cbc` solver.
Called as `use_Cbc(true)` causes the solver to be verbose in
its working.

The `Gurobi` solver may used instead. Since this module is not
dependent on `Gurobi`, do this:
```
julia> using Gurobi
julia> use_Gurobi()
```
Alternatively, `use_Gurobi(true)` for extensive output as the
solver does its work.

To switch back to the `Cbc` solver, do `use_Cbc()`.

These functions rely on my `ChooseOptimizer` module.
