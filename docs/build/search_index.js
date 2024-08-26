var documenterSearchIndex = {"docs":
[{"location":"#SimpleGraphAlgorithms","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"This module provides additional functions for the SimpleGraphs module that rely on integer programming. In addition to requiring the SimpleGraphs module, it also requires JuMP  and ChooseOptimizer to  select an integer programming solver (defaults to HiGHS).","category":"page"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"New in version 0.6: The functions use_Gurobi and use_Cbc have been removed. The default solver is now HiGHS. To change solvers, use the set_solver function from ChooseOptimizer.","category":"page"},{"location":"#Cliques-and-Independent-Sets","page":"SimpleGraphAlgorithms","title":"Cliques and Independent Sets","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"max_indep_set(G) returns a maximum size independent set of an UndirectedGraph.\nmax_clique(G) returns a maximum size clique of an UndirectedGraph.\nmax_matching(G) returns a maximum size matching of an UndirectedGraph.\nfractional_matching(G) returns a (maximum) fractional matching of the  graph G. This is presented a dictionary mapping edges of G to rational values  in {0, 1/2, 1}. \nkfactor(G,k) returns a k-factor of G. This is a set of edges with the property that every vertex of G is incident with exactly k edges of the set. An error is thrown if no such set exists. (When k==1 this returns a perfect matching.)","category":"page"},{"location":"#Covering-and-Domination","page":"SimpleGraphAlgorithms","title":"Covering and Domination","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"min_dom_set(G) returns a smallest dominating set of an UndirectedGraph. That is, a smallest set S with the property that every vertex of G either is in S or is adjacent to a vertex of S.\nmin_vertex_cover(G) returns a smallest vertex cover of G. This is a set of vertices S such that every edge of G has at least one end point in S.\nmin_edge_cover(G) returns a smallest edge cover of G. This is a set of edges F such that every vertex of G is the end point of at least one edge in F. Note: If G has an isolated vertex, then no edge cover is possible and error is generated.","category":"page"},{"location":"#Isomorphism","page":"SimpleGraphAlgorithms","title":"Isomorphism","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"iso(G,H) finds an isomorphism between graphs G and H. Specifically, it finds a Dict mapping the vertices of G to the vertices of H that gives the isomorphism. If the graphs are not isomorphic, an error is raised. \niso2(G,H) has the same functionality as iso but omits various  preliminary checks. This may be faster for highly symmetric graphs  (e.g., for vertex transitive graphs).\nis_iso(G,H) checks if the two graphs are isomorphic.\nis_iso(G,H,d) checks if the dictionary d is an isomorphism from G to H.\niso_matrix(G,H) finds an isomorphism between graphs G and H. Specifically, it finds a permutation matrix P such that A*P==P*B where A and B are the adjacency matrices of the graphs G and H, respectively. If the graphs are not isomorphic, an error is raised.\nfrac_iso(G,H) finds a fractional isomorphism between the graphs. Specifically,  if A and B are the adjacency matrices of the two graphs, then produce a doubly stochastic matrix S such that A*S == S*B, or throw an error if no such matrix exists.\nis_frac_iso(G,H) returns true if the graphs are fractionally isomorphic and false if not. \nhom(G,H) finds a graph homomorphism from G to H. This is a mapping f (dictionary) with the property that if {u,v} is an edge of G then {f[u],f[v]} is an edge of H. If no homomorphism exists an error is raised.\nhom_check(G,H,d) determines if d is a homomorphism from G to H.\ninfo_map(G) creates a mapping from the vertices of G to 128-bit integers. If there is an automorphism between a pair of vertices, then they will map to the same value, and the converse is likely to be true. \nuhash(G) creates a hash value for the graph G with the property  that isomorphic graphs have the same hash value.","category":"page"},{"location":"#Coloring","page":"SimpleGraphAlgorithms","title":"Coloring","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"vertex_color(G,k) returns a k-coloring of G (or throws an error if no such coloring exists). If k is omitted, the number of colors is χ(G)  (chromatic number).\nvertex_color(G,a,b) returns an a:b-coloring of G (or throws an error if no such coloring exists). An a:b-coloring is a mapping from the vertices of G to b-element subsets of {1,2,...,a} such that adjacent vertices are  assigned disjoint sets. \nchromatic_number(G) returns the least k such that G is k-colorable.\nchromatic_poly(G) computes the chromatic polynomial of G. (See the help message for more information.)\nedge_color(G,k) returns a k-edge-coloring of G.\nedge_chromatic_number(G) returns the edge chromatic number of G.","category":"page"},{"location":"#Connectivity","page":"SimpleGraphAlgorithms","title":"Connectivity","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"min_cut(G) returns a minimum size (vertex) cut set. min_cut(G,s,t) return a smallest set of vertices that separate s and t.\nconnectivity(G) or connectivity(G,s,t) returns the size of such a cut set.\nmin_edge_cut(G) returns a minimum size edge cut set. min_edge_cut(G,s,t) returns a minimum set of edges that separate vertices s and t.\nedge_connectivity(G) or edge_connectivity(G,s,t) returns the size of such an edge cut set.","category":"page"},{"location":"#Maximum-Average-Degree","page":"SimpleGraphAlgorithms","title":"Maximum Average Degree","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"ad(G) returns the average degree of G.\nmad(G) returns the maximum average degree of G.\nmad_core(G) returns a subgraph H of G for which ad(H)==mad(G).","category":"page"},{"location":"#Examples","page":"SimpleGraphAlgorithms","title":"Examples","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"julia> using SimpleGraphs; using SimpleGraphAlgorithms; using ChooseOptimizer; using ShowSet\n\njulia> set_solver_verbose(false)\n[ Info: Setting verbose option for Cbc to false\n\njulia> G = Paley(17)\nPaley (n=17, m=68)\n\njulia> max_indep_set(G)\n{1,4,7}\n\njulia> max_clique(G)\n{3,4,5}\n\njulia> min_dom_set(G)\n{3,6,9}\n\njulia> max_matching(G)\n{(1, 16),(2, 4),(3, 12),(5, 9),(6, 15),(7, 8),(10, 11),(13, 14)}\n\njulia> vertex_color(G,6)\nDict{Int64,Int64} with 17 entries:\n  2  => 3\n  16 => 1\n  11 => 4\n  0  => 4\n  7  => 6\n  9  => 2\n  10 => 1\n  8  => 3\n  6  => 4\n  4  => 6\n  3  => 5\n  5  => 3\n  13 => 1\n  14 => 5\n  15 => 2\n  12 => 2\n  1  => 6","category":"page"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"Petersen's graph can be described as either the 5,2-Kneser graph or as the complement of the line graph of K(5).","category":"page"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"julia> G = Kneser(5,2);\n\njulia> H = complement(line_graph(Complete(5)));\n\njulia> iso(G,H)\nDict{Set{Int64},Tuple{Int64,Int64}} with 10 entries:\n  {1,4} => (1, 5)\n  {2,4} => (1, 4)\n  {2,5} => (3, 4)\n  {1,3} => (2, 5)\n  {3,4} => (1, 2)\n  {1,2} => (4, 5)\n  {3,5} => (2, 3)\n  {4,5} => (1, 3)\n  {2,3} => (2, 4)\n  {1,5} => (3, 5)\n\njulia> iso_matrix(G,H)\n10×10 Array{Int64,2}:\n 0  0  0  0  0  0  0  1  0  0\n 0  0  0  0  0  0  0  0  1  0\n 0  0  0  1  0  0  0  0  0  0\n 0  0  0  0  0  1  0  0  0  0\n 0  0  0  0  1  0  0  0  0  0\n 0  0  0  0  0  0  0  0  0  1\n 1  0  0  0  0  0  0  0  0  0\n 0  1  0  0  0  0  0  0  0  0\n 0  0  0  0  0  0  1  0  0  0\n 0  0  1  0  0  0  0  0  0  0","category":"page"},{"location":"#Using-SimpleGraphAlgorithms-with-Graphs","page":"SimpleGraphAlgorithms","title":"Using SimpleGraphAlgorithms with Graphs","text":"","category":"section"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"SimpleGraphAlgorithms is built to work with UndirectedGraph objects as defined in SimpleGraphs.  To apply these functions to graphs from Julia's Graphs module, you can use SimpleGraphConverter like this:","category":"page"},{"location":"","page":"SimpleGraphAlgorithms","title":"SimpleGraphAlgorithms","text":"julia> using Graphs, SimpleGraphAlgorithms, SimpleGraphs, SimpleGraphConverter\n\njulia> g = circular_ladder_graph(9)\n{18, 27} undirected simple Int64 graph\n\njulia> chromatic_number(UG(g))\n3","category":"page"}]
}