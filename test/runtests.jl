using Test
using SimpleGraphs, SimpleGraphAlgorithms

set_solver(:Cbc)
G = Paley(17)
A = max_clique(G)
B = max_indep_set(G)
@test length(A) == length(B)
