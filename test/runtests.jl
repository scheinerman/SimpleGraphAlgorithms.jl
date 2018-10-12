using Test
using SimpleGraphs, SimpleGraphAlgorithms

G = Paley(17)
A = max_clique(G)
B = max_indep_set(G)
@test length(A) == length(B)
