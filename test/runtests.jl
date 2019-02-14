using Test
using SimpleGraphs, SimpleGraphAlgorithms, Polynomials

G = Paley(17)

# Cliques and independent sets
A = max_clique(G)
B = max_indep_set(G)
@test length(A) == length(B)

# Isomorphism
f = iso(G,G')
@test length(f) == 17
f = iso2(G,G')
@test length(f) == 17
@test is_iso(G,G')

# Edge matching
M = max_matching(G)
@test length(M) == 8

# Vertex coloring
f = color(G,6)
@test length(f) == 17

# Edge coloring
@test edge_chromatic_number(G) == 9
f = edge_color(G,9)
@test length(f) == NE(G)

# Domination
A = min_dom_set(G)
@test length(A) == 3

# Covering
A = min_vertex_cover(G)
@test length(A) == 14
A = min_vertex_cover(G,10)
@test length(A) == 2
A = min_edge_cover(G)
@test length(A) == 9

# MAD
@test mad(G) == maximum(deg(G))

# Chromatic Polynomial
f = chrome_poly(Cycle(5))
@test f(2) == 0
@test f(3) == 30

# Chromatic number
@test chromatic_number(Petersen()) == 3

G = Spindle()
add!(G,4,0)
H = Spindle()
add!(H,7,0)
@test fast_iso_test(G,H)
f = iso2(G,H)
@test f[0] == 0
g = iso(G,H)
@test g[0] == 0
