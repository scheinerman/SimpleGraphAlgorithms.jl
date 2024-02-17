using Test
using SimpleGraphs, SimpleGraphAlgorithms, Polynomials

# use_GLPK()
# G = Paley(17)

# Cliques and independent sets
@testset "Clique/Indep" begin
    G = Paley(17)
    A = max_clique(G)
    B = max_indep_set(G)
    @test length(A) == length(B)
end

# Isomorphism
@testset "Isomorphism" begin
    G = Paley(17)
    f = iso(G, G')
    @test length(f) == 17
    @test is_iso(G, G')

    # fractional Isomorphism
    G = Dodecahedron()
    H = RandomRegular(20, 3)
    @test is_frac_iso(G, H)
    H = RandomRegular(20, 4)
    @test !is_frac_iso(G, H)

    d = hom(Cube(3), Complete(2))
    @test d["000"] != d["001"]
end

# Edge matching
@testset "Matching" begin
    G = Paley(17)
    M = max_matching(G)
    @test length(M) == 8


    d = fractional_matching(Cycle(5))
    @test sum(values(d)) == 5 // 2

    H = Dodecahedron()
    X = kfactor(H, 3)
    @test length(X) == 30
end

# Vertex coloring
@testset "Coloring" begin
    G = Cycle(7)
    f = vertex_color(G, 6)
    @test length(f) == NV(G)

    # Edge coloring
    @test edge_chromatic_number(G) == 3
    f = edge_color(G, 3)
    @test length(f) == NE(G)
end

# Domination
@testset "Domination/Covering" begin
    G = Paley(17)
    A = min_dom_set(G)
    @test length(A) == 3

    # Covering
    A = min_vertex_cover(G)
    @test length(A) == 14
    A = min_vertex_cover(G, 10)
    @test length(A) == 2
    A = min_edge_cover(G)
    @test length(A) == 9

    # MAD
    @test mad(G) == maximum(deg(G))
end

@testset "Coloring" begin
    # Chromatic Polynomial
    f = chromatic_poly(Cycle(5))
    @test f(2) == 0
    @test f(3) == 30

    # Chromatic number
    @test chromatic_number(Petersen()) == 3

    G = Spindle()
    add!(G, 4, 0)
    H = Spindle()
    add!(H, 7, 0)
    @test fast_iso_test(G, H)

    d = vertex_color(Spindle(), 7, 2)
    @test length(d[1] ∩ d[2]) == 0
end

# edge and vertex connectivity
@testset "Connectivity" begin
    G = Complete(8, 8)'
    add!(G, 1, 9)
    add!(G, 1, 10)
    @test edge_connectivity(G) == 2
    @test connectivity(G) == 1
end

