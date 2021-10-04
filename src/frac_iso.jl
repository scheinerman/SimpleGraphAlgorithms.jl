export frac_iso, is_frac_iso

"""
    frac_iso(G,H)

Find a doubly stochastic matrix `S` so that `AS=SB` where `A` and `B`
are the adjacency matrices of `G` and `H` respectively.
"""
function frac_iso(G::SimpleGraph, H::SimpleGraph)
    if deg(G) != deg(H) || degdeg(G) != degdeg(H) # quick basic check before LP
        error("The graphs are not fractionally isomorphic")
    end
    A = adjacency(G)
    B = adjacency(H)
    return frac_iso(A, B)
end

function frac_iso(A::Matrix, B::Matrix)
    n, c = size(A)
    nn, cc = size(B)
    @assert n == c && nn == cc "Matrices must be square"
    @assert n == nn "Matrices must have the same size"

    m = Model(get_solver())

    @variable(m, S[1:n, 1:n])

    for i = 1:n
        for k = 1:n
            @constraint(
                m,
                sum(A[i, j] * S[j, k] for j = 1:n) == sum(S[i, j] * B[j, k] for j = 1:n)
            )
        end
        @constraint(m, sum(S[i, j] for j = 1:n) == 1)
        @constraint(m, sum(S[j, i] for j = 1:n) == 1)
    end

    for i = 1:n
        for j = 1:n
            @constraint(m, S[i, j] >= 0)
        end
    end

    optimize!(m)
    status = Int(termination_status(m))

    if status != 1
        error("The graphs are not fractionally isomorphic")
    end

    return value.(S)

end

"""
    is_frac_iso(G,H)
Test if two graphs are fractionally isomorphic.
"""
function is_frac_iso(G::SimpleGraph, H::SimpleGraph)::Bool
    try
        f = frac_iso(G, H)
        return true
    catch
        return false
    end
    false
end
