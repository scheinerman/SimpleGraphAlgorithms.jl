export kcolor

"""
`kcolor(G,k)`: Return a `k`-coloring of `G` (or error if none exists).
"""
function kcolor(G::SimpleGraph, k::Int)
    VV = vlist(G)
    EE = elist(G)
    n = NV(G)
    m = NE(G)

    if k<1
        error("Number of colors must be positive")
    end

    err_msg = "This graph is not " * string(k) * " colorable"

    # Special cases that don't require integer programming

    result = Dict{vertex_type(G),Int}()

    if k==1
        if NE(G) > 0
            error(err_msg)
        end
        for v in VV
            result[v] = 1
        end
        return result
    end

    if k==2
        return two_color(G)
    end

    MOD = Model()

    @variable(MOD, x[VV,1:k], Bin)

    for v in VV
        @constraint(MOD, sum{x[v,i], i=1:k} == 1)
    end

    for e in EE
        u = e[1]
        v = e[2]
        for i=1:k
            @constraint(MOD, x[u,i] + x[v,i] <= 1)
        end
    end

    status = solve(MOD, suppress_warnings=true)

    if status != :Optimal
        error(err_msg)
    end
    
    
    X = getValue(x)

    for v in VV
        for c = 1:k
            if X[v,c] > 0
                result[v] = c
            end
        end
    end

    return result
   
end
