export color, chromatic_number

"""
`color(G,k)`: Return a `k`-coloring of `G` (or error if none exists).
"""
function color(G::SimpleGraph, k::Int)
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

    MOD = Model(solver=_SOLVER())

    @variable(MOD, x[VV,1:k], Bin)

    for v in VV
        @constraint(MOD, sum(x[v,i] for i=1:k) == 1)
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


    X = getvalue(x)

    for v in VV
        for c = 1:k
            if X[v,c] > 0
                result[v] = c
            end
        end
    end

    return result
end




function chromatic_number(G::SimpleGraph, verb::Bool=false)::Int
    if cache_check(G,:chromatic_number)
        return cache_recall(G,:chromatic_number)
    end

    if NV(G) == 0
        return 0
    end
    if NE(G) == 0
        return 1
    end

    n = NV(G)
    alf = length(max_indep_set(G))
    lb1 = Int(floor(n/alf))
    lb2 = length(max_clique(G))

    lb = max(lb1,lb2)

    k =  chromatic_number_work(G,lb,n+1,verb)
    cache_save(G,:chromatic_number,k)
    return k
end



function chromatic_number_work(G::SimpleGraph, lb::Int, ub::Int, verb::Bool)::Int
    if verb
        print("$lb <= chi(G) <= $ub")
    end
    if lb == ub
        if verb
            println("\tand we're done")
        end
        return lb
    end

    mid = Int(floor((ub+lb)/2))

    # if verb
    #     print("\ttry k = $mid\t")
    # end

    try
        f = color(G,mid)  # success
        if verb
            println("\tgraph is $mid-colorable")
        end
        return chromatic_number_work(G,lb,mid,verb)
    catch
    end
    if verb
        println("\tgraph is not $mid-colorable")
    end
    return chromatic_number_work(G,mid+1,ub,verb)
end
