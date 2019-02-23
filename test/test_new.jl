using Cbc, MathProgBase, SimpleGraphs, SimpleGraphAlgorithms, Test

function old_max_matching(G::SimpleGraph)
    m = NE(G)
    if m == 0
      T = vertex_type(G)
      S = Tuple{T,T}
      return Set{S}()
    end
    A = incidence(G,false)
    c = ones(m)

    X = mixintprog(-c,A,'<',1,:Int,0,1,CbcSolver())

    indices = findall(x->x!=0,round.(Int,X.sol))
    EE = elist(G)
    result = Set(EE[indices])
    # cache_save(G,:max_matching,result)
    return result
end

function old_max_indep_set(G::SimpleGraph)
    if NV(G)==0
      T = vertex_type(G)
      return Set{T}()
    end
    # if cache_check(G,:max_indep_set)
    #   return cache_recall(G,:max_indep_set)
    # end
    n = NV(G)
    A = incidence(G,false)'
    c = ones(n)

    X = mixintprog(-c,A,'<',1,:Int,0,1,CbcSolver())

    indices = findall(x->x!=0,round.(Int,X.sol))
    VV = vlist(G)
    result = Set(VV[indices])
    # cache_save(G,:max_indep_set,result)
    return result
end


function old_min_edge_cover(G::SimpleGraph)
    # if cache_check(G,:min_edge_cover)
    #   return cache_recall(G,:min_edge_cover)
    # end

    if NV(G) == 0
      T = vertex_type(G)
      S = Tuple{T,T}
      return Set{S}()
    end

    if minimum(deg(G))==0
        error("Graph has an isolated vertex; no minimum edge cover exists.")
    end
    m = NE(G)
    M = incidence(G,false)
    c = ones(m)

    X = mixintprog(c,M,'>',1,:Int,0,1,CbcSolver())

    indices = findall(x->x!=0,round.(Int,X.sol))
    EE = elist(G)
    result = Set(EE[indices])
    # cache_save(G,:min_edge_cover,result)
    return result
end


function old_min_dom_set(G::SimpleGraph)
    # if cache_check(G,:min_dom_set)
    #   return cache_recall(G,:min_dom_set)
    # end
    n = NV(G)
    if n==0
      T = vertex_type(G)
      return Set{T}()
    end

    # A = adjacency(G)+eye(Int,n)

    A = adjacency(G)
    for j=1:n
        A[j,j] = 1
    end
    c = ones(n)

    X = mixintprog(c,A,'>',1,:Int,0,1,CbcSolver())

    indices = findall(x->x!=0,round.(Int,X.sol))
    VV = vlist(G)
    result = Set(VV[indices])
    # cache_save(G,:min_dom_set,result)
    return result
end

function func_test(G::SimpleGraph, old_f::Function, new_f::Function)::Bool
    A = old_f(G)
    B = new_f(G)
    return length(A)==length(B)
end

function quad_test(G::SimpleGraph)
    func_test(G,max_matching,old_max_matching) &&
    func_test(G,max_indep_set,old_max_indep_set) &&
    func_test(G,old_min_edge_cover,min_edge_cover) &&
    func_test(G,old_min_dom_set,min_dom_set)
end

@test quad_test(BuckyBall())
@test quad_test(Paley(17))
@test quad_test(RandomTree(20))
@test quad_test(RandomRegular(20,3))
