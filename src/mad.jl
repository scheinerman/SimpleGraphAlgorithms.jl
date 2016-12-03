export ad, mad, mad_core

"""
`mad_model(G)` returns the solved linear program whose optimum
value is `mad(G)`.
"""
function mad_model(G::SimpleGraph)
  EE = elist(G)
  VV = vlist(G)

  n = length(VV)
  m = length(EE)

  MOD = Model()

  # variables
  @variable(MOD, z >= 0)
  @variable(MOD, x[i=1:n,j=1:m; in(VV[i],EE[j])] >= 0 )

  # vertex constraints
  for i=1:n
    v = VV[i]
    # get all edges that have v as an end point
    j_iterator = (j for j=1:m if EE[j][1]==v || EE[j][2]==v)
    j_list = collect(j_iterator)
    @constraint(MOD, sum{x[i,j],j = j_list}<=z)
  end

  # edge constraints
  for j=1:m
    v,w = EE[j]
    i1 = first(find( x -> x==v, VV))
    i2 = first(find( x -> x==w, VV))
    @constraint(MOD, x[i1,j]+x[i2,j] == 2)
  end

  @objective(MOD, :Min, z)
  solve(MOD)
  return MOD, x, VV, EE
end


"""
`mad(G)` computes the maximum average degree of `G`.
"""
function mad(G::SimpleGraph)
  MOD, x, VV, EE = mad_model(G)
  return getobjectivevalue(MOD)
end

"""
`ad(G)` is the average degree of a vertex in `G`.
"""
ad(G::SimpleGraph) = 2*NE(G)/NV(G)

"""
`mad_core(G)` returns a subgraph `H` of `G` for which
`ad(H)==mad(G)`.
"""
function mad_core(G::SimpleGraph)
  T = vertex_type(G)
  GG = deepcopy(G)
  while true
    if NV(GG)==0
      error("Bad graph or something went wrong")
    end
    # solve the LP for GG
    MOD, x, VV, EE = mad_model(GG)
    n = length(VV)
    m = length(EE)
    err = 0.1/n
    mval = getobjectivevalue(MOD)

    # if balanced, return
    if abs(ad(GG)-mval) <= err
      return GG
    end

    # otherwise, find and destroy defective vertices
    for i=1:n
      v = VV[i]
      total = 0.0
      for j=1:m
        e = EE[j]
        if in(v,e)
          val = getvalue(x[i,j])
          total += val
        end
      end
      if total < mval - err   # this is a defective vertex
        delete!(GG,v)
      end
    end
  end # end while
  error("This can't happen (but it did)")
  return GG
end



# For debugging purposes only ...
function madness(G::SimpleGraph)
  MOD, x, VV, EE = mad_model(G)
  opt = getobjectivevalue(MOD)
  println("mad(G) = $opt\n")
  n = length(VV)
  m = length(EE)
  for i=1:n
    v = VV[i]
    println("Vertex $v info")
    total = 0.0
    for j=1:m
      e = EE[j]
      if in(v,e)
        val = getvalue(x[i,j])
        total += val
        println("Value at edge $e is $val")
      end
    end
    println("Total value = $total: $(opt==total) \n")
  end
  nothing
end
