export mad
"""
`mad(G)` computes the maximum average degree of `G`.
"""
function mad(G::SimpleGraph)
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
  # println(getvalue(x))

  return getobjectivevalue(MOD)
end


# @variable(m, x[i=1:10,j=1:10; isodd(i+j)] >= 0)
