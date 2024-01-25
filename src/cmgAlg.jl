## define structures
Base.@kwdef struct lPreconditioner{T<:Function}
  f::T
  function lPreconditioner(pfunc)
    f = pfunc
    return new{typeof(f)}(f)
  end
end

Base.@kwdef struct CholT
  ld::SparseMatrixCSC{Float64,Int64}
  ldT::SparseMatrixCSC{Float64,Int64}
  d::Vector{Float64}
  p::Vector{Int64}
  invp::Vector{Int64}
end

Base.@kwdef struct HierarchyLevel
  sd::Bool
  islast::Bool
  iterative::Bool
  A::SparseMatrixCSC{Float64,Int64}
  invD::Vector{Float64}
  cI::Vector{Int64}
  nc::Int64
  n::Int64
  nnz::Int64
  #chol::CholT
  chol::LDLFactorizations.LDLFactorization{Float64,Int64,Int64,Int64}
end

Base.@kwdef struct Hierarchy
  A::SparseMatrixCSC{Float64,Int64}
  invD::Vector{Float64}
  cI::Vector{Int64}
  #chol::CholT
  chol::LDLFactorizations.LDLFactorization{Float64,Int64,Int64,Int64}
end

Base.@kwdef mutable struct LevelAux
  fwd::Bool
  rc::Int32
  repeat::Int32
  #
  sd::Bool
  islast::Bool
  iterative::Bool
  n::Int64
  nc::Int64
end

Base.@kwdef struct Workspace
  x::Vector{Float64}
  b::Vector{Float64}
  tmp::Vector{Float64}
end

## Hierarchy Structures Initialization

function init_LevelAux(H::Vector{HierarchyLevel})
  local n_levels = length(H)
  local LevelAux_ = Vector{LevelAux}(undef, n_levels)

  @inbounds for j = 1:n_levels
    local repeat::Int32 = 0
    if j == 1
      repeat = 1
    elseif j == n_levels
      repeat = 0
    else
      repeat = max(floor(nnz(H[j-1].A) / nnz(H[j].A) - 1), 1)
    end

    LevelAux_[j] = LevelAux(
      fwd = true,
      rc = 1,
      repeat = repeat,
      sd = H[j].sd,
      islast = H[j].islast,
      iterative = H[j].iterative,
      n = H[j].n,
      nc = H[j].nc,
    )
  end

  return LevelAux_
end

function init_Workspace(H::Vector{HierarchyLevel})
  local n_levels = length(H)
  local Workspace_ = Vector{Workspace}(undef, n_levels)

  @inbounds for j = 1:n_levels
    local n_ = size(H[j].A, 1)
    if j == n_levels && !H[j].iterative
      n_ = n_ + 1
    end
    Workspace_[j] = Workspace(x = zeros(n_), b = zeros(n_), tmp = zeros(n_))
  end

  return Workspace_
end


function init_Hierarchy(H::Vector{HierarchyLevel})
  return [Hierarchy(A = h.A, invD = h.invD, cI = h.cI, chol = h.chol) for h in H]
end

##

function cmg_preconditioner_lap(A_lap::SparseMatrixCSC)
  local A_lap_ = validateInput!(A_lap)  # throws if not valid
  cmg_!(A_lap, A_lap_)
end

function cmg_preconditioner_adj(A::SparseMatrixCSC)
  cmg_preconditioner_lap(lap(A))
end

function cmg_!(A::T, A_::T) where {T<:SparseMatrixCSC}
  A_o = A
  flag = Int64(0)
  sd = true
  j = Int64(0)
  h_nnz = Int64(0)
  n = Int64(0)
  H = Vector{HierarchyLevel}()
  local sflag::Bool = true
  local sd::Bool = size(A_, 1) > size(A, 1)

  # build up H
  while true
    n = size(A_, 1)
    dA_ = Array(diag(A_))
    local cI, nc = steiner_group(A_, dA_)
    islast = false
    A = A_ # !
    invD = 1 ./ (2 * dA_) # !
    R = sparse(cI, 1:n, ones(n), nc, n) # added for efficiency

    if nc == 1
      islast = true
      flag = 1
    end

    # check for hierarchy stagnation for potentially bad reasons
    h_nnz = h_nnz + nnz(A_)
    if (nc >= n - 1) || (h_nnz > 5 * nnz(A_o))
      islast = true
      flag = 3 # indicates stagnation
      @warn "CMG convergence may be slow due to matrix density. Future versions of CMG will eliminate this problem."
      break
    end

    Rt = sparse(cI, 1:n, 1, nc, n) # ! take out double
    A_ = (Rt * A) * Rt'
    push!(
      H,
      HierarchyLevel(
        sd = sd,
        islast = islast,
        iterative = true,
        A = A,
        invD = invD,
        cI = cI,
        nc = nc,
        n = n,
        nnz = nnz(A),
        chol = ldl([1.0 0; 0 1.0]),
      ),
    )

    if sflag
      sd = true
      sflag = false
    end

    if nc == 1
      break
    end
  end

  # code for last hierarchy level
  if flag == 0
    j += 1
    B = A_[1:n-1, 1:n-1]
    ldlt = ldl(B)
    push!(
      H,
      HierarchyLevel(
        sd = true,
        islast = true,
        iterative = false,
        A = B,
        invD = Float64[],
        cI = Int64[],
        nc = 0,
        n = n,
        nnz = nnz(ldlt.L),
        chol = ldlt,
      ),
    )
  end

  X = init_LevelAux(H)
  W = init_Workspace(H)
  M = init_Hierarchy(H)


  # create precondition function
  local pfunc = make_preconditioner(M, W, X)

  return (pfunc, H)
end

function findRowSumAndDominance(
  A::SparseMatrixCSC,
)::Union{Nothing,Tuple{Vector{Float64},Vector{Int64},Vector{Int64},Vector{Float64}}}
  local colptr::Vector{Int64} = A.colptr
  local nzval::Vector = A.nzval
  local rowval::Vector = A.rowval
  local il = similar(nzval, Int64)
  local jl = similar(nzval, Int64)
  local vl = similar(nzval, Float64)
  local ncols = size(A, 2)
  local sumR = zeros(Float64, ncols)

  s = Int64(1)
  k = Int64(0)
  @inbounds for i = 1:ncols
    l = colptr[i+1] - colptr[i]
    t = s + l - 1
    sum = 0
    for j = s:t
      val = nzval[j]
      sum += val
      row = rowval[j]
      k += 1
      il[k] = row
      jl[k] = i
      vl[k] = val
      if row != i && val > 0
        return nothing
      end
    end
    sumR[i] = sum
    s += l
  end

  return sumR, il, jl, vl
end

function validateInput!(A::SparseMatrixCSC)::SparseMatrixCSC
  # check symmetry
  if !issymmetric(A)
    throw(ArgumentError("Input Matrix Must Be Symmetric!"))
    return A
  end
  # detect strict dominance && positive off diagonals
  local n = size(A, 1)
  local sAp = Vector{Float64}(undef, n)
  local sd = Vector{Int64}(undef, n)
  local dA = diag(A)

  local res = findRowSumAndDominance(A)
  if isnothing(res)
    throw(ArgumentError("Current Version of CMG Does Not Support Positive Off-Diagonals!"))
    return A
  end
  local sA, i, j, v = res

  @inbounds @simd for i = 1:length(sA)
    sAp[i] = (sA[i] + abs(sA[i])) / 2
    sd[i] = (sAp[i] / dA[i]) > 1e-13 ? 1 : 0
  end

  # augment by extra coordinate if strictly dominant
  if maximum(sd) > 0.0
    local ex_v = -sAp[sd]
    local ex_v_sum = -sum(ex_v)
    local exd = length(ex_v)
    local exd_c = findall(!iszero, sd) #nonzeros(sd) # get coordinates
    local i_ = ones(Int64, exd) * (n + 1)
    local i = vcat(i, i_, exd_c, n + 1)
    local j = vcat(j, exd_c, i_, n + 1)
    local v = vcat(v, ex_v, ex_v, ex_v_sum)
    A = sparse(i, j, v, n + 1, n + 1)
  end
  return A
end

"""
    function(cI, nc) = steiner_group(A, dA_)
    Steiner groups
    # Arguments
    -`A`: Laplacian
"""

function steiner_group(A::SparseMatrixCSC, dA_::Vector{Float64})
  local C, M = findMinSparse(A)
  split_forest_!(C)
  local efd = abs.(M ./ dA_)
  if minimum(efd) < 1.0 / 8 # low effective degree nodes found
    # TODO(pratyai): Had to disable this because it somehow causes an index-out-of-bounds error.
    # C = update_groups_(A, C, dA_)
  end
  #return C, efd
  local cI, nc, _ = forest_components_(C)
  return cI, nc
end

"""
    function C1 = update_groups_(A, C, dA_)
    update groups based on nodes with low effective degree
    # Arguments
    -`A`: Laplacian
"""
function update_groups_(A::SparseMatrixCSC, C::Vector{Int64}, dA_::Vector{Float64})
  n = length(C)
  B = zeros(Float64, n)
  # B[j] is the total tree weight incident to node j
  @inbounds for i = 1:n
    if C[i] != i
      B[i] = A[i, C[i]] + B[i]
      B[C[i]] = A[i, C[i]] + B[C[i]]
    end
  end
  ndx = findall(x -> x > -0.125, B ./ dA_)
  @inbounds @simd for i = 1:length(ndx)
    C[ndx[i]] = Int32(ndx[i])
  end
  return C
end

"""
    function[cI, nc, csizes] = forest_components_(C)
    forest components, connected components in unimodal forest
    # Arguments
    -`C`: unimodal tree
"""

function forest_components_(C::Vector{Int64})
  local n = length(C)
  local cI = zeros(Int64, n)  # the connected component ID for a node
  local cSizes = zeros(Int64, n)  # the size of a connnected component for an ID
  local ccI = 1  # the next available connected component ID

  local buffer = zeros(Int64, n)
  @inbounds for j = 1:n
    local bufferI = 1
    local jwalk = j
    # tree walk the path toward ancestors
    while cI[jwalk] == 0  # until we have found an already labeled node (or have reached a root)
      cI[jwalk] = ccI  # tentatively label the current node with the next available connected component ID
      buffer[bufferI] = jwalk  # push the current to the ancestry queue
      bufferI = bufferI + 1
      jwalk = C[jwalk]  # move to the parent
    end
    bufferI = bufferI - 1  # and we want an inclusive range
    local en = C[jwalk] # end node

    if cI[en] != ccI  # if our tentative labeling was wrong
      cI[buffer[1:bufferI]] .= cI[en]  # then fix those mislabelings
    else
      ccI = ccI + 1  # otherwise, make a new connected component ID available
    end
    cSizes[en] = cSizes[en] + bufferI  # update the subtree size of that terminating node.
  end

  ccI -= 1  #the current value was the next unused one
  cSizes = cSizes[1:ccI]
  return cI, ccI, cSizes
end

"""
   function pfun = make_preconditioner(H)
Make preconditioner
"""
function make_preconditioner(
  H::Vector{Hierarchy},
  W::Vector{Workspace},
  X::Vector{LevelAux},
)::Function
  local fi = b -> preconditioner_i(H, W, X, b)
  local fsd = b -> preconditioner_sd(b, H, X, W)
  return X[1].sd ? fsd : fi
end

"""
   function x = preconditioner_sd(b,H)
preconditioner sd
"""
function preconditioner_sd(
  b::Array{Float64},
  H::Array{Hierarchy},
  X::Vector{LevelAux},
  W::Vector{Workspace},
)
  local n = length(b)
  local bt = [b; -sum(b)]
  local x = preconditioner_i(H, W, X, bt)
  return x[1:n] .- x[n+1]
end

"""
   function C = split_forest_(C1::Array{Int})
decompose unimodal forest into low conductance components
"""

function split_forest_!(C1::Vector{Int64})
  n = length(C1)
  C = C1
  new_front = Int64(0)
  removed_ancestors = Int64(0)
  k = Int64(0)
  jwalk = Int64(0)
  jwalka = Int64(0)
  jwalkb = Int64(0)
  cut_mode = false
  startwalk = false
  ancestors_in_path = Int64(0)
  ancestors = zeros(Int64, n)
  indegree = zeros(Int64, n + 2)
  visited = falses(n) # logical sparse array
  walkbuffer = zeros(Int64, 20)
  newancestorbuff = zeros(Int64, 20)

  # compute indegrees
  for i = 1:n
    indegree[C[i]] = indegree[C[i]] + 1
  end

  # partition into clusters of small diameter

  for j = 1:n
    jwalk = j
    startwalk = true

    while (startwalk && (indegree[jwalk] == 0) && !visited[jwalk])
      startwalk = false
      ancestors_in_path = 0 # change over C-CMG
      k = 1
      walkbuffer[k] = jwalk
      newancestorbuff[k] = 0
      while (k <= 6 || visited[jwalk]) # +1 for c-indexing adjust
        jwalk = C[jwalk]
        walkterminated = (jwalk == walkbuffer[k]) || ((k > 1) && (jwalk == walkbuffer[k-1]))
        if walkterminated
          break # while
        end
        k += 1
        walkbuffer[k] = jwalk
        if visited[jwalk]
          newancestorbuff[k] = ancestors_in_path
        else
          ancestors_in_path = ancestors_in_path + 1
          newancestorbuff[k] = ancestors_in_path
        end
      end

      if k > 6 # large diameter - cut
        middlek = Int64(ceil(k / 2))
        C[walkbuffer[middlek]] = walkbuffer[middlek] # cut middle edge
        indegree[walkbuffer[middlek+1]] = indegree[walkbuffer[middlek+1]] - 1 # update indegree

        for ik = (middlek+1):k
          ancestors[walkbuffer[ik]] =
            ancestors[walkbuffer[ik]] - ancestors[walkbuffer[middlek]]
        end

        # update ancestors and visited flag
        for ik = 1:middlek
          visited[walkbuffer[ik]] = true
          ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]] + newancestorbuff[ik]
        end

        # set first vertex in new walk
        jwalk = walkbuffer[middlek+1]
        startwalk = true
      end # end cut procedure

      # commit walk changes
      if !startwalk
        for ik = 1:k
          ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]] + newancestorbuff[ik]
          visited[walkbuffer[ik]] = true
        end
      end
    end # outer while
  end

  # tree partition into clusters of high conductance
  for j = 1:n
    jwalk = j
    startwalk = true

    while startwalk && (indegree[jwalk] == 0)
      startwalk = false
      jwalkb = jwalk
      cut_mode = false
      # initialize new_front
      new_front = 0
      removed_ancestors = 0

      while true
        jwalka = C[jwalk]
        walkterminated = (jwalka == jwalk) || (jwalka == jwalkb)
        if walkterminated
          break # while
        end

        if (
          !cut_mode && (ancestors[jwalk] > 2) && (ancestors[jwalka] - ancestors[jwalk] > 2)
        ) # possibly low conductance - make cut
          C[jwalk] = jwalk # cut edge
          indegree[jwalka] = indegree[jwalka] - 1
          removed_ancestors = ancestors[jwalk]
          new_front = jwalka
          cut_mode = true
        end # end making cut

        jwalkb = jwalk
        jwalk = jwalka
        if cut_mode
          ancestors[jwalk] = ancestors[jwalk] - removed_ancestors
        end
      end
      if cut_mode
        startwalk = true
        jwalk = new_front
      end
    end
  end
end # split_forest_

function findMinSparse(A::SparseMatrixCSC)
  local colptr = A.colptr
  local rowval = A.rowval
  local nzval = A.nzval

  local m = size(A, 2)
  local min_cols = zeros(Int64, m)
  local min_vals = zeros(Float64, m)

  for i = 1:m
    local rb, re = colptr[i], colptr[i+1]  # range: [rb, re) in `rowval`
    if rb == re  # empty range
      continue
    end

    local minval::Float64, row_::Int64 = Inf, 0
    for j = rb:(re-1)
      local val, row = nzval[j], rowval[j]
      if val < minval
        minval, row_ = val, row
      end
    end
    min_cols[i], min_vals[i] = row_, minval
  end

  return min_cols, min_vals
end


@inline function interpolate!(x::Vector{Float64}, cI::Vector{Int64}, z::Vector{Float64})
  x .= 0.0
  @inbounds @simd for i = 1:length(z)
    x[cI[i]] += z[i]
  end
end

function preconditioner_i(
  H::Vector{Hierarchy},
  W::Vector{Workspace},
  X::Vector{LevelAux},
  b::Vector{Float64},
)
  level = Int64(1)
  #@inbounds W[1].b = b
  BLAS.blascopy!(length(b), b, 1, W[1].b, 1)

  while level > 0
    x = W[level].x
    b = W[level].b
    tmp = W[level].tmp

    invD = H[level].invD
    A = H[level].A
    cI = H[level].cI

    if X[level].islast && !X[level].iterative
      @inbounds ldiv!(x, H[level].chol, b)
      level -= 1
    elseif X[level].islast && X[level].iterative
      W[level].x .= b .* invD
      level -= 1
    elseif !X[level].islast && X[level].fwd
      repeat = X[level].repeat   #number of 'recursive' calls

      if X[level].rc > repeat
        X[level].rc = 1
        level -= 1
      else
        if X[level].rc == 1
          W[level].x .= b .* invD
        else
          mul!(tmp, A, x)
          tmp .-= b
          tmp .*= -invD
          W[level].x .+= tmp
        end

        mul!(tmp, A, x)
        tmp .-= b
        X[level].fwd = false
        interpolate!(W[level+1].b, cI, -tmp)
        level += 1
      end
    elseif !X[level].islast && !X[level].fwd
      z = W[level+1].x
      W[level].x .+= z[cI]
      mul!(tmp, A, x)
      tmp .-= b
      tmp .*= -invD
      W[level].x .+= tmp
      X[level].rc += 1
      X[level].fwd = true
    end
  end

  return W[1].x
end
