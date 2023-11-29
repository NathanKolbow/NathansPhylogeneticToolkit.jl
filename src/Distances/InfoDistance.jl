# This code was originally written for the `R` package `TreeDist` by authors
# Martin R. Smith, Roy Jonker, Yong Yang, and Yi Cao. Code is copied more or
# less one-to-one, so all credit goes to these original authors.

# TODO: organize bipartite split methods out of this file
using PhyloNewtorks, Assignment


global const __LOG2::Vector{Float64} = Vector{Float64}([log2(j) for j=1:1e6])     # using this instead of calling log saves a few milliseconds

# TODO: implement tree-like checks
function steelDistance(tree1::HybridNetwork, tree2::HybridNetwork)
    s1, leafmap = @inbounds getSplitBits(tree1)
    s2 = @inbounds getSplitBits(tree2, leafmap=leafmap, return_leafmap=false)

    return clusteringInfoDistance(s1, s2, normalize=true)
end


function steelDistance(s1::BitMatrix, s2::BitMatrix)
    @inbounds return clusteringInfoDistance(s1, s2, normalize=true)
end


function clusteringInfoDistance(s1::BitMatrix, s2::BitMatrix; normalize=false)
    mci = calculateMutualClusteringInfo(s1, s2)[1]
    info1 = clusteringEntropy(s1)
    info2 = clusteringEntropy(s2)
    
    indepInfo = info1 + info2
    ret = indepInfo - mci - mci
    ret = normalize ? ret / indepInfo : ret

    if ret > -1e-10 && ret < 1e-10
        return 0
    end
    return ret
end


@inbounds function calculateMutualClusteringInfo(a::BitMatrix, b::BitMatrix)
    max_score = 1e8
    
    nTip::Int64 = size(a)[2]
    if nTip > 1000
        error("More than 1000 tips is not supported by default. Contact the maintainer if you need to accomodate more than 1000 tips.")
    end
    asplits::Int64 = size(a)[1]
    bsplits::Int64 = size(b)[1]

    most_splits::Int64 = max(asplits, bsplits)
    a_extra_splits::Int64 = max(asplits-bsplits, 0)
    b_extra_splits::Int64 = max(bsplits-asplits, 0)

    if most_splits == 0 || nTip == 0
        return 0
    end

    scores::Array{Float64} = Array{Float64}(undef, most_splits, most_splits)
    scores .= 0

    exact_match_score = 0.
    exact_matches = 0

    a_match::Array{Int64} = Array{Int64}(undef, asplits)
    a_match .= 0
    b_match::Array{Int64} = Array{Int64}(undef, bsplits)
    b_match .= 0

    nas::Array{Int64} = Array{Int64}(undef, asplits)
    nas .= 0
    nbs::Array{Int64} = Array{Int64}(undef, bsplits)
    nbs .= 0
    a_and_bs::Array{Int64} = Array{Int64}(undef, asplits, bsplits)
    a_and_bs .= 0

    for bit_idx=1:nTip
        for ai=1:asplits
            if !a[ai, bit_idx] continue end

            for bi=1:bsplits
                if b[bi, bit_idx]
                    a_and_bs[ai, bi] += 1
                end
            end
        end
    end

    for ai=1:asplits
        nas[ai] = 0
        for bit_idx=1:nTip
            if a[ai,bit_idx]
                nas[ai] += 1
            end
        end
    end
    for bi=1:bsplits
        nbs[bi] = 0
        for bit_idx=1:nTip
            if b[bi,bit_idx]
                nbs[bi] += 1
            end
        end
    end

    for ai = 1:asplits
        if a_match[ai] != 0 continue end

        leavesinA = nas[ai]
        for bi = 1:bsplits
            if b_match[bi] != 0 continue end

            leavesinB = nbs[bi]
            a_and_b::Int64 = a_and_bs[ai, bi]

            if (leavesinA == a_and_b && leavesinB == a_and_b) ||
                (a_and_b == 0 && nTip-leavesinA == leavesinB)

                exact_match_score += ic_matching(leavesinA, nTip-leavesinA, nTip)
                exact_matches += 1
                a_match[ai] = bi
                b_match[bi] = ai
                break
            elseif (leavesinB == leavesinA == a_and_b + a_and_b && nTip == leavesinA + leavesinB)
                scores[ai, bi] = max_score
            else
                ic_sum = ic_element(a_and_b, leavesinA, leavesinB, nTip) +
                        ic_element(leavesinA - a_and_b, leavesinA, nTip-leavesinB, nTip) +
                        ic_element(leavesinB - a_and_b, nTip-leavesinA, leavesinB, nTip) +
                        ic_element(nTip-leavesinA-leavesinB+a_and_b, nTip-leavesinA, nTip-leavesinB, nTip)
                @assert ic_sum >= 0 && ic_sum <= nTip

                scores[ai, bi] = max_score - (max_score * (ic_sum / nTip))
            end
        end

        for bi=(bsplits+1):most_splits
            scores[ai, bi] = max_score
        end
    end

    if exact_matches == bsplits || exact_matches == asplits
        return exact_match_score / nTip
    end

    lap_dim = most_splits - exact_matches

    if exact_matches != 0
        nomatch_scores = scores[[j for j=1:most_splits if j ∉ b_match], [j for j=1:most_splits if j ∉ a_match]]

        score = ((max_score * lap_dim) - find_best_assignment(nomatch_scores).cost) / max_score
        score += exact_match_score / nTip
        return score
    else
        for ai=(asplits+1):most_splits
            for bi=1:most_splits
                scores[ai, bi] = max_score
            end
        end

        score = ((max_score*lap_dim) - find_best_assignment(scores).cost) / max_score
        return score
    end
end


function ic_element(nkK::Integer, nk::Integer, nK::Integer, n::Integer)
    global __LOG2
    if nkK == 0 || nk == 0 || nK == 0
        return 0
    end

    if nkK == nk && nkK == nK && nkK + nkK == n return nkK end
    
    num = nkK * n
    denom = nk * nK

    if num == denom return 0 end

    # return nkK * (__LOG2[num] - __LOG2[denom])  below is faster b/c it saves garbage collection time
    return nkK * (__LOG2[nkK] + __LOG2[n] - __LOG2[nk] - __LOG2[nK])
end


@inline function quickLog2(val::Integer)
    return val == 0 ? 0 : __LOG2[val]
end


@inline function ic_matching(a::Integer, b::Integer, n::Integer)
    return (a+b) * quickLog2(n) - a*quickLog2(a) - b*quickLog2(b)
end


function clusteringEntropy(tree::Union{HybridNetwork, <:AbstractString}; p=1, summation=true)
    typeof(tree) == HybridNetwork || (tree = readTopology(tree))
    inSplitCounts = tipsInSplits(tree)
    ret = [bitentropy([tree.numTaxa - t, t] ./ tree.numTaxa) for t in inSplitCounts]
    return summation ? sum(ret) : ret
end


@inline function clusteringEntropy(bits::BitMatrix)
    total = 0
    for row in eachrow(bits)
        rowsum = sum(row)
        total += bitentropy([size(bits)[2] - rowsum, rowsum] ./ size(bits)[2])
    end
    return total
end


@inline function bitentropy(vec)
    return -sum(vec .* log2.(vec))
end


"""
Gets the properly aligned split bits for each tree in the vector.
"""
Base.@propagate_inbounds function getSplitBits(trees::AbstractVector{HybridNetwork}; return_leafmap::Bool=false)
    if length(trees) == 0
        return nothing
    elseif length(trees) == 1
        return @inbounds getSplitBits(trees[1], return_leafmap=return_leafmap)
    else
        bitmats = Array{BitMatrix}(undef, length(trees))
        leafmap = nothing

        for (i, _) in enumerate(trees)
            if leafmap === nothing
                bits, leafmap = @inbounds getSplitBits(trees[i])
                bitmats[i] = bits
            else
                bitmats[i] = @inbounds getSplitBits(trees[i], leafmap=leafmap, return_leafmap=false)
            end
        end
        return return_leafmap ? (bitmats, leafmap) : bitmats
    end
end


"""
Returns a BitMatrix with all bipartite graphs splits in `tree`,
and, by default, a Dict mapping tip labels to bits.
"""
Base.@propagate_inbounds function getSplitBits(tree::Union{<:AbstractString, HybridNetwork};
    leafmap::Union{Dict{<:AbstractString, <:Integer}, Nothing}=nothing, return_leafmap::Bool=true, unsafe::Bool=false)
    typeof(tree) == HybridNetwork || (tree = readTopology(tree))

    if !unsafe
        for node in tree.node
            if length(node.edge) != 1 && node.number >= 0
                @warn "Re-reading tree in `getSplitBits` because internal node had non-negative number."
                tree = readTopology(writeTopology(tree))
                break
            end
        end
    end

    # Allocate the most bits we could possibly need, then free up 
    # unnecessary space at the end
    ret = BitMatrix(undef, length(tree.edge), tree.numTaxa)
    split_dict = Dict{Int64, Int64}()
    
    leafmap_predefined = !(leafmap === nothing)

    if !leafmap_predefined
        leafmap = Dict{String, Int64}()
    end

    logged_root = false
    bit_idx = 1
    for leaf in tree.leaf
        if !leafmap_predefined
            leafmap[leaf.name] = bit_idx
            split_dict[leaf.edge[1].number] = length(split_dict) + 1

            ret[length(split_dict),:] .= 0
            ret[length(split_dict),bit_idx] = 1

            bit_idx += 1
        else
            split_dict[leaf.edge[1].number] = leafmap[leaf.name]

            ret[leafmap[leaf.name],:] .= 0
            ret[leafmap[leaf.name],leafmap[leaf.name]] = 1
        end

        if leaf.edge[1] == tree.node[tree.root].edge[1] || leaf.edge[1] == tree.node[tree.root].edge[2]
            logged_root = true
        end
    end

    ret_idx = length(leafmap) + 1
    for edge in tree.edge
        if edge.node[1].number > 0
            # edge leads to a leaf, so it can't be a bipartite split
            continue
        end

        # only log a root edge once
        if edge == tree.node[tree.root].edge[1] || edge == tree.node[tree.root].edge[2]
            if !logged_root
                logged_root = true
            else
                continue
            end
        end

        # This function fills in the provided BitMatrix
        splitBitHelper(edge, split_dict, ret)
        ret_idx += 1
    end

    if return_leafmap
        return (ret[(length(leafmap)+1):(ret_idx-1), :], leafmap)
    else
        return ret[(length(leafmap)+1):(ret_idx-1), :]
    end
end


function splitBitHelper(edge::PhyloNetworks.EdgeT{PhyloNetworks.Node}, split_dict::Dict{Int64, Int64}, data::BitMatrix)
    if edge.number in keys(split_dict)
        return data[split_dict[edge.number],:]
    end

    leftleaves = edge.node[1].edge[1] in keys(split_dict) ? data[split_dict[edge.node[1].edge[1].number],:] : splitBitHelper(edge.node[1].edge[1], split_dict, data)
    rightleaves = edge.node[1].edge[2] in keys(split_dict) ? data[split_dict[edge.node[1].edge[2].number],:] : splitBitHelper(edge.node[1].edge[2], split_dict, data)
    
    split_dict[edge.number] = length(split_dict) + 1
    data[length(split_dict),:] = leftleaves .|| rightleaves
end


@inline function tipsInSplits(bits::BitMatrix)
    return [sum(col) for col in eachcol(bits)]
end


"""

"""
function tipsInSplits(tree::Union{<:AbstractString, HybridNetwork})
    typeof(tree) == HybridNetwork || (tree = readTopology(tree))

    ret = []
    d = Dict()

    for leaf in tree.leaf
        d[leaf.edge[1]] = 1
    end

    for edge in tree.edge
        if edge == tree.node[tree.root].edge[1]
            continue
        end
        push!(ret, tipsHelper(edge, d))
    end

    return ret[ret .!= 1 .&& ret .!= tree.numTaxa-1]
end


function tipsHelper(edge::PhyloNetworks.EdgeT{PhyloNetworks.Node}, d::Dict)
    if edge in keys(d)
        return d[edge]
    end

    leftcount = edge.node[1].edge[1] in keys(d) ? d[edge.node[1].edge[1]] : tipsHelper(edge.node[1].edge[1], d)
    rightcount = edge.node[1].edge[2] in keys(d) ? d[edge.node[1].edge[2]] : tipsHelper(edge.node[1].edge[2], d)

    return d[edge] = leftcount + rightcount
end


function getPairScores(a::BitMatrix, b::BitMatrix)
    max_splits = max(size(a)[1], size(b)[1])
    mat = Array{Float64}(undef, max_splits, max_splits)
    for i=1:max_splits
        for j=1:max_splits
            mat[i, j] = calculateMutualClusteringInfo(a[[i],:], b[[j],:])
        end
    end

    return mat    
end