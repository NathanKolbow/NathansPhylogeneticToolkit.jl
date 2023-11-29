# Author: Nathan Kolbow (nkolbow@wisc.edu)
# GitHub: https://github.com/NathanKolbow
#
# Source code adapated from https://github.com/ms609/TreeDist


using PhyloNetworks, LinearAlgebra

@inline function KendallColijn(vector1::AbstractVector{<:Real}, vector2::AbstractVector{<:Real})
    return norm(vector1 .- vector2)
end

function KendallColijn(tree1::Union{AbstractString, HybridNetwork}, tree2::Union{AbstractString, HybridNetwork}, Vector; unsafe::Bool=true)
    typeof(tree1) == HybridNetwork || (tree1 = readTopology(tree1))
    typeof(tree2) == HybridNetwork || (tree2 = readTopology(tree2))
    
    FunValue(nTip) = nTip * (nTip-1) / 2
    
    if !unsafe
        if length(tipLabels(tree1)) != length(tipLabels(tree2)) || 
            length(setdiff(tipLabels(tree1), tipLabels(tree2))) > 0 ||
            length(setdiff(tipLabels(tree2), tipLabels(tree1))) > 0

            error("Leaves must have identical labels.")
        end
    end
    return norm(Vector(tree1) .- Vector(tree2))
end


function PathVector(tree::HybridNetwork)
    edges = tree.edge
    nTip = tree.numTaxa
    nNode = tree.numNodes
    tipOrder = order(tipLabels(tree))

    tip_idxs = Array{Int64}(undef, nTip)
    for (j, leaf) in enumerate(tree.leaf)
        for (i, node) in enumerate(tree.node)
            if node.number == leaf.number
                tip_idxs[j] = i
                break
            end
        end
    end

    rowcoltoidx(row, col) = (col-1)*(col-2)÷2+row

    ancestors = findancestors(tree)
    pathLength = Array{Int64}(undef, nTip*(nTip-1)÷2)   # nTip choose 2

    for i=1:nTip
        for j=1:nTip
            if tipOrder[i] >= tipOrder[j]
                continue
            end

            anc1 = ancestors[i]
            anc2 = ancestors[j]

            plsum = findpathlengthsum(anc1, anc2)
            # mrca = minimum(intersect(anc1, anc2))
            # @assert plsum == sum(anc1 .< mrca) + sum(anc2 .< mrca)

            pathLength[rowcoltoidx(tipOrder[i], tipOrder[j])] = plsum
        end
    end

    return pathLength
end


function order(arr::AbstractArray)
    sortperm_ret = sortperm(arr)
    ret_arr = Array{Int64}(undef, length(arr))
    for (i, idx) in enumerate(sortperm_ret)
        ret_arr[idx] = i
    end
    return ret_arr
end


function findancestors(tree::HybridNetwork)
    retarr = Array{Vector{Integer}}(undef, tree.numTaxa)
    root_num = tree.node[tree.root].number
    for (i, leaf) in enumerate(tree.leaf)
        retarr[i] = Vector{Integer}()
        node = leaf
        while node.number != root_num
            node = node.edge[length(node.edge)].node[2]
            push!(retarr[i], node.number)
        end
    end
    return retarr
end


"""
Finds the mrca of `anc1` and `anc2`. Saves significant time over
minimum(intersect(anc1, anc2)).
"""
function findpathlengthsum(anc1, anc2)
    idx1 = 1
    idx2 = 1
    lanc1 = length(anc1)
    lanc2 = length(anc2)

    while idx1 <= lanc1 && idx2 <= lanc2
        if @inbounds anc1[idx1] == anc2[idx2]
            return idx1 + idx2 - 2
        elseif @inbounds anc1[idx1] < anc2[idx2]
            idx1 += 1
        else
            idx2 += 1
        end
    end
    
    error("Loop terminated without finding common ancestors.")
end