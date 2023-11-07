export reAUCluster_kernel

using Random, DelimitedFiles, SparseArrays, HypothesisTests, LinearAlgebra

include("roc.jl")

include("process_expr.jl")

function cluster_inter_auc(nmat::AbstractMatrix,
                          group_cell::Vector{Int64},
                          replace_cell::Int64,
                          gs_num::Int64,
                          positives::BitVector)
    alter_cluster_sets = reduce(hcat, [[[group_cell[setdiff(1:end,i)];replace_cell]] for i in 1:gs_num])
    new_group_auc = mapreduce(x-> roc_kernel(reshape(sum.(eachrow(nmat[:, x])),:,1), positives), hcat, alter_cluster_sets)
    return new_group_auc
end

"""
Aside each cell into clusters.
"""
function classify_cell_cluster(nmat::AbstractMatrix,
                      max_group::Vector{Int64},
                      max_group_auc::Float64,
                      min_group::Vector{Int64},
                      min_group_auc::Float64,
                      all_group_sample::Vector{Int64},
                      diff_cell::Vector{Int64},
                      positives::BitVector)
    gs_num = size(max_group)[1]
    group2 = copy(min_group)
    group1 = copy(max_group)
    # # 每列是一个细胞依次替换初始簇后的各个均值
    # println("计算时间")
    # @time max_group_cell = reduce(hcat, [cluster_inter_auc(nmat,max_group,diff_cell[j],gs_num,positives) for j in 1:size(diff_cell)[1]])
    # @time min_group_cell = reduce(hcat, [cluster_inter_auc(nmat,min_group,diff_cell[j],gs_num,positives) for j in 1:size(diff_cell)[1]])

    # d_recell_gmax = sum(max_group_cell, dims = 2)/gs_num .- max_group_auc
    # d_recell_gmin = sum(min_group_cell, dims = 2)/gs_num .- min_group_auc

    # group1_local = (((d_recell_gmax .>= 0) .&& (d_recell_gmin .>= 0)) .|| (.!((d_recell_gmax .< 0) .&& (d_recell_gmin .< 0)) .&& (abs.(d_recell_gmax) .< abs.(d_recell_gmin))))

    # @time max_group_cell = mapreduce(x-> cluster_inter_auc(nmat,max_group,x,gs_num,positives), hcat, diff_cell)



    for j in 1:size(diff_cell)[1]
        max_group_cell = cluster_inter_auc(nmat,max_group,diff_cell[j],gs_num,positives)
        min_group_cell = cluster_inter_auc(nmat,min_group,diff_cell[j],gs_num,positives)
        # 计算替换初始簇细胞后的平均AUC
        d_recell_gmax = sum(max_group_cell)/gs_num - max_group_auc
        d_recell_gmin = sum(min_group_cell)/gs_num - min_group_auc
        # 当都大于0则划分到max簇；当都小于0则划分到min簇；其余直接比较大小，分类到距离小的
        group1_local = (((d_recell_gmax .>= 0) .&& (d_recell_gmin .>= 0)) .|| (.!((d_recell_gmax .< 0) .&& (d_recell_gmin .< 0)) .&& (abs.(d_recell_gmax) .< abs.(d_recell_gmin))))
        if group1_local
            group1 = vcat(group1,diff_cell[j,:])
            max_auc = findmax(max_group_cell)
            if(max_auc[1] > max_group_auc)
                max_group[max_auc[2][1]] = diff_cell[j]
                max_group_auc = max_auc[1]
            end
        else
            group2 = vcat(group2,diff_cell[j,:])
            min_auc = findmin(min_group_cell)
            if (min_auc[1] < min_group_auc)
                min_group[min_auc[2][1]] = diff_cell[j]
                min_group_auc = min_auc[1]
            end
        end
    end
    return group1,group2
end

"""
    pathway_cluster mode. `pathway_cluster` is subgroups based on pathway activation.


# Examples
```jldoctest
julia> mode_pathway_cluster([[1,5,7] [6,4,3] [8,5,2]],[[1,5,7] [6,4,3] [8,5,2]],BitVector([0,1,1]),reshape([1,2,3],:,1),["sample1","sample2","sample3"])
  0.005570 seconds (2.62 k allocations: 150.455 KiB, 97.78% compilation time)
1×5 Matrix{Any}:
 
["sample1"]  ["sample2", "sample3"]  NaN  NaN  [0.5 0.0 0.0]
```
 """
function mode_pathway_cluster(nmat::AbstractMatrix,
                            nmat_p::AbstractMatrix,
                            positives::BitVector,
                            it_group::AbstractMatrix,
                            barcodes::AbstractVector;  
                            decreasing::Bool = true,
                              auc_only::Bool = true,# must be true
                               verbose::Bool = false)
    expr_auc = roc_kernel(nmat_p, positives)
    g_max = findmax(expr_auc)
    g_min = findmin(expr_auc)
    (g_max != g_min) || throw(ArgumentError("The AUC of all samples in this pathway remained consistent at $(g_max[1])"))
    max_group = it_group[g_max[2][2],:]
    min_group = it_group[g_min[2][2],:]
    all_group_sample = vcat(max_group,min_group)
    diff_cell = setdiff([1:size(nmat)[2]...],all_group_sample)
    # group1为AUC最大的初始簇获得，group2由AUC最小的初始簇获得
    @time group1,group2 = classify_cell_cluster(nmat,max_group,g_max[1],min_group,g_min[1],all_group_sample,diff_cell,positives)
    # 返回分簇的样本结果，第一列为重分簇1，第二列为重分簇2，第三列为各样本的AUC值（顺序按照表达谱样本顺序）。
    if (size(nmat) == size(nmat_p))
        ttest_result = EqualVarianceTTest(expr_auc[:,group1][:],expr_auc[:,group2][:])
        return [[barcodes[group1]] [barcodes[group2]] [[ttest_result.t]] [[pvalue(ttest_result)]] [expr_auc]]
    else
        expr_auc_all = roc_kernel(nmat, positives)
        ttest_result = EqualVarianceTTest(expr_auc_all[:,group1][:],expr_auc_all[:,group2][:])
        return [[barcodes[group1]] [barcodes[group2]] [[ttest_result.t]] [[pvalue(ttest_result)]] [expr_auc_all]]
    end
end

function reAUCluster_kernel(
		mat::AbstractMatrix, # Expression profiles matrix 
		features::AbstractVector, # Row names for the profiles
        barcodes::AbstractVector, # Col names for the profiles
		gene_set::AbstractVector; # A set of genes (features) to be tested or a group of gene-set (Vector of Vectors)
		np::Int = 0, # If np > 0, pseudobulk mode will be turned on
		remove_zeros::Bool = true # Whether clear those features with 0 expression values in all profiles
	)
	r, c = size(mat)
	r == length(features) ||  throw(DimensionMismatch("'mat' and `features' do not have equal number of rows."))
	fea = features
	nmat = mat
	if remove_zeros
		keep = reshape(sum(mat, dims = 2), :) .> 0
		nmat = nmat[keep,:]
		fea = fea[keep]
	end
	if np > 0
		nmat_p,it_group = generate_pseudobulk(nmat, np)
    else
        it_group = reshape([1:size(nmat)[2]...],:,1)
        nmat_p = nmat
    end
	if typeof(first(gene_set)) <: Vector
		return mapreduce(x-> mode_pathway_cluster(nmat, nmat_p, fea .∈ (x, ), it_group, barcodes), vcat, gene_set)
	else
		pos = fea .∈ (gene_set,)
		return mode_pathway_cluster(nmat, nmat_p, fea .∈ (gene_set, ), it_group, barcodes)
	end
end
