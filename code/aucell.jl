export aucell_kernel

using Random, DelimitedFiles, SparseArrays

include(joinpath(@__DIR__, "..", "code", "process_expr.jl"))

include(joinpath(@__DIR__, "..", "code", "roc.jl"))

"""
    aucell_kernel(mat, features, gene_set)

Calculate the AUC  for each `gene_set` in each profile of `mat` and row names of `mat` are stored as `features` which should be the same types with those in `gene_set`.

# Examples
```jldoctest
julia> mat = [1 2 3;4 5 6;0 1 2;7 8 0]
4×3 Matrix{Int64}:
 1  2  3
 4  5  6
 0  1  2
 7  8  0

julia> fea = ["a", "b", "c", "d"]

julia> gene_set = ["b","c", "e"]

julia> aucell_kernel(mat, fea, gene_set)
1×3 Matrix{Float64}:
 0.125  0.125  0.5

julia> gene_sets = [["a", "b", "e"], ["b", "d", "e"]]

julia> aucell_kernel(mat, fea, gene_sets)
2×3 Matrix{Float64}:
 0.25  0.25  0.75
 0.75  0.75  0.375
```
"""
function aucell_kernel(
		mat::AbstractMatrix, # Expression profiles matrix 
		features::AbstractVector, # Row names for the profiles
		gene_set::AbstractVector; # A set of genes (features) to be tested or a group of gene-set (Vector of Vectors) 
		np::Int = 0, # If np > 0, pseudobulk mode will be turned on
		x_threshold::Number = 1, # threshold for calculating AUC, (default:1, common AUC)  
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
		nmat = generate_pseudobulk(nmat, np)
	end
	if typeof(first(gene_set)) <: Vector
		return mapreduce(x-> roc_kernel(nmat, fea .∈ (x, ), x_threshold = x_threshold), vcat, gene_set)
	else
		pos = fea .∈ (gene_set,)
		return roc_kernel(nmat, pos, x_threshold = x_threshold)
	end
end

