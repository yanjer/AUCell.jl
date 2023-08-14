export aucell_kernel

using Random, DelimitedFiles, SparseArrays

include(joinpath(@__DIR__, "code", "process_expr.jl"))

include(joinpath(@__DIR__, "code", "roc.jl"))

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


"""
    cell_marker_score(mat, features, barcodes, gene_set, group)

Given a single-cell RNA expression matrix `mat` with row-names of `features` and column-names of `barcodes`, calculate the relative cell type marker scores (0-1) for the ` gene_set`; the grouping information is specified in the `group` (vector of vectors, which store the cell barcodes in each group).

# Examples
```jldoctest

julia> mat = rand(0:32, 12, 8)

julia> features = 1:12

julia> gene_set = [1,5,6,8]

julia> barcodes = ["a", "b", "c", "d", "e", "f", "g", "h"]

julia> group = [["a", "b", "g", "h"], ["c", "d", "e", "f"]]
2-element Vector{Vector{String}}:
 ["a", "b", "g", "h"]
 ["c", "d", "e", "f"]

julia> cell_marker_score(mat, features, barcodes, gene_set, group)
4 genes are found among 4 genes.
1×2 Matrix{Float64}:
 0.476227  0.523773

```
"""
function cell_marker_score(
		mat::AbstractMatrix, # Expression profiles matrix 
		features::AbstractVector, # Row names for the profiles
		bar::AbstractVector, # Row names for the profiles
		gene_set::AbstractVector, # A set of cell-specific genes (expressed only in certain kind of cells, or upregulated)
		group::AbstractVector
	)
	r, c = size(mat)
	r == length(features) ||  throw(DimensionMismatch("`mat` and `features` do not have equal number of rows."))
	c == length(bar)    ||  throw(DimensionMismatch("`mat` and `bar` do not have equal number of columns."))
	gen = unique(gene_set)
	ind = filter(.!=(nothing), indexin(gen, features))
	isnothing(ind) && throw("None of marker genes are found in the features.")
	dat = mat[ind,:]
	println(size(dat, 1), " genes are found among ", length(gen), " genes.")
	if size(dat,1) == 0
		return zeros(Float64, 1, length(group))
	else
		res = mapreduce(i -> mean(dat[:, filter(.!=(nothing),indexin(i, bar))], dims=2), hcat, group)
		res = mapslices(x -> (sum(x) == 0 ? normalize(x .+ 1) : normalize(x, 1)), res, dims = 2)
		return normalize(sum(res, dims = 1), 1)
	end
end
