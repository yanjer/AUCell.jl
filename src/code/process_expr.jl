export filter_expr_matrix, generate_pseudobulk

using Random, DelimitedFiles, SparseArrays, HypothesisTests, LinearAlgebra

"""
	filter_expr_matrix(mat, feature_threshold, cell_threshold)

Filter an expression matrix `mat`, only keep those genes expressed in greater than `feature_threshold` cells and cells expressing greater than `cell_threshold` features.
Return the filtered matrix and the bit vectors for keeping features and cells.

# Examples
```jldoctest

julia> @time mat, fea, bar = read_mtx("matrix.mtx", "features.tsv", "barcodes.tsv")

julia> size(mat)
(36601, 5744)

julia> @time mat2, kf, kb = filter_expr_matrix(mat)
 26.438175 seconds (978.08 k allocations: 1.320 GiB, 0.52% gc time)
(sparse([2, 12, 15, 25, 26, 27, 29, 32, 34, 37  …  21104, 21105, 21106, 21107, 21108, 21109, 21110, 21111, 21113, 21116], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  5728, 5728, 5728, 5728, 5728, 5728, 5728, 5728, 5728, 5728], Int32[1, 1, 5, 1, 4, 1, 1, 1, 1, 1  …  287, 8, 239, 124, 32, 8, 145, 41, 99, 2], 21121, 5728), Bool[0, 0, 0, 1, 0, 0, 1, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 1, 0], Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

julia> size(mat2)
(21121, 5728)

julia> fea2 = fea[kf]; bar2 = bar[kb];

julia> length(fea2)
21121

julia> length(bar2)
5728

```

# Arguments
- `mat::AbstractMatrix`: expression matrix (either dense or sparse).
- `feature_threshold::Int`: the least number of cells that a feature must express in, in order to be kept. Default: 30.
- `cell_threshold::Int`: the least number of genes that a cell must express, in order to be kept. Default: 200.

"""
function filter_expr_matrix(mat::AbstractMatrix, feature_threshold::Int=30, cell_threshold::Int=200)
	feature_threshold > 0 || throw("`feature_threshold` must be a positive integer.")
	cell_threshold > 0 || throw("`cell_threshold` must be a positive integer.")
	local nf, nc
	if typeof(mat) <: SparseMatrixCSC
		nc = mapslices(nnz, mat, dims = 1) # 1xc
		nf = mapslices(nnz, mat, dims = 2) # rx1
	else
		nc =count(==(0), mat, dims = 1) # 1xc
		nf =count(==(0), mat, dims = 2) # rx1
	end
	kf = reshape(nf .> feature_threshold, :)
	kc = reshape(nc .> cell_threshold, :)
	return (mat[kf, kc], kf, kc)
end

"""
    generate_pseudobulk(mat, np)

Generate a matrix of pseudobulk profiles from `mat` which stores single-cell RNA profiles. Each column represents a cell's profile. Each pseudobulk profile is generated from `np` (default: 10) single-cell profiles.


# Examples
```jldoctest
julia> generate_pseudobulk(rand(0:32, 10, 6), 3)
10×2 Matrix{Int64}:
 59  30
 66  34
 37  26
 58  70
 83  86
 15  11
 58  62
 38  62
 62  35
 15  51
```
"""
function generate_pseudobulk(mat::AbstractMatrix, # Each column is a profile
							  np::Int = 10 # Number of profiles in each pseudobulk profile
	)
	r,c = size(mat)
	np > 0  || throw("`np` must be a positive integer.")
	c >= np || throw("There are not enough profiles to generate a pseduobulk profle.")
	ind = randperm(c) # random permutation
	ns  = floor(Int, c/np) # Number of pseudobulk profiles generated
	ind = reshape(ind[1:(ns*np)], np, :)
	#TODO: `mapslices` and `sum` (with dims) cannot be combined together.
	reduce(hcat, [sum(mat[:,i], dims = 2) for i in eachcol(ind)]),ind
end