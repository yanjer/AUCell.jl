export roc_kernel


#y-axis:   Sensitivity, TPR: TP/(TP+FN)
#x-axis:  1-Specificity, FPR: FP/(TN+FP)
"""
    roc_kernel(scores, positives)

Calculate the receiver operating curve for `scores` where the positive samples are marked as `true` in `positives`. If the keyword argument `decreasing` is set as `false`, it will assume that lower scores mean higher probability as a positive. The return value is a 3*n matrix, where the first column is '1-specificitiy' (false positive rate), the second column is 'sensitivity' (true positive rate) and the third column is the area under the curve (AUC) to the current point. However, if the keyword argument `auc_only` is set true, only the AUC is returned and `x_threshod` can be set to obtain the area under the curve less than the given 'FPR'.

If `scores` is a matrix, each column is assumed to be the score vector. The AUC will be returned these scores with the same ground truth specified by `positives`.


# Examples
```jldoctest
julia> roc_kernel(rand(10), BitVector(rand(Bool, 10)))
5×3 Matrix{Float64}:
 0.0   0.0       0.0
 0.25  0.166667  0.0208333
 0.5   0.166667  0.0625
 0.75  0.166667  0.104167
 1.0   0.5       0.1875

julia> roc_kernel(rand(10), BitVector(rand(Bool, 10)), decreasing = false)
5×3 Matrix{Float64}:
 0.0   0.0       0.0
 0.25  0.0       0.0
 0.5   0.166667  0.0208333
 0.75  0.333333  0.0833333
 1.0   1.0       0.25

julia> roc_kernel(rand(10, 5), BitVector(rand(Bool, 10)))
1×5 Matrix{Float64}:
 0.190476  0.309524  0.357143  0.642857  0.357143

julia> roc_kernel(rand(10, 5), BitVector(rand(Bool, 10)), x_threshold = 0.5)
1×5 Matrix{Float64}:
 0.666667  0.25  0.583333  0.375  0.416667

julia> roc_kernel(rand(10, 5), BitVector(rand(Bool, 10)), x_threshold = 0.5, verbose = true)
x=	1.0	y=	0.25
x=	1.0	y=	0.875
x=	1.0	y=	1.0
x=	1.0	y=	0.75
x=	1.0	y=	1.0
1×5 Matrix{Float64}:
 0.125  0.46875  0.625  0.5  0.3125
```

# Arguments
- `scores::AbstractVector` or `AbstractMatrix`: the scores vector or matrix.
- `positives::BitVector`: the ground-truth vector.
- `decreasing::Bool = true`: score's direction.
- `auc_only::Bool = false`: whether to run the calculation in the AUC-only mode.
- `verbose::Bool = false`: the verbosity of output.

"""
function roc_kernel(
					scores::AbstractVector, # Scores for each sample
			     positives::BitVector; # Whether a sample is a positive, (1 true, 0 false) 
		        decreasing::Bool = true,  # By default, higher scores mean more likely to be a positive; set it to `false` if lower scores mean higher probability
				  auc_only::Bool = false, # If `true`, return `auc` only
			   x_threshold::Number = 1, # return the Area under the curve to this FPR
				   verbose::Bool = false # whether output extra information
			)
	n = length(scores)
	n == length(positives) ||  throw(DimensionMismatch("'scores' and 'positives' do not have equal number of rows."))
	m = sum(positives)   # total number of true positives
	o = sortperm(scores, rev = decreasing)
	tp = cumsum(positives[o]) # True positives
	# TP + FP = 1:n
	fp = (1:n) .- tp  # False positives
	# TP + FN = m 
	#/ fn = m .- fp 
	# FN + TN = n .- (1:n)
	#/ tn = (n - m) .- tp
	y  = vcat([0], (1/m) .* tp)
	x  = vcat([0], (1/(n-m)) .* fp)
	ind= vcat([true], (x[1:n] .!= x[2:n+1]))
	x = x[ind]
	y = y[ind]
	a  = vcat([0], 0.5 * (diff(x) .* (y[1:end-1] .+ y[2:end]))) # Area under each interval
	auc = cumsum(a)
	if auc_only
		ind = findfirst(x -> x >= x_threshold, auc)
		if ind == nothing
			verbose && println("x=\t", x[end], "\ty=\t", y[end])
			return auc[end]
		else
			verbose && println("x=\t", x[ind], "\ty=\t", y[ind])
			return auc[ind]
		end
	else
		return hcat(x, y, auc)
	end
end

# Matrix model for scores 
function roc_kernel(
					scores::AbstractMatrix,
			     positives::BitVector;  
		        decreasing::Bool = true,
				  auc_only::Bool = true,# must be true
				x_threshold::Number = 1, # return the Area under the curve to this FPR
				   verbose::Bool = false # whether output extra information
			)
	r, c = size(scores)
	r == length(positives) ||  throw(DimensionMismatch("'scores' and 'positives' do not have equal number of rows."))
	auc_only ||  throw("If `scores` is a matrix, it can only run in the `auc_only` mode.")
	mapslices(x -> roc_kernel(x, positives, decreasing = decreasing, auc_only = true, x_threshold = x_threshold, verbose = verbose), scores, dims = 1)
end

