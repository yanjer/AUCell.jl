module AUCell

using Random, DelimitedFiles, SparseArrays, HypothesisTests, LinearAlgebra

export pathway_AUC_main,
       reAUCluster_kernel, aucell_kernel,
       filter_expr_matrix,
       read_mtx, read_gmt, read_meta

include("code/read_file.jl")

include("code/process_expr.jl")

include("code/auc_pathway_recluster.jl")

include("code/aucell.jl")

"""

# Examples
## default: reAUCluster mode
```jldoctest
julia> pathway_AUC_main(use_testdata = "yes")
1.632442 seconds (8.44 M allocations: 279.642 MiB, 7.95% gc time, 73.87% compilation time)
[ Info: INFO: The size of expression profile was (36602, 8).
1.779532 seconds (4.95 M allocations: 260.557 MiB, 4.14% gc time, 97.91% compilation time)
[ Info: INFO: The filtered of expression profile size was (7549, 8).
0.000320 seconds (27 allocations: 34.672 KiB)
[ Info: INFO: There are 1 pathways to be analyzed.
0.768511 seconds (1.50 M allocations: 99.943 MiB, 2.46% gc time, 95.16% compilation time)
2×5 Matrix{Any}:
"pathways_name"                     ["cluster1"]                                                                                                       …  ["t"]      ["pvalue"]
"HALLMARK_TNFA_SIGNALING_VIA_NFKB"  Any["AAACCCAAGGGTTAAT-1", "AAACCCAAGAAACCAT-1", "AAACCCAAGCAACAAT-1", "AAACCCAAGCCAGAGT-1", "AAACCCACAGCAGATG-1"]     [4.92654]  [0.00263937]
```
## aucell mode
```jldoctest
julia> pathway_AUC_main(use_testdata = "yes", mode = "aucell")
  1.557316 seconds (8.44 M allocations: 279.659 MiB, 3.27% gc time, 78.85% compilation time)
[ Info: INFO: The size of expression profile was (36602, 8).
  1.771720 seconds (4.95 M allocations: 260.557 MiB, 3.69% gc time, 97.39% compilation time)
[ Info: INFO: The filtered of expression profile size was (7549, 8).
  0.000329 seconds (27 allocations: 34.672 KiB)
[ Info: INFO: There are 1 pathways to be analyzed.
  0.667055 seconds (1.75 M allocations: 87.598 MiB, 3.82% gc time, 99.79% compilation time)
[ Info: INFO: According to the meta information, there are 8 groups of data and each group will be analyzed with the rest of the sample.
  3.153389 seconds (6.62 M allocations: 421.960 MiB, 3.39% gc time, 80.62% compilation time)
2×65 Matrix{Any}:
 "GeneSet"                            "AAACCCAAGAAACCAT-1"   "AAACCCAAGAAACCAT-1"   "AAACCCAAGAAACCAT-1"  …   "AAACCCAGTACGGGAT-1"   "AAACCCAGTACGGGAT-1"   "AAACCCAGTACGGGAT-1"
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB"  0.506962               0.500821               0.515332                  0.512858               0.482078               0.440029
```
"""
function pathway_AUC_main(fn_expr::AbstractString = "matrix.mtx",
                          rn_expr::AbstractString = "features.tsv",
                          cn_expr::AbstractString = "barcodes.tsv",
                       fn_feature::AbstractString = "fn_feature.gmt",
                          fn_meta::AbstractString = "fn_meta.txt";
                      fn_meta_delim::AbstractChar = '\t',
                    fn_meta_group::AbstractString = "group",
                 file_format_expr::AbstractString = "read_mtx", # There are two input modes "read_mtx" and "read_expr_matrix" for the expression profile file format.
                                          T::Type = Int32,
                                 feature_col::Int = 2,
                                 barcode_col::Int = 1,
                          rem_delim::AbstractChar = ' ',
                           feature_threshold::Int = 30, # Include features (genes) detected in at least this many cells
                              cell_threshold::Int = 200, # Include profiles (cells) where at least this many features are detected
              file_format_feature::AbstractString = "read_gmt", # There are two input modes "read_gmt" and "read_gsf" for the file format of the features contained in the pathways.
                   fn_feature_delim::AbstractChar = ' ',
            #  use_HALLMARK_pathway::AbstractString = "no",
                             mode::AbstractString = "reAUCluster",   # "reAUCluster" is subgroups based on pathway activation. "aucell" is an optional mode to calculate AUC based on characteristic genes for two groups. 
                                ncell_pseudo::Int = 0, # ncell_pseudo is the number of pseudobulk combined cells in each group
                         auc_x_threshold::Float64 = 1.0,
                               remove_zeros::Bool = true,
                     use_testdata::AbstractString = "no",
                         work_dir::AbstractString = "./")
    cd(work_dir)
    if use_testdata == "yes"
        fn_expr = joinpath(@__DIR__, "..", "test", "matrix.mtx")
        rn_expr = joinpath(@__DIR__, "..", "test", "features.tsv")
        cn_expr = joinpath(@__DIR__, "..", "test", "barcodes.tsv")
		fn_feature = joinpath(@__DIR__, "..", "test", "fn_feature.gmt")
        fn_meta = joinpath(@__DIR__, "..", "test", "fn_meta.txt")
        feature_threshold = 1
        cell_threshold = 1
    end
    # (use_HALLMARK_pathway == "yes") ? fn_feature = joinpath(@__DIR__, "..", "HALLMARK_pathway", "h_all_v2023_1_Hs_symbols.gmt") : fn_feature
    @time mat, fea, bar = (file_format_expr == "read_mtx") ? read_mtx(fn_expr, rn_expr, cn_expr; T, feature_col, barcode_col) : read_expr_matrix(fn_expr, rn_expr, cn_expr; matrix_delim = rem_delim)
    @info "INFO: The size of expression profile was $(size(mat))."
    @time mat, kf, kb = filter_expr_matrix(mat, feature_threshold, cell_threshold)
    @info "INFO: The filtered of expression profile size was $(size(mat))."
    fea = fea[kf]
    bar = bar[kb]
    @time pathway_name, pathway_genes = (file_format_feature == "read_gmt") ? read_gmt(fn_feature) : read_gsf(fn_feature; delim = fn_feature_delim)
    @info "INFO: There are $(length(pathway_name)) pathways to be analyzed."
    if mode == "reAUCluster"
        # 重分簇的结果，第一列为通路，第二列为重分簇1，第三列为重分簇2，，第四列为配对t检验t值，第五列为配对t检验p值，第六列为各样本的AUC值（顺序按照表达谱样本顺序）
        recluster_bar = reAUCluster_kernel(mat, fea, bar, pathway_genes; np = ncell_pseudo, remove_zeros = remove_zeros)
        # 存储重分簇的结果，第一列为通路，第二列为重分簇1，第三列为重分簇2，第四列为配对t检验t值，第五列为配对t检验p值
        recluster_bar = vcat(hcat("pathways_name",[[["cluster1"]] [["cluster2"]] [["t"]] [["pvalue"]] [["sample1","sample2","sample3"]]]),hcat(pathway_name,recluster_bar))
        writedlm("reAUCluster_result.tsv", recluster_bar[:,1:5], "\t")
        #存储各样本在各通路中的AUC，行为通路，列为样本
        writedlm("reAUCluster_AUC.tsv", recluster_bar[:,[:,setdiff([1:end]...,[2:end-1]...)]], "\t")
        return recluster_bar[:,1:5]
    else
        @time grp, nam = read_meta(fn_meta, fn_meta_group; delim = fn_meta_delim)
        @info "INFO: According to the meta information, there are $(length(grp)) groups of data and each group will be analyzed with the rest of the sample."
        result = []
        cluster = []
        #writedlm("pathway_names.tsv", pathway_name, '\t')
        @time for i in 1:length(nam)
            profiles = mat[:, indexin(grp[i], bar)]
            res = aucell_kernel(mat, fea, pathway_genes, np = ncell_pseudo, x_threshold = auc_x_threshold, remove_zeros = remove_zeros)
            _, c =size(res)
            push!(cluster, fill(nam[i],c))
            #writedlm(nam[i]*"_hallmark.tsv", res, '\t')
            push!(result, res) 
        end
        #/ using BSON: @save
        #/ comb = (result, cluster, pathway_name)
        #/ @save "aucell_result.bson"  comb

        result = hcat(pathway_name, Matrix{Any}(reduce(hcat, result)))
        cluster = vcat(["GeneSet"], reduce(vcat, cluster))
        # Vectory转Matrix
        result = vcat(reshape(cluster, 1, :), result)
        writedlm("aucell_result.tsv", result, '\t')
        return result
    end
end

end
