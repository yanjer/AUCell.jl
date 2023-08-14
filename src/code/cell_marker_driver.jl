include("/public/xwang/doc/aucell/auc/src/aucell.jl")

@time grp, nam = read_meta("meta.tsv", "group")
#@time pathway_name, pathway_genes = read_gmt("/public/xwang/doc/ref/human_motor_cortex_subclass.gmt")
@time pathway_name, pathway_genes = read_gmt("/public/xwang/doc/ref/pbmc_celltype_l2.gmt")

@time mat, fea, bar = read_mtx("matrix.mtx", "features.tsv", "barcodes.tsv")
@time mat, kf, kb = filter_expr_matrix(mat)
fea = fea[kf]
bar = bar[kb]

result = []
@time for i in 1:length(pathway_name)
	gene_set = pathway_genes[i]
	score = cell_marker_score(mat, fea, bar, gene_set, grp)
	push!(result, score) 
end

result = hcat(pathway_name, reduce(vcat, result))
cluster = vcat(["GeneSet"], nam)
result = vcat(reshape(cluster, 1, :), result)
#writedlm("human_motor_cortex_subclass_score.tsv", result, '\t')
writedlm("pbmc_celltype_l2_score.tsv", result, '\t')

