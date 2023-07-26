# AUCell

APRC is an algorithm for cell reclassification based on AUC values by feature genes in the pathway.

## 1 Used in the Julia language

### 1.1 Installation

The algorithm is implemented in Julia. Version 1.7 or later is recommended. The simpliest way to install is using the `Pkg` facility in Julia. 

```julia
using Pkg
Pkg.add("AUCell")
```

### 1.2 Examples

#### 1.2.1 Quick Start

Run a test job with the input files distributed with the package.

```julia
julia> using AUCell
# Use the default values for the following other parameters. If you need to modify the parameters, add them directly.
julia> result = pathway_AUC_main(use_testdata="yes")
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `pathway_AUC_main` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```julia
julia> result
  1.595452 seconds (8.44 M allocations: 279.644 MiB, 4.83% gc time, 77.52% compilation time)
[ Info: INFO: The size of expression profile was (36602, 8).
  1.945127 seconds (4.95 M allocations: 260.557 MiB, 11.17% gc time, 96.92% compilation time)
[ Info: INFO: The filtered of expression profile size was (7549, 8).
  0.000401 seconds (27 allocations: 34.641 KiB)
[ Info: INFO: There are 1 pathways to be analyzed.
  0.660084 seconds (1.75 M allocations: 87.597 MiB, 3.11% gc time, 99.78% compilation time)
[ Info: INFO: According to the meta information, there are 2 groups of data and each group will be analyzed with the rest of the sample.
  2.731819 seconds (6.61 M allocations: 365.662 MiB, 3.77% gc time, 94.64% compilation time)
2×17 Matrix{Any}:
 "GeneSet"                            "group1"   "group1"   "group1"   "group1"   "group1"  …   "group2"   "group2"   "group2"   "group2"   "group2"   "group2"   "group2"   "group2"
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB"  0.506962   0.500821   0.515332   0.529347   0.453294      0.506962   0.500821   0.515332   0.529347   0.453294   0.512858   0.482078   0.440029
```

#### 1.2.2 Run your own AUCell analysis

You need to prepare two input files before the analysis: pathway features gene file and expression profile file. 

##### pathway features gene file  

| Funtion                                                 | Format                 | Describe                                                     |
| ------------------------------------------------------- | ---------------------- | ------------------------------------------------------------ |
| read_gmt(fn::AbstractString)                            |                        | Read in a GMT file (MSigDB gene set format), where `fn` is the file path. |
| read_gsf(fn::AbstractString; delim::AbstractChar = ' ') | `.csv`, `.txt`, `.tsv` | Read in a general gene set file, where `fn` is the file path and the fields are separated by the `delim` character (default: white space). Each row represents a gene set and the first column is the name of the set and the rest are the genes in the set. |

##### 1.2.2.1 pathway features gene file

1) `read_gmt`: Read in a GMT file (MSigDB gene set format), where `fn` is the file path. `.gmt`
1) `read_gsf`: Read in a general gene set file, where `fn` is the file path and the fields are separated by the `delim` character (default: white space). Each row represents a gene set and the first column is the name of the set and the rest are the genes in the set. `.csv`, `.txt` and `.tsv` are supported. 

##### 1.2.2.2 expression profile file  

1) `read_mtx`: Read in the common 10X single-cell RNA expression file in the MTX format (unzipped).
2) `read_expr_matrix`: Read in an expression matrix stored in `fn` where its row names are stored in `rn` and column names are stored in `cn`.

```julia
julia> using AUCell
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia> pathway_AUC_main("matrix.mtx",
                        "features.tsv",
                        "barcodes.tsv",
                     	"fn_feature.gmt")
```

Other parameters can be set by passing the value to the corresponding keyword. 

```julia
pathway_AUC_main("matrix.mtx",
                 "features.tsv",
                 "barcodes.tsv",
              	 "fn_feature.gmt",
                 "fn_meta.txt";
           fn_meta_delim = '\t',
           fn_meta_group = "group",
        file_format_expr = "read_mtx",
             		   T = Int32,
             feature_col = 2,
             barcode_col = 1,
       feature_threshold = 30,
          cell_threshold = 200,
     file_format_feature = "read_gmt",
        fn_feature_delim = ' ',
    use_HALLMARK_pathway = "no",
                	mode = "AUCell",
           ncell_pseudo: = 0,
         auc_x_threshold = 1.0,
            remove_zeros = true,
         	use_testdata = "no",
             	work_dir = "./")
```

#### 1.2.3 Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `ncell_pseudo > 0`. 

```julia
julia> pathway_AUC_main("matrix.mtx",
                        "features.tsv",
                        "barcodes.tsv",
                     	"fn_feature.gmt",
						ncell_pseudo = 10)
```

 `ncell_pseudo` is the number of pseudobulk combined cells in each group. By default, profiling does not use the pseudo-bulk method (`ncell_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate how many cells are merged into a sample.


### 1.3 Optional Parameters

Below lists the optional keyword parameters and their default values.

| Parameter            | Parameter types | Default value    | Parameters to describe                                       |
| -------------------- | --------------- | ---------------- | ------------------------------------------------------------ |
| fn_expr              | AbstractString  | "matrix.mtx"     | MTX file path. (required).                                   |
| rn_expr              | AbstractString  | "features.tsv"   | features file path. (required)                               |
| cn_expr              | AbstractString  | "barcodes.tsv"   | barcodes file path. (required)                               |
| fn_feature           | AbstractString  | "fn_feature.gmt" |                                                              |
| fn_meta              | AbstractString  | "fn_meta.txt"    | Grouping information file path. Read in a meta data file with the first row assumed to be the header and the row names assumed to be the profile names (cell barcodes). |
| fn_meta_delim        | AbstractChar    | '\t'             | Delimiter of the metadata file data.                         |
| fn_meta_group        | AbstractString  | "group"          | Grouping information is specified by the column with the header name of `group`. If `group` is not found, the second column will be used. |
| file_format_expr     | AbstractString  | "read_mtx"       | There are two input modes "read_mtx" and "read_expr_matrix" for the expression profile file format. |
| T                    | Type            | Int32            | Express the storage format of the spectrum input variable.   |
| feature_col          | Int             | 2                | feature in the column.                                       |
| barcode_col          | Int             | 1                | barcode in the column.                                       |
| feature_threshold    | Int             | 30               | Include features (genes) detected in at least this many cells. |
| cell_threshold       | Int             | 200              | Include profiles (cells) where at least this many features are detected. |
| file_format_feature  | AbstractString  | "read_gmt"       | There are two input modes "read_gmt" and "read_gsf" for the file format of the features contained in the pathways. |
| fn_feature_delim     | AbstractChar    | ' '              | Delimiter of the pathway features file data.                 |
| use_HALLMARK_pathway | AbstractString  | "no"             | Whether to use the built-in HALLMARK pathways.               |
| mode                 | AbstractString  | "AUCell"         | "AUCell" is an optional mode to calculate AUC based on characteristic genes for two groups. "pathway_recluster" is subgroups based on pathway activation. |
| ncell_pseudo         | Int             | 0                | ncell_pseudo is the number of pseudobulk combined cells in each group. By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate how many cells are merged into a sample. |
| auc_x_threshold      | Float64         | 1.0              | Threshold for the X-axis (1-specificity) in the auc calculation, 0~auc_x_threshold. |
| remove_zeros         | Bool            | true             | Whether to remove all cells with zero gene expression values. |
| work_dir             | AbstractString  | "./"             | Working Directory.                                           |
| use_testdata         | AbstractString  | "no"             | Whether to use the default provided test data for analysis, yes or no. |







## 先写到这























### 1.4 Example output file

#### 1.4.1 result

- The file content is the pathways AUC value for each group sample. Behavioral pathways, listed as samples  (See )



#### 1.4.2 log file

- [aucell_result.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata_output/aucell_result.tsv)

## 2 Used in the R language

### 2.1 Installation

##### 2.1.1 You can install just like any other R packages by `JuliaCall`

```R
install.packages("JuliaCall")
```

##### 2.1.2 To use you must have a working installation of Julia. This can be easily done via: `JuliaCall`

```R
library(JuliaCall)
install_julia()
```

##### 2.1.3 which will automatically install and setup a version of Julia specifically for use with `JuliaCall`. Or you can do

```R
library(JuliaCall)
julia <-julia_setup()
```

##### 2.1.4 Download RankCompV3

```julia
julia_install_package_if_needed("RankCompV3")
```

### 2.2 Examples

#### 2.2.1 Quick Start

Run a test job with the input files distributed with the package.

```R
julia_library("RankCompV3")
result <- julia_do.call("reoa",list(use_testdata="yes"),need_return="Julia",show_value=FALSE)
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `reoa` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```R
> result
Julia Object of type Tuple{DataFrames.DataFrame, DataFrames.DataFrame, DataFrames.DataFrame}.
(19999×16 DataFrame
   Row │ Name     pval        padj        n11      n12      n13      n21       ⋯
       │ String   Float64     Float64     Float64  Float64  Float64  Float64   ⋯
───────┼────────────────────────────────────────────────────────────────────────
     1 │ DE1      0.237557    0.719072     1531.0     75.0      0.0   1213.0   ⋯
     2 │ DE2      0.280253    0.762245     2276.0    285.0      0.0    978.0
     3 │ DE3      0.0578256   0.358368     1574.0     37.0      0.0   1375.0
     4 │ DE4      0.155619    0.595079     1539.0    139.0      0.0    760.0
     5 │ DE5      0.0432114   0.30634      1166.0     12.0      0.0   1582.0   ⋯
     6 │ DE6      0.0572817   0.356322     1583.0     70.0      0.0    909.0
     7 │ DE7      0.00925067  0.123914     2182.0     64.0      0.0   3608.0
     8 │ DE8      0.036411    0.278039     1927.0     69.0      0.0   2155.0
   ⋮   │    ⋮         ⋮           ⋮          ⋮        ⋮        ⋮        ⋮      ⋱
 19993 │ EE19994  0.0833578   0.438012     1719.0    353.0      0.0   1143.0   ⋯
 19994 │ EE19995  0.0457029   0.316136     1840.0    215.0      0.0   1323.0
 19995 │ EE19996  0.843871    0.999171     1103.0    276.0      0.0    440.0
 19996 │ EE19997  6.70797e-5  0.00465755   1248.0     25.0      0.0   2058.0
 19997 │ EE19998  0.909702    1.0          1112.0    136.0      0.0    931.0   ⋯
 19998 │ EE19999  0.348734    0.827127     1475.0    304.0      0.0   1752.0
 19999 │ EE20000  0.699539    0.981932     1576.0    313.0      0.0    970.0
                                                9 columns and 19984 rows omitted, 19999×6 DataFrame
   Row │ Name     Sample1  Sample2  Sample3  Sample4  Sample5
       │ String   Int64    Int64    Int64    Int64    Int64
───────┼──────────────────────────────────────────────────────
     1 │ DE1           10       53       18       84       18
     2 │ DE2            5       22       59       11       33
     3 │ DE3           62       39       18       19        8
     4 │ DE4            6      116      131        3       49
     5 │ DE5           28       31       27       58       16
     6 │ DE6           91       81       26        3       29
     7 │ DE7            8        4       18        4        5
     8 │ DE8            3       31        5       32       11
   ⋮   │    ⋮        ⋮        ⋮        ⋮        ⋮        ⋮
 19993 │ EE19994        6      176       36        0        1
 19994 │ EE19995        3      103        8        0       48
 19995 │ EE19996       60       27      405      132        0
 19996 │ EE19997       23        4        0        0      116
 19997 │ EE19998      172       28       31        2      128
 19998 │ EE19999      102       49        3        0        6
 19999 │ EE20000       96        2        2       27      112
                                            19984 rows omitted, 19999×6 DataFrame
   Row │ Name     Sample76  Sample77  Sample78  Sample79  Sample80
       │ String   Int64     Int64     Int64     Int64     Int64
───────┼───────────────────────────────────────────────────────────
     1 │ DE1           200       180       176        26        24
     2 │ DE2            26        37        36        40        79
     3 │ DE3           157        81        89       220        97
     4 │ DE4            15       205        63        75       228
     5 │ DE5           505        92       119       157       261
     6 │ DE6           191       131        99       147        68
     7 │ DE7           116        38        34        58        64
     8 │ DE8            53        39       183        22        73
   ⋮   │    ⋮        ⋮         ⋮         ⋮         ⋮         ⋮
 19993 │ EE19994        88        15       192         9        35
 19994 │ EE19995       174        73        14        16        97
 19995 │ EE19996       138        61       407       182         0
 19996 │ EE19997       325        94       231       281        41
 19997 │ EE19998        24         0       116        43       656
 19998 │ EE19999         0       144         4       393         0
 19999 │ EE20000       294         0        26        98        18
                                                 19984 rows omitted)
```

#### 2.2.2 Run your own DEG analysis

You need to prepare two input files before the analysis: metadata file and expression matrix. You need to prepare two input files before the analysis: metadata file and expression matrix. `.rds`, `.csv`, `.txt`, `.tsv` and `.RData` files are supported.    

- **metadata file (required).**

 The metadata file contains at least two columns. The first column is the sample names, and the second column is the grouping information. Only two groups are supported at present, therefore, do not include more than two groups. 

 Column names for a metadata should be `Name` and `Group`. 

 See an example metadata file, [fn_metadata.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_metadata.txt).


- **expression matrix file (required).**

 The first column is the gene name and the column header should be `Name` and the rest columns are profiles for each cell or each sample. Each column header should be the sample name which appears in the metadata file.

 See an example expression matrix file, [fn_expr.txt](https://github.com/yanjer/RankCompV3.jl/blob/master/test/fn_expr.txt).

 Raw counts or expression values are recommended to use. Other values, e.g, FPKM, RPKM, TPM, log(counts) and log(normalized counts), can also be used, though normalization and batch effect removal are neither necessary nor recommended. 

Once the files are ready, you can carry out the DEG analysis with the default settings as follows. 

```R
julia_library("RankCompV3")
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia_do.call("reoa",list("/public/yanj/data/fn_expr.txt",
						"/public/yanj/data/fn_meta.txt"),need_return="Julia",show_value=FALSE)
```

Other parameters can be set by passing the value to the corresponding keyword. 

```R
julia_do.call("reoa",list("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_meta.txt",
    	expr_threshold = 0,
    	min_profiles = 0,
    	min_features = 0,
    	pval_reo = 0.01,
     	pval_deg = 0.05,
     	padj_deg = 0.05,
    	n_pseudo = 0,
    	use_hk_genes = "yes"
    	hk_file = "HK_genes_info.tsv",
    	gene_name_type = "ENSEMBL",
    	ref_gene_max = 3000,
    	ref_gene_min = 100
    	n_iter = 128,
    	n_conv = 5,
    	work_dir = "./",
    	use_testdata = "no"),need_return="Julia",show_value=FALSE)
```

#### 2.2.3 Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `n_pseudo = 1`. 

```R
julia_do.call("reoa",list("scRNA_expr.txt",
						"scRNA_metadata.txt",
        				n_pseudo = 1),need_return="Julia",show_value=FALSE)
```

By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined.


### 2.3 Optional Parameters

See 1.3 Optional Parameters.

### 2.4 Example output file

See 1.4 Example output file.

## 3 Specific experimental designs

In this chapter, we will present some examples of the use of experimental designs.

### 3.1 Datasets from different sources

Differential expression analysis often involves datasets from different sources, such as different sequencing platforms, different sample sources, or different sample processing methods. RankCompV3 supports direct column concatenation of expression profiles for integrated analysis of different datasets.

```julia
julia> expr
5×11 DataFrame
 Row │ gene_name  s1     s2     s3     c1     c8     c6     a1     a2     b1     b6
     │ String7    Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64
─────┼─────────────────────────────────────────────────────────────────────────────────
   1 │ DE1           10     53     18     84     18    200    180    176     26     24
   2 │ DE2            5     22     59     11     33     26     37     36     40     79
   3 │ DE3           62     39     18     19      8    157     81     89    220     97
   4 │ DE4            6    116    131      3     49     15    205     63     75    228
   5 │ DE5           28     31     27     58     16    505     92    119    157    261
```

Among them, s, c, a and b samples are data from different sources. Samples s and c are in the same group, and samples a and b are in the same group.

```julia
julia> meta
10×2 DataFrame
 Row │ sample_name  group
     │ String15     String7
─────┼──────────────────────
   1 │ s1           group1
   2 │ s2           group1
   3 │ s3           group1
   4 │ c1           group1
   5 │ c8           group1
   6 │ c6           group1
   7 │ a1           group2
   8 │ a2           group2
   9 │ b1           group2
  10 │ b6           group2
```

### 3.2 Two or more groups

When multiple cell types are included, the characteristics of each cell class need to be analyzed and each cell class compared to the rest. Suppose there are three types of cells, such as cell types s, c, a, and b.

```julia
julia> meta
10×2 DataFrame
 Row │ sample_name  group
     │ String       String
─────┼────────────────────────
   1 │ s1           celltype1
   2 │ s2           celltype1
   3 │ s3           celltype1
   4 │ c1           celltype2
   5 │ c8           celltype2
   6 │ c6           celltype2
   7 │ a1           celltype3
   8 │ a2           celltype3
   9 │ b1           celltype4
  10 │ b6           celltype4
```

By default, conditions are listed in order of precedence.

```julia
julia> unique(meta.group)
4-element Vector{String}:
 "celltype1"
 "celltype2"
 "celltype3"
 "celltype4"
```

### 3.3 Combined sample

For analyzing single-cell data, drop-out phenomenon often exists. The `n_pseudo` parameter of pseudobulk mode of RankCompV3 can be used to combine different cells to reduce the effect of drop-out.

```julia
julia> reoa("/public/yanj/data/fn_expr.txt",
    	"/public/yanj/data/fn_meta.txt",
    	n_pseudo = 10)
```

By default, profiling does not use the pseudo-bulk method (`n_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate the number of samples in each group after pseudo-bulk is combined.























