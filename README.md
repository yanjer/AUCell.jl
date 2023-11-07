# AUCell.jl

AUCell.jl is an algorithm for cell reclassification based on AUC values by feature genes in the pathway.

## 1 Used in the Julia language

### 1.1 Installation

The algorithm is implemented in Julia. Version 1.7 or later is recommended. The simpliest way to install is using the `Pkg` facility in Julia. 

```julia
using Pkg
Pkg.add("AUCell.jl")
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

1) `read_gmt`: Read in a GMT file (MSigDB gene set format), where `fn` is the file path. `.gmt` (See [fn_feature.gmt](https://github.com/yanjer/AUCell/blob/master/HALLMARK_pathway/h_all_v2023_1_Hs_symbols.gmt))
1) `read_gsf`: Read in a general gene set file, where `fn` is the file path and the fields are separated by the `delim` character (default: white space). Each row represents a gene set and the first column is the name of the set and the rest are the genes in the set. `.csv`, `.txt` and `.tsv` are supported. (See `.csv`: [fn_feature.csv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.csv) or `.txt`:  [fn_feature.txt](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.txt) or `.tsv`:  [fn_feature.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.tsv))

##### 1.2.2.2 expression profile file  

1. `read_mtx`: Read in the common 10X single-cell RNA expression file in the MTX format (unzipped). (See `fn`: [matrix.mtx](https://github.com/yanjer/AUCell/blob/master/test/matrix.mtx), `rn`: [features.tsv](https://github.com/yanjer/AUCell/blob/master/test/features.tsv), `cn`: [barcodes.tsv](https://github.com/yanjer/AUCell/blob/master/test/barcodes.tsv))

2. `read_expr_matrix`: Read in an expression matrix stored in `fn` where its row names are stored in `rn` and column names are stored in `cn`.   (See `fn`: [matrix.csv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.csv) (`.csv`) or [matrix.txt](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.txt) (`.txt`) or [matrix.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.tsv) (`.tsv`); `rn`: [features.tsv](https://github.com/yanjer/AUCell/blob/master/test/features.tsv), `cn`: [barcodes.tsv](https://github.com/yanjer/AUCell/blob/master/test/barcodes.tsv))

##### 1.2.2.3 pathway features gene file

   `read_meta`: Read in a meta data file with the first row assumed to be the header and the row names assumed to be the profile names (cell barcodes). Grouping information is specified by the column with the header name of `group`. If `group` is not found, the second column will be used. It returns the grouped profile names (vector of vectors) and group names. (See [fn_meta.txt](https://github.com/yanjer/AUCell/blob/master/test/fn_meta.txt))

```julia
julia> using AUCell
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia> pathway_AUC_main("matrix.mtx",
                        "features.tsv",
                        "barcodes.tsv",
                     	"fn_feature.gmt",
                 		"fn_meta.txt")
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
         	   rem_delim = ' ',
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
                 		"fn_meta.txt";
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
| fn_feature           | AbstractString  | "fn_feature.gmt" | Pathway feature gene set file path. (required)               |
| fn_meta              | AbstractString  | "fn_meta.txt"    | Grouping information file path. Read in a meta data file with the first row assumed to be the header and the row names assumed to be the profile names (cell barcodes). |
| fn_meta_delim        | AbstractChar    | '\t'             | Delimiter of the metadata file data.                         |
| fn_meta_group        | AbstractString  | "group"          | Grouping information is specified by the column with the header name of `group`. If `group` is not found, the second column will be used. |
| file_format_expr     | AbstractString  | "read_mtx"       | There are two input modes "read_mtx" and "read_expr_matrix" for the expression profile file format. |
| T                    | Type            | Int32            | Express the storage format of the spectrum input variable.   |
| feature_col          | Int             | 2                | feature in the column.                                       |
| barcode_col          | Int             | 1                | barcode in the column.                                       |
| rem_delim            | AbstractChar    | ' '              | Enter the file separator when file_format_expr is "read_expr_matrix". |
| feature_threshold    | Int             | 30               | Include features (genes) detected in at least this many cells. |
| cell_threshold       | Int             | 200              | Include profiles (cells) where at least this many features are detected. |
| file_format_feature  | AbstractString  | "read_gmt"       | There are two input modes "read_gmt" and "read_gsf" for the file format of the features contained in the pathways. |
| fn_feature_delim     | AbstractChar    | ' '              | Delimiter of the pathway features file data.                 |
| use_HALLMARK_pathway | AbstractString  | "no"             | Whether to use the built-in HALLMARK pathways.               |
| mode                 | AbstractString  | "AUCell"         | "AUCell" is an optional mode to calculate AUC based on characteristic genes for two groups. "pathway_recluster" is subgroups based on pathway activation. |
| ncell_pseudo         | Int             | 0                | ncell_pseudo is the number of pseudobulk combined cells in each group. By default, profiling does not use the pseudo-bulk method (`ncell_pseudo= 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate how many cells are merged into a sample. |
| auc_x_threshold      | Float64         | 1.0              | Threshold for the X-axis (1-specificity) in the auc calculation, 0~auc_x_threshold. |
| remove_zeros         | Bool            | true             | Whether to remove all cells with zero gene expression values. |
| work_dir             | AbstractString  | "./"             | Working Directory.                                           |
| use_testdata         | AbstractString  | "no"             | Whether to use the default provided test data for analysis, yes or no. |

### 1.4 Example output file

#### 1.4.1 result

- The file content is the pathways AUC value for each group sample. Behavioral pathways, listed as samples.  (See [aucell_result.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata_output/aucell_result.tsv)

#### 1.4.2 log file

- [AUCell_testdata.log](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata_output/AUCell_testdata.log)

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

##### 2.1.4 Download AUCell

```R
julia_install_package_if_needed("AUCell")
```

### 2.2 Examples

#### 2.2.1 Quick Start

Run a test job with the input files distributed with the package.

```R
julia_library("AUCell")
result <- julia_do.call("pathway_AUC_main",list(use_testdata="yes"),need_return="Julia",show_value=FALSE)
```

The analysis results and a few plots will be generated and saved in the current work directory. They are also returned by the `pathway_AUC_main` function and can be captured by assign the returned values to a variable,  e.g., `result` in the above example.  

The first return value is a DataFrame, where rows are genes and columns are statistical values for each gene. All the genes passing the basic preprocessing step are retained. 


```R
> result
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

#### 2.2.2 Run your own DEG analysis

You need to prepare four input files before the analysis: metadata file and expression matrix. You need to prepare two input files 

##### 2.2.2.1 pathway features gene file 

1) read_gmt`: Read in a GMT file (MSigDB gene set format), where `fn` is the file path. `.gmt` (See [fn_feature.gmt](https://github.com/yanjer/AUCell/blob/master/HALLMARK_pathway/h_all_v2023_1_Hs_symbols.gmt))
1) `read_gsf`: Read in a general gene set file, where `fn` is the file path and the fields are separated by the `delim` character (default: white space). Each row represents a gene set and the first column is the name of the set and the rest are the genes in the set. `.csv`, `.txt` and `.tsv` are supported. (See `.csv`: [fn_feature.csv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.csv) or `.txt`:  [fn_feature.txt](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.txt) or `.tsv`:  [fn_feature.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/fn_feature.tsv))

##### 2.2.2.2 expression profile file  

1. `read_mtx`: Read in the common 10X single-cell RNA expression file in the MTX format (unzipped). (See `fn`: [matrix.mtx](https://github.com/yanjer/AUCell/blob/master/test/matrix.mtx), `rn`: [features.tsv](https://github.com/yanjer/AUCell/blob/master/test/features.tsv), `cn`: [barcodes.tsv](https://github.com/yanjer/AUCell/blob/master/test/barcodes.tsv))

2. `read_expr_matrix`: Read in an expression matrix stored in `fn` where its row names are stored in `rn` and column names are stored in `cn`.   (See `fn`: [matrix.csv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.csv) (`.csv`) or [matrix.txt](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.txt) (`.txt`) or [matrix.tsv](https://github.com/yanjer/testdata-output/blob/master/AUCell_testdata/matrix.tsv) (`.tsv`); `rn`: [features.tsv](https://github.com/yanjer/AUCell/blob/master/test/features.tsv), `cn`: [barcodes.tsv](https://github.com/yanjer/AUCell/blob/master/test/barcodes.tsv))

##### 2.2.2.3 pathway features gene file

   `read_meta`: Read in a meta data file with the first row assumed to be the header and the row names assumed to be the profile names (cell barcodes). Grouping information is specified by the column with the header name of `group`. If `group` is not found, the second column will be used. It returns the grouped profile names (vector of vectors) and group names. (See [fn_meta.txt](https://github.com/yanjer/AUCell/blob/master/test/fn_meta.txt))

Once the files are ready, you can carry out the AUCell analysis with the default settings as follows. 

```R
julia_library("AUCell")
# Use the default values for the following other parameters. If you want to modify the parameters, add them directly.
julia_do.call("reoa",list("matrix.mtx",
                         "features.tsv",
                         "barcodes.tsv",
                         "fn_feature.gmt",
                         "fn_meta.txt"),need_return="Julia",show_value=FALSE)
```

Other parameters can be set by passing the value to the corresponding keyword. 

```R
julia_do.call("reoa",list("matrix.mtx",
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
         			   rem_delim = ' ',
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
                        work_dir = "./"),need_return="Julia",show_value=FALSE)
```

#### 2.2.3 Pseudobulk method

For scRNA-seq data, one can carry out a pseudobulk analysis. Rather than using the original single-cell profiles, pseudobulk profiles can be generated and used for DEG analysis. In this method, a random subset of cells from a group is aggregated into a pseudo-bulk profile. 

The pseudobulk method can be turned on by setting `ncell_pseudo > 0`. 

```R
julia_do.call("reoa",list("matrix.mtx",
                        "features.tsv",
                        "barcodes.tsv",
                     	"fn_feature.gmt",
                 		"fn_meta.txt";
						ncell_pseudo = 10),need_return="Julia",show_value=FALSE)
```

 `ncell_pseudo` is the number of pseudobulk combined cells in each group. By default, profiling does not use the pseudo-bulk method (`ncell_pseudo = 0`). 0 indicates that the pseudo-bulk mode is not used, and other values indicate how many cells are merged into a sample.


### 2.3 Optional Parameters

**See** 1.3 Optional Parameters.

### 2.4 Example output file

**See** 1.4 Example output file.
