## Scripts, results and other files

### Data files
cts_all.csv - raw counts for all samples

samples_all.csv - bat metadata for all samples 

### R session info

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] flashClust_1.01-2           WGCNA_1.73                  fastcluster_1.2.6           dynamicTreeCut_1.63-1       lubridate_1.9.4            
 [6] forcats_1.0.0               stringr_1.5.1               purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                
[11] tibble_3.2.1                tidyverse_2.0.0             ggupset_0.4.1               pheatmap_1.0.12             ggplotify_0.1.2            
[16] ggbeeswarm_0.7.2            ggridges_0.5.6              ggrepel_0.9.6               org.Hs.eg.db_3.19.1         AnnotationDbi_1.66.0       
[21] clusterProfiler_4.12.6      sva_3.52.0                  genefilter_1.86.0           mgcv_1.9-1                  nlme_3.1-167               
[26] variancePartition_1.34.0    BiocParallel_1.38.0         gridExtra_2.3               patchwork_1.3.0             reshape2_1.4.4             
[31] ggplot2_4.0.1               edgeR_4.2.2                 limma_3.60.6                dplyr_1.1.4                 DESeq2_1.44.0              
[36] SummarizedExperiment_1.34.0 Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.5.0           GenomicRanges_1.56.2       
[41] GenomeInfoDb_1.40.1         IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.50.0        

loaded via a namespace (and not attached):
  [1] fs_1.6.5                bitops_1.0-9            enrichplot_1.24.4       httr_1.4.7              RColorBrewer_1.1-3     
  [6] doParallel_1.0.17       numDeriv_2016.8-1.1     tools_4.4.2             backports_1.5.0         R6_2.5.1               
 [11] lazyeval_0.2.2          withr_3.0.2             preprocessCore_1.66.0   cli_3.6.3               scatterpie_0.2.4       
 [16] mvtnorm_1.3-3           S7_0.2.1                yulab.utils_0.2.0       gson_0.1.0              foreign_0.8-88         
 [21] DOSE_3.30.5             R.utils_2.12.3          rstudioapi_0.17.1       impute_1.78.0           RSQLite_2.3.9          
 [26] generics_0.1.3          gridGraphics_0.5-1      gtools_3.9.5            GO.db_3.19.1            Matrix_1.7-2           
 [31] abind_1.4-8             R.methodsS3_1.8.2       lifecycle_1.0.4         gplots_3.2.0            qvalue_2.36.0          
 [36] SparseArray_1.4.8       grid_4.4.2              blob_1.2.4              crayon_1.5.3            lattice_0.22-6         
 [41] cowplot_1.1.3           annotate_1.82.0         KEGGREST_1.44.1         pillar_1.10.1           knitr_1.49             
 [46] fgsea_1.30.0            boot_1.3-31             corpcor_1.6.10          codetools_0.2-20        fastmatch_1.1-6        
 [51] glue_1.8.0              ggfun_0.1.8             data.table_1.16.4       vctrs_0.6.5             png_0.1-8              
 [56] treeio_1.28.0           Rdpack_2.6.2            gtable_0.3.6            cachem_1.1.0            xfun_0.50              
 [61] rbibutils_2.3           S4Arrays_1.4.1          tidygraph_1.3.1         reformulas_0.4.0        survival_3.8-3         
 [66] iterators_1.0.14        statmod_1.5.0           pbkrtest_0.5.3          ggtree_3.12.0           bit64_4.6.0-1          
 [71] EnvStats_3.0.0          vipor_0.4.7             KernSmooth_2.23-26      rpart_4.1.24            colorspace_2.1-1       
 [76] DBI_1.2.3               Hmisc_5.2-2             nnet_7.3-20             tidyselect_1.2.1        bit_4.5.0.1            
 [81] compiler_4.4.2          httr2_1.1.0             htmlTable_2.4.3         DelayedArray_0.30.1     shadowtext_0.1.4       
 [86] checkmate_2.3.2         scales_1.4.0            caTools_1.18.3          remaCor_0.0.18          rappdirs_0.3.3         
 [91] digest_0.6.37           minqa_1.2.8             rmarkdown_2.29          aod_1.3.3               XVector_0.44.0         
 [96] RhpcBLASctl_0.23-42     htmltools_0.5.8.1       pkgconfig_2.0.3         base64enc_0.1-3         lme4_1.1-36            
[101] fastmap_1.2.0           rlang_1.1.5             htmlwidgets_1.6.4       UCSC.utils_1.0.0        farver_2.1.2           
[106] jsonlite_1.8.9          GOSemSim_2.30.2         R.oo_1.27.0             RCurl_1.98-1.16         magrittr_2.0.3         
[111] Formula_1.2-5           GenomeInfoDbData_1.2.12 Rcpp_1.1.0              ape_5.8-1               viridis_0.6.5          
[116] stringi_1.8.4           ggraph_2.2.1            zlibbioc_1.50.0         MASS_7.3-64             plyr_1.8.9             
[121] parallel_4.4.2          Biostrings_2.72.1       graphlayouts_1.2.2      splines_4.4.2           hms_1.1.3              
[126] locfit_1.5-9.11         igraph_2.1.4            XML_3.99-0.18           evaluate_1.0.3          nloptr_2.1.1           
[131] tzdb_0.4.0              foreach_1.5.2           tweenr_2.0.3            polyclip_1.10-7         ggforce_0.4.2          
[136] broom_1.0.7             xtable_1.8-4            fANCOVA_0.6-1           tidytree_0.4.6          viridisLite_0.4.2      
[141] lmerTest_3.1-3          aplot_0.2.4             memoise_2.0.1           beeswarm_0.4.0          cluster_2.1.8          
[146] timechange_0.3.0  
