# iDINGO: Integrative Differential Network Analysis in Genomics

iDINGO is a pathway-based method for estimating group-specific conditional dependencies and inferring differential networks between groups, based on genomic data. This can be done in a single-platform framework (for example, RNA-Seq data) or an integrative multi-platform framework (microRNA -> RNA -> Proteomics, where data from all three platforms are available for every sample). 

## Using iDINGO

We recommend filtering genomic data to fewer than 300 genes, generally filtered using a pathway/pathways of interest. Single-platform analyses are run using `dingo` with an nxp matrix, where n is the number of samples. Multi-platform analyses are run using `idingo`, with up to 3 separate data matrices containing the same n samples. For both `dingo` and `idingo`, the number of bootstraps is specified by `B` (we recommend at least 100). Parallel computing can speed this step up significantly, by setting the number of `cores`.  Finally, the `plotNetwork` function plots the differential network identified by `dingo` or `idingo`, based on a user-specified p-value or differential score threshold.
