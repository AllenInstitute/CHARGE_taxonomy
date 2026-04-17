# CHARGE_taxonomy
Cell HierARchy Gene Explorer (CHARGE) uses fast, cluster-centric approaches to find global or local marker genes, differentially expressed genes, and genes following user selected gradients. For more details on CHARGE, please visit [the CHARGE GitHub page](https://github.com/alleninstitute/CHARGE).

<img width="150" alt="CHARGE logo" src="https://github.com/AllenInstitute/CHARGE/blob/main/R/www/CHARGE_logo.png">

**`CHARGE_taxonomy` is the GitHub repository and associated R library for generating RData files for CHARGE.**  Starting either from any taxonomy in [Allen Institute Taxonomy (AIT) format](https://github.com/AllenInstitute/AllenInstituteTaxonomy/) or from clustered and annotated cell x gene expression data, `CHARGE_taxonomy` packages the relevant statistics and other variables into a .RData file required for CHARGE. 

As of April 2026, a second hybrid implementation utilizing both R and python (`chargeTaxonomyHybrid`) extends functionality to larger taxonomies (e.g., > ~200,000 cells)

## Installation

CHARGE_taxonomy should be run in the Docker environment `"docker://alleninst/scrattch:1.1.4.1"`. It may fail if run elsewhere. From the docker environment, install this R library with:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)
remotes::install_github("AllenInstitute/CHARGE_taxonomy")
```

## Example

CHARGE_taxonomy can be run largely with a single function. 

For example if you want to use the R-only version:

```R
AIT.file    <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SMART_seq_08082025.h5ad"
AIT.anndata <- loadTaxonomy(AIT.file)
chargeTaxonomy(AIT.anndata)
```

See the [manual file](https://github.com/AllenInstitute/CHARGE_taxonomy/blob/main/man/chargeTaxonomy.Rd) or type `?chargeTaxonomy` in R for more details. 

If you want to use the R + python hybrid version:

```R
AIT.url  <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SMART_seq_08082025.h5ad"
AIT.file <- "Human_MTG_SMART_seq_08082025.h5ad"
options(timeout = 9000)
download.file(AIT.url, destfile = AIT.file, mode = "wb", timeout = 9000)
chargeTaxonomyHybrid(AIT.file)
```

See the [manual file](https://github.com/AllenInstitute/CHARGE_taxonomy/blob/main/man/chargeTaxonomyHybrid.Rd) or type `?chargeTaxonomyHybrid` in R for more details. 

[This R script](https://github.com/AllenInstitute/CHARGE_taxonomy/blob/main/vignettes/create_taxonomies_08152025.R) generates CHARGE RData files for [all taxonomies available in AIT format as of 15 August 2025](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/taxonomies.md).  Here is a [version better formatted for web viewing](https://github.com/AllenInstitute/CHARGE_taxonomy/blob/main/vignettes/create_taxonomies_08152025.md).

## Reporting issues or suggestions

Please share any comments, suggestions, bugs, or any other thoughts using the button in the app, or by [submitting an issue](https://github.com/AllenInstitute/CHARGE_taxonomy/issues).

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/CHARGE_taxonomy/blob/master/LICENSE

## Level of Support

We only plan updates to this tool as needed to ensure compatability with the main [CHARGE repo](https://github.com/alleninstitute/CHARGE).
