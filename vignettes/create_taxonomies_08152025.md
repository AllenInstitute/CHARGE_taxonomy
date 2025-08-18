# Create taxonomies

This script generates CHARGE RData files for [all remaining taxonomies available in AIT format as of 15 August 2025](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/taxonomies.md). This script can be updated with additional taxonomies as needed.


**NOTE:** This script was run in the Docker environment `"docker://alleninst/scrattch:1.1.2"`. It will likely fail if run elsewhere.

***

## Prerequisites

Prior to running this script, download `"chargeTaxonomy.R"`, `"utils.R"`, and `"constellation.R"` from [CHARGE_taxonomy](https://github.com/AllenInstitute/CHARGE_taxonomy) into your working directory. Then enter the above Docker environment and navigate to the same working directory.

Alternatively, install this R library with:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)
remotes::install_github("AllenInstitute/CHARGE_taxonomy")
```

## Set up the environment

```R
## Load Relevant Libraries
library(scrattch.taxonomy)  # Note, this loads all other relevant libraries.

# Option 1: source CHARGE_taxonomy functions
source("chargeTaxonomy.R")
source("utils.R")
source("constellation.R")
# Option 2: load library
library(CHARGE.taxonomy) # Note _ to . change

## Set-up global variables
# Where are all the files located?
S3.folder <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/"
# Subsample to keep a lot of cells for the stats, but not necessarily all of them 
subsample <- 1000  
```

-----

## Run on the updated MTG 2018 taxonomy

This update addresses an error in the UMAP coordinates from the original release.  We also are shortening the subclass names in this taxonomy for convenience.

```r
## Read in taxonomy
file.name   <- "Human_MTG_SMART_seq_08082025.h5ad"
AIT.anndata <- loadTaxonomy(paste0(S3.folder,file.name))

## Shorten subclasses
levs = levels(AIT.anndata$obs$subclass)
subs = as.character(AIT.anndata$obs$subclass)
rems = c("MTG ", "human ", "cerebral cortex ", "GABAergic ", "-expressing", "Human ", "FGFR3 ")
for (rem in rems){
  subs = gsub(rem,"",subs)
  levs = gsub(rem,"",levs)
}
AIT.anndata$obs$subclass <- factor(subs,levels=levs)
AIT.anndata$uns$cluster_info <- NULL

## Run chargeTaxonomy script
charge.file <- gsub(".h5ad","_CHARGE.RData",file.name)
chargeTaxonomy(AIT.anndata = AIT.anndata, subsample = subsample, charge.file.name = charge.file)
```

-----

## Repeat for all available taxonomies

This section processes all remaining taxonomies as of 15 August 2025.

```r
file.names = c(
  "Mouse_VISp_ALM_SMART_seq_04042025.h5ad",
  "Human_LGN_SMART_seq_04042025.h5ad",
  "Mouse_LGN_SMART_seq_04042025.h5ad",
  "Macaque_LGN_SMART_seq_04042025.h5ad",
  "Human_neocortex_SMART_seq_04042025.h5ad",
  "Human_M1_10X_seq_04042025.h5ad",
  "Mouse_cortex_hippocampus_10X_seq_04042025.h5ad",
  "Mouse_cortex_hippocampus_SMART_seq_04042025.h5ad",
  "Human_MTG_SEAAD_04042025.h5ad"
)

for (file.name in file.names[2:9]){
  print(paste0("############################ RUNNING ",file.names," #################"))
  AIT.anndata <- loadTaxonomy(paste0(S3.folder,file.name))
  charge.file <- gsub(".h5ad","_CHARGE.RData",file.name)
  chargeTaxonomy(AIT.anndata = AIT.anndata, subsample = subsample, charge.file.name = charge.file)
}
```
