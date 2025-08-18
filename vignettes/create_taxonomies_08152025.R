#######################################################
## NOTE:  This was run in the docker environment     ##
##        "docker://alleninst/scrattch:1.1.2".       ##
##        It will likely fail if run elsewhere.      ##
#######################################################


#######################################################
## Prior to running this script, download "chargeTaxonomy.R", "utils.R", and "constellation.R" from https://github.com/AllenInstitute/CHARGE_taxonomy into your working directory or install this R library with `remotes::install_github("AllenInstitute/CHARGE_taxonomy")`.  Then enter the above docker environment and navigate to the same working directory.

## Load relevant libraries

library(scrattch.taxonomy)  # Note, this loads all other relevant libraries.

# Option 1: source CHARGE_taxonomy functions
source("chargeTaxonomy.R")
source("utils.R")
source("constellation.R")
# Option 2: load library
library(CHARGE.taxonomy)

## Set-up global variables

S3.folder <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/"
subsample <- 1000  # Keep a lot of cells for the stats, but not necessarily all of them 


#######################################################
## Run on the updated MTG 2018 taxonomy with some minor modifications

file.name   <- "Human_MTG_SMART_seq_08082025.h5ad"
AIT.anndata <- loadTaxonomy(paste0(S3.folder,file.name))

## For convenience, we are shortening the subclass names in this taxonomy
levs = levels(AIT.anndata$obs$subclass)
subs = as.character(AIT.anndata$obs$subclass)
rems = c("MTG ", "human ", "cerebral cortex ", "GABAergic ", "-expressing", "Human ", "FGFR3 ")
for (rem in rems){
  subs = gsub(rem,"",subs)
  levs = gsub(rem,"",levs)
}
AIT.anndata$obs$subclass <- factor(subs,levels=levs)
AIT.anndata$uns$cluster_info <- NULL

charge.file <- gsub(".h5ad","_CHARGE.RData",file.name)
chargeTaxonomy(AIT.anndata = AIT.anndata, subsample = subsample, charge.file.name = charge.file)


#######################################################
## Finally, repeat for all available taxonomies as of 15 August 2025 on
#    https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/taxonomies.md
#    This can be updated with additional taxonomies as needed.

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
