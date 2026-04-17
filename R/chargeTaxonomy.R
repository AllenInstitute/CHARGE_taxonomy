# Helper function for reading stats matrix
# 
.read_stats_matrix <- function(path) {
  df <- data.table::fread(path, data.table = FALSE)
  rn <- df[[1]]
  mat <- as.matrix(df[,-1, drop = FALSE])
  rownames(mat) <- rn
  mode(mat) <- 'numeric'
  mat
}

# Helper function for reading counts
# 
.read_count_n <- function(path) {
  df <- data.table::fread(path, data.table = FALSE)
  out <- df$n_cells
  names(out) <- df$group
  as.numeric(out) -> out
  names(out) <- df$group
  out
}

# Helper function for reading all stats
# 
.read_stats_bundle <- function(stats.dir) {
  list(
    sums   = .read_stats_matrix(file.path(stats.dir, 'sums.csv')),
    counts = .read_stats_matrix(file.path(stats.dir, 'counts.csv')),
    props  = .read_stats_matrix(file.path(stats.dir, 'props.csv')),
    means  = .read_stats_matrix(file.path(stats.dir, 'means.csv')),
    sds    = .read_stats_matrix(file.path(stats.dir, 'sds.csv')),
    count_n = .read_count_n(file.path(stats.dir, 'count_n.csv'))
  )
}



#' Save an RData object for use with CHARGE app
#' 
#' Write an RData object with the following variables
#'  - counts_n - vector of total number of cells per cluster
#'  - counts - cell x cluster matrix where each value represents the number of cells in a cluster with at least one read for each gene 
#'  - props - cell x cluster matrix where each value represents the fraction of cells in a cluster with at least one read for each gene
#'  - sums - cell x cluster matrix where each value represents the total number of counts of that gene for that cluster 
#'  - means - cell x cluster matrix where each value represents the average counts of that gene for that cluster 
#'  - sds - cell x cluster matrix where each value represents the standard deviation of counts of that gene for that cluster
#'  - hierarchy - cell type hierarchy stored in uns$hierarchy, but converted to a corrected ordered character vector
#'  - cluster_info - a table of cluster information that encodes the hierarchy for building sunburst plots
#'  - constellation - a list of plotly objects holding the constellation diagrams for each levels of the hierarchy
#'  **NOTE: This function is tested on "docker://alleninst/scrattch:1.1.4.1". We strongly encourage using this docker environment.**
#'
#' @param AIT.anndata A reference taxonomy anndata object.  If provided, must contain counts in X or raw$X, uns$hierarchy, and obs (with columns for the hierarchy), and can optionally contain both X and raw$X and default embeddings and variable genes.  See https://github.com/AllenInstitute/AllenInstituteTaxonomy for details on expected formatting. If provided, this variable takes priority over the variables for ingesting these items separately.
#' @param cell_counts A sparse matrix of read counts where ROWS are cells and COLUMNS are genes, and both row names containing gene names (symbols or IDs) and column names containing unique IDs are included. Ignored if AIT.anndata is provided.
#' @param norm_counts OPTIONAL. If provided, a sparse matrix of log-normalized read counts where ROWS are cells and COLUMNS are genes, and both row names containing gene names (symbols or IDs) and column names containing unique IDs are included.  If not provided log2(CPM+1) is used. Ignored if AIT.anndata is provided.
#' @param metadata A matrix of cell metadata where ROWS are cells, COLUMNS are any metadata but that must include columns for each cell type level in the hierarchy. Row names in metadata MUST match row names in cell_counts. Ignored if AIT.anndata is provided.
#' @param hierarchy A character vector listing the columns corresponding to the cell hierarchy, with the highest resolution type (e.g., cluster) listed first and the lowest resolution type (e.g., class) listed last. Ignored if AIT.anndata is provided.
#' @param embedding OPTIONAL. A matrix where ROWS are cells and the first two COLUMNS correspond to x and y dimensions of a two dimensional embedding (e.g., UMAP or tSNE), with values corresponding to the specific X,Y coordinates for the embedding. Any other columns are ignored. Row names in metadata MUST match row names in cell_counts. Ignored if AIT.anndata is provided.
#' @param variable.genes OPTIONAL. A character vector of genes to use for calculating the embedding (if not provided) and the constellation diagram. Typically this corresponds to a set of variable or differentially expressed genes. Must be a subset of gene names included as column names for cell_counts. Ignored if AIT.anndata is provided.
#' @param check.taxonomy Should the function check if the input variable is a valid AIT taxonomy (default = TRUE)
#' @param subsample How much subsampling should be done before running statistics. Default (recommended!) is none, but subsampling can be done to speed up calculations.
#' @param stats.dir A directory that includes precomputed statistics from chargeTaxonomyHybrid. Typically this should be left as default (NULL) with chargeTaxonomyHybrid run directly if python is to be used.
#' @param weight.by Should statistics for higher levels of the hierarchy be weighted by "cell" (default), whereby each statistic is recalculated on all the cells in a given group for each level, or by "cluster", whereby statistics are averaged, counting each item from the first level of the hierarchy (e.g., cluster) evenly.
#' @param charge.file.name File name (and path) to write CHARGE file to
#' @param seed The seed to use for reproducibility
#' 
#' @import Matrix
#' 
#' @examples
#' \dontrun{
#'   AIT.file    <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SMART_seq_08082025.h5ad"
#'   AIT.anndata <- loadTaxonomy(AIT.file)
#'   chargeTaxonomy(AIT.anndata)
#' }
#'
#' @export
chargeTaxonomy <- function(AIT.anndata = NULL,
                           cell_counts = NULL,      # ignored if AIT.anndata is provided
                           norm_counts = NULL,      # OPTIONAL; ignored if AIT.anndata is provided
                           metadata    = NULL,      # ignored if AIT.anndata is provided
                           hierarchy   = NULL,      # ignored if AIT.anndata is provided
                           embedding   = NULL,      # OPTIONAL; ignored if AIT.anndata is provided
                           variable.genes   = NULL, # OPTIONAL; ignored if AIT.anndata is provided
                           subsample   = 100000000,
						   stats.dir = NULL,
                           weight.by   = "cell",
                           charge.file.name = "CHARGE.RData",
                           seed = 42)
{
  #############################################################
  print("===== Setting up variables =====")
  
  ##########################
  # metadata
  print("... read metadata and underlying hierarchy.")
  if(!is.null(AIT.anndata)){
    library(anndata)
    metadata <- AIT.anndata$obs
  } else {
    if(is.null(metadata)) error("AIT.anndata or metadata must be provided to run chargeTaxonomy.")
  }
  
  ##########################
  # hierarchy
  if(!is.null(AIT.anndata)){
    hierarchy = names(AIT.anndata$uns$hierarchy)[order(-as.numeric(AIT.anndata$uns$hierarchy))]
  } else {
    if(is.null(hierarchy)) error("AIT.anndata or hierarchy must be provided to run chargeTaxonomy.")
  }
  
  ##########################
  # cluster_info
  print("... define cluster information.")
  cluster_vector = metadata[,hierarchy[1]]
  all_clusters   = unique(cluster_vector)
  if(is.factor(all_clusters)) 
    all_clusters = levels(all_clusters)
  cluster_info = metadata[match(all_clusters,metadata[,hierarchy[1]]),hierarchy]
  
  # Rename anything called "_id" or "_label" to avoid breaking some scripts
  hierarchy_old <- hierarchy
  hierarchy <- gsub("_id","",hierarchy)
  hierarchy <- gsub("_label","",hierarchy)
  colnames(cluster_info) <- hierarchy
  names(hierarchy_old) <- hierarchy
  
  ##########################
  # cluster colors/ids
  print("... add cluster ids and colors ")
  ## NOTE: I'LL NEED TO UPDATE THIS PART ONCE I SORT OUT HOW TO EMBED CLUSTER 
  ##       COLORS (AND CLUSTER ORDER?) WITHIN THE CLUSTER_INFO DATA FRAME
  cluster_info$sample_name = paste0("i",1:dim(cluster_info)[1])
  for (i in 1:dim(cluster_info)[2]){
    cluster_info[,i] <- factor(cluster_info[,i],levels = unique(cluster_info[,i]))
  }    
  cluster_info <- auto_annotate(cluster_info)
  
  ##########################
  # Subsample vector
  print("... define subsampling, if any.")
  set.seed(seed)
  keep_sample <- subsampleCells(cluster_vector,subsample)
  
  #############################################################
  ## Heavy-data block: either read from Python outputs or fall back to original R logic
  if(!is.null(stats.dir)){
    print('===== Reading Python-precomputed statistics =====')
    stats.bundle <- .read_stats_bundle(stats.dir)
    sums   <- stats.bundle$sums
    counts <- stats.bundle$counts
    props  <- stats.bundle$props
    means  <- stats.bundle$means
    sds    <- stats.bundle$sds
    count_n <- stats.bundle$count_n
    
    all.genes <- rownames(sums)
    
    obs.info.path <- file.path(stats.dir, 'obs_info.csv')
    if(file.exists(obs.info.path)){
      obs.df <- data.table::fread(obs.info.path, data.table = FALSE)
      sample.names <- obs.df[[1]]
      metadata <- obs.df
      rownames(metadata) <- sample.names
      metadata <- metadata[,-1, drop = FALSE]
      cluster_vector <- metadata[,hierarchy_old[hierarchy[1]], drop = TRUE]
    } else {
      sample.names <- rownames(metadata)
    }
    
    cell_counts <- NULL
    norm_counts <- NULL
    t_cell_counts <- NULL
    t_norm_counts <- NULL
  } else {
    ##########################
    ## cell_counts
    print('... read and format cell counts.')
    if(!is.null(AIT.anndata)){
      cell_counts <- AIT.anndata$raw$X
      if(is.null(cell_counts)){
        cell_counts <- AIT.anndata$X
        if(is.null(cell_counts)){
          stop('chargeTaxonomy requires raw counts in the raw$X slot or the X to run.')
        }
      }
      rownames(cell_counts) <- rownames(AIT.anndata)
      colnames(cell_counts) <- colnames(AIT.anndata)
    } else {
      if(is.null(cell_counts)) stop('AIT.anndata or cell_counts must be provided to run chargeTaxonomy.')
    }
    all.sample.names = rownames(cell_counts)
    all.genes        = colnames(cell_counts)
    
    ##########################
    ## norm_counts
    print('... read and format normalized expression matrix.')
    if(!is.null(AIT.anndata)){
      norm_counts <- AIT.anndata$X
      if(is.null(norm_counts)){
        print('Creating normalized count matrix from cell_counts.')
        norm_counts = log2CPM_byRow(cell_counts)
      }
      if(max(norm_counts) > 100){
        print('Counts do not appear to be log-normalized; assuming X holds counts.')
        norm_counts = log2CPM_byRow(cell_counts)
      }
    } else {
      if(is.null(norm_counts))
        norm_counts = log2CPM_byRow(cell_counts)
    }
    rownames(norm_counts) <- all.sample.names
    colnames(norm_counts) <- all.genes
    
    #########################
    ## subset variables and define cluster_factor
    if(mean(keep_sample) < 1 ){
      print('... subsample count matrix')
      cell_counts  <- cell_counts[keep_sample,]
      norm_counts  <- norm_counts[keep_sample,]
      sample.names <- all.sample.names[keep_sample]
      cluster_vector <- cluster_vector[keep_sample]
    } else {
      sample.names <- all.sample.names
    }
    cluster_factor <- factor(cluster_vector, levels = all_clusters)
    names(cluster_factor) <- sample.names
    
    ##########################
    ## define transposes
    print('... transpose count matrix')
    t_cell_counts  <- as(Matrix::t(cell_counts),'dgCMatrix')
    rownames(t_cell_counts) <- all.genes
    colnames(t_cell_counts) <- sample.names
    
    print('... transpose logCPM matrix')
    t_norm_counts  <- as(Matrix::t(norm_counts),'dgCMatrix')
    rownames(t_norm_counts) <- all.genes
    colnames(t_norm_counts) <- sample.names
    
    #############################################################
    print('===== Building cluster-level statistics =====')
    
    count_n <- as.numeric(table(cluster_factor))
    names(count_n) <- all_clusters
    print('... sums')
    sums    <- round(get_cl_sums(t_cell_counts,cluster_factor))
    print('... counts')
    counts  <- round(get_cl_sums(t_cell_counts>0,cluster_factor))
    print('... props (calculated above)')
    props   <- t(t(counts)/as.numeric(count_n))
    print('... means')
    means   <- get_cl_means(t_norm_counts,cluster_factor)
    print('... sds')
    sds     <- sqrt(get_cl_vars(t_norm_counts,cluster_factor,means))
  }
    
  #############################################################
  print("===== Calculate distances and embeddings =====")
  
  ##########################
  # variable.genes
  print("... read in or calculate variable genes.")
  if(!is.null(AIT.anndata)){
    if(is.null(AIT.anndata$var$highly_variable_genes_standard)){
      betaScore      <- getBetaScore_fast(props[rowMaxs(props)>0.5,1:length(all_clusters)],returnScore=FALSE)
      betaScore      <- sort(betaScore)
      variable.genes <- names(betaScore)[1:min(1200,length(betaScore))]
    } else {
      variable.genes <- all.genes[AIT.anndata$var$highly_variable_genes_standard]
    }
  } else {
    if(is.null(variable.genes)) {
      betaScore      <- getBetaScore_fast(props[rowMaxs(props)>0.5,1:length(all_clusters)],returnScore=FALSE)
      betaScore      <- sort(betaScore)
      variable.genes <- names(betaScore)[1:min(1200,length(betaScore))]
    }
  }
  variable.genes <- intersect(variable.genes,all.genes)
  if(length(variable.genes)<50) error("<50 valid variable.genes, potentially due to misalignment between count matrix column names and variable.genes input. chargeTaxonomy cannot run with so few valid genes.")
  
 ##########################
  ## principal components
  print('... calculate principal components.')
  if(!is.null(stats.dir) && file.exists(file.path(stats.dir, 'rd_dat.csv'))){
    rd.df <- data.table::fread(file.path(stats.dir, 'rd_dat.csv'), data.table = FALSE)
    rn <- rd.df[[1]]
    rd.dat <- as.matrix(rd.df[,-1, drop = FALSE])
    rownames(rd.dat) <- rn
    if(!is.null(sample.names)) {
      rd.dat <- rd.dat[sample.names, , drop = FALSE]
    }
  } else {
    rd.dat = rd_PCA(t_norm_counts,
                    select.genes=variable.genes,
                    select.cells=sample.names,
                    max.pca = 50,
                    sampled.cells=sample.names,
                    th=0.5)
    rd.dat <- rd.dat$rd.dat
    rownames(rd.dat) <- sample.names
  }
  
  ##########################
  # embedding
  print('... read in or calculate UMAP')
  if(!is.null(stats.dir) && file.exists(file.path(stats.dir, 'umap.csv'))){
    umap.raw <- data.table::fread(file.path(stats.dir, 'umap.csv'), data.table = FALSE)
    rn <- umap.raw[[1]]
    umap.df <- as.data.frame(umap.raw[,-1, drop = FALSE])
    rownames(umap.df) <- rn
    umap.df <- umap.df[sample.names, , drop = FALSE]
  } else if(!is.null(AIT.anndata)){
    embedding <- NULL
    if(length(AIT.anndata$obsm)>0){
      embedding <- AIT.anndata$uns$default_embedding[[1]]
      if(is.null(embedding)) embedding = names(AIT.anndata$obsm)[1]
    }
    if(length(embedding)==1){
      umap.df <- AIT.anndata$obsm[[embedding]][sample.names,]
    } else {
      umap.df <- umap(rd.dat)$layout
    }
  } else {
    if(is.null(embedding)) {
      umap.df <- umap(rd.dat)$layout
    } else {
      umap.df <- embedding[sample.names,]
    }
  }
  umap.df <- as.data.frame(umap.df)
  rownames(umap.df) <- sample.names
  
  
  #############################################################
  print('===== Building statistics for the rest of the hierarchy =====')
  
  if(!is.null(stats.dir) && weight.by == 'cell'){
    print('... hierarchy stats pre-computed by Python, skipping.')
  } else if(!is.null(stats.dir) && weight.by == 'cluster') {
    print('Building sums:')
    sums   <- addHierarchyToStat(sums,hierarchy,cluster_info,'sum')
    print('Building counts:')
    counts <- addHierarchyToStat(counts,hierarchy,cluster_info,'sum')
    print('Building means:')
    means  <- addHierarchyToStat(means,hierarchy,cluster_info)
    print('Building sds:')
    sds    <- addHierarchyToStat(sds,hierarchy,cluster_info)
    print('Building props:')
    props  <- addHierarchyToStat(props,hierarchy,cluster_info)
    print('Building counts:')
    count_n2 <- addHierarchyToStat(rbind(count_n,count_n),hierarchy,cluster_info,'sum')
    count_n  <- setNames(as.numeric(count_n2[1,]),colnames(count_n2))
  } else if(weight.by != 'cluster'){
    # NOTE: THIS TREATS EACH CELL WITH THE SAME WEIGHT (e.g., BIGGER CLUSTER GET WEIGHTED HIGHER)
    if(weight.by != 'cell') warning("weight.by is not set to 'cell' or 'cluster'; defaulting to 'cell'.")
    for (i in 2:length(hierarchy)){
      print(paste(hierarchy[i],'........',i,'of',length(hierarchy)))
      cluster_vector2 = AIT.anndata$obs[,hierarchy[i]][keep_sample]
      all_clusters2   = unique(cluster_vector2)
      if(is.factor(all_clusters2))
        all_clusters2 = levels(all_clusters2)
      cluster_factor2 <- factor(cluster_vector2,levels=all_clusters2)
      names(cluster_factor2) <- names(cluster_factor)
      
      nm      <- names(count_n)
      count_n <- c(count_n,as.numeric(table(cluster_factor2)))
      names(count_n) <- c(nm,all_clusters2)
      print('... sums')
      sums    <- cbind(sums,round(get_cl_sums(t_cell_counts,cluster_factor2)))
      print('... counts')
      counts  <- cbind(counts,round(get_cl_sums(t_cell_counts>0,cluster_factor2)))
      print('... props')
      props   <- cbind(props,t(t(counts)/as.numeric(count_n)))
      print('... means')
      means2  <- get_cl_means(t_norm_counts,cluster_factor2)
      means   <- cbind(means,means2)
      print('... sds')
      sds     <- cbind(sds,sqrt(get_cl_vars(t_norm_counts,cluster_factor2,means2)))
    }
  } else {
    # NOTE: THIS SUMMARIZES EVERYTHING BY CLUSTER, TREATING EACH CLUSTER WITH THE SAME WEIGHT
    print('Building sums:')
    sums   <- addHierarchyToStat(sums,hierarchy,cluster_info,'sum')
    print('Building counts:')
    counts <- addHierarchyToStat(counts,hierarchy,cluster_info,'sum')
    print('Building means:')
    means  <- addHierarchyToStat(means,hierarchy,cluster_info)
    print('Building sds:')
    sds    <- addHierarchyToStat(sds,hierarchy,cluster_info)
    print('Building props:')
    props  <- addHierarchyToStat(props,hierarchy,cluster_info)
    print('Building counts:')
    count_n2 <- addHierarchyToStat(rbind(count_n,count_n),hierarchy,cluster_info,'sum')
    count_n  <- setNames(as.numeric(count_n2[1,]),colnames(count_n2))
  }
  
  
  #############################################################
  print("===== Build the constellation plots =====")
  
  constellation <- list()
  for (level in hierarchy){
    print(paste("... creating constellation for",level))
    cl.cl <- metadata[rownames(rd.dat), hierarchy_old[level]]
    cl.cl <- as.character(cl.cl)
    names(cl.cl) <- rownames(rd.dat)
    result = get_knn_graph(rd.dat, cl=cl.cl, k =50) 
    
    ## Select robust edges for plotting
    knn.cl.df = result$knn.cl.df 
    knn.cl.df = knn.cl.df %>% group_by(cl.from) %>% mutate(cl.from.rank = rank(-Freq))
    knn.cl.df = knn.cl.df %>% group_by(cl.to) %>% mutate(cl.to.rank = rank(-Freq))
    select.knn.cl.df = with(knn.cl.df, knn.cl.df[odds > 1 & pval.log < log(1/100000) & (frac > 0.1 | frac > 0.03 & Freq > 100) & (cl.from.rank < 4| cl.to.rank < 4),])
    
    ## Reorganize the cluster_info matrix for the information required for plotConstellation
    ## --- THIS SECTION NEEDS TO BE EDITED!
    prefix = level
    
    cl.center.df = as.data.frame(get_RD_cl_center(umap.df,cl.cl)) 
    cl.center.df$x <- cl.center.df$x + runif(length(cl.center.df$x), -0.000001, 0.000001)  # Small jitter
    cl.center.df$y <- cl.center.df$y + runif(length(cl.center.df$y), -0.000001, 0.000001)  # Small jitter
    
    ## Define cl.df
    types <- rownames(cl.center.df)
    cl.df <- cluster_info[,paste0(level,c("_id","_label","_color"))]
    cl.df <- cl.df[match(unique(cl.df[,1]),cl.df[,1]),]
    rownames(cl.df) <- cl.df[,2]
    cl.df <- cl.df[types,]
    cl.df$cluster_size <- count_n[types]
    
    ## Define cl.center.df (centroids for clusters)
    cl.center.df$cluster_id    <- cl.df[,1]  # id column
    cl.center.df$cluster_label <- rownames(cl.df)
    cl.center.df$cluster_color <- cl.df[,3]  # color column
    cl.center.df$cluster_size  <- cl.df$cluster_size
    rownames(cl.center.df)     <- rownames(cl.df)
    
    # Set cl as cluster_id since that was used to summarise the edges
    cl.center.df$cl <- cl.df[,1]  # The _id column
    
    # Convert labels to ids 
    select.knn.cl.df$cl.from <- cl.df[,1][match(select.knn.cl.df$cl.from,cl.df[,2])]
    select.knn.cl.df$cl.to   <- cl.df[,1][match(select.knn.cl.df$cl.to,cl.df[,2])]
    
    # Define the edges (I think that is what tmp.knn.cl.df does???)
    tmp.cl = cl.center.df$cluster_id
    tmp.knn.cl.df = select.knn.cl.df %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
    
    # Create the plot
    c.plot=try(plot_constellation(tmp.knn.cl.df, 
                                  cl.center.df=cl.center.df, 
                                  out.dir=NULL,
                                  prefix=prefix,
                                  node.label="cluster_label",
                                  exxageration=2,
                                  plot.parts=FALSE,
                                  return.list = T,
                                  node.dodge = F,
                                  label_repel = TRUE,
                                  label.size = 3,
                                  plot.height = 15,
                                  plot.width = 15,
                                  max_size = 5,
                                  enable_plotly = TRUE,
                                  plotly_labels_on_plot = TRUE))
    if(class(c.plot)[1]=="try-error") {
      constellation[[level]] <- plot_ly() %>%
        add_annotations(
          text = paste("No constellation diagram for",level),
          x = 0.5, y = 0.5,          # Center coordinates
          xref = "paper", yref = "paper", # Relative to plot area
          showarrow = FALSE,
          font = list(size = 36, color = "black") # Basic font styling
        ) %>%
        layout(
          xaxis = list(visible = FALSE), # Hide X-axis
          yaxis = list(visible = FALSE), # Hide Y-axis
          # Optional: make background transparent if embedding or don't want default gray
          plot_bgcolor = 'rgba(0,0,0,0)',
          paper_bgcolor = 'white'
        )
    } else {
      constellation[[level]] <- c.plot$constellation
    }
  } 
  
  #############################################################
  print("===== Save CHARGE file =====")
  save(
    hierarchy,
    cluster_info,
    count_n,
    counts,
    sums,
    means,
    sds,
    props,
    constellation,
    file = charge.file.name
  )
  
}



#' Convenience wrapper: read a backed h5ad plus Python-precomputed stats and write CHARGE.RData
#' 
chargeTaxonomyFromStats <- function(AIT.file,
                                    stats.dir,
                                    charge.file.name = 'CHARGE.RData',
                                    subsample = 100000000,
                                    weight.by = 'cell',
                                    seed = 42) {
  library(anndata)
  ad <- anndata::read_h5ad(AIT.file, backed = 'r')
  invisible(chargeTaxonomy(
    AIT.anndata = ad,
    stats.dir = stats.dir,
    subsample = subsample,
    weight.by = weight.by,
    charge.file.name = charge.file.name,
    seed = seed
  ))
}

#' Save an RData object for use with CHARGE app
#' 
#' This implementation of chargeTaxonomy uses python for all steps requiring reading data from the the X or raw.X slots (e.g., computing cluster statistics) and therefore should be compatible with large data sets.  **Note that python needs to be properly install and referenced in your working environment in addition to the R docker environment below.**
#' 
#' Write an RData object with the following variables
#'  - counts_n - vector of total number of cells per cluster
#'  - counts - cell x cluster matrix where each value represents the number of cells in a cluster with at least one read for each gene 
#'  - props - cell x cluster matrix where each value represents the fraction of cells in a cluster with at least one read for each gene
#'  - sums - cell x cluster matrix where each value represents the total number of counts of that gene for that cluster 
#'  - means - cell x cluster matrix where each value represents the average counts of that gene for that cluster 
#'  - sds - cell x cluster matrix where each value represents the standard deviation of counts of that gene for that cluster
#'  - hierarchy - cell type hierarchy stored in uns$hierarchy, but converted to a corrected ordered character vector
#'  - cluster_info - a table of cluster information that encodes the hierarchy for building sunburst plots
#'  - constellation - a list of plotly objects holding the constellation diagrams for each levels of the hierarchy
#'  **NOTE: This function is tested on "docker://alleninst/scrattch:1.1.4.1". We strongly encourage using this docker environment.**
#'
#' @param AIT.file The location of a reference taxonomy anndata object in AIT format.  Note that for this function the h5ad file must be local.
#' @param stats.dir A directory where python should write precomputed statistics from chargeTaxonomyHybrid for reading back in. The default is typically fine as thse statistics are not needed if the function runs properly.
#' @param python.script The default (NULL) looks in the R library for the correct script. In almost all cases, this should not be changed.
#' @param charge.file.name File name (and path) to write CHARGE file to
#' @param subsample How much subsampling should be done before running statistics. Default (recommended!) is none, but subsampling can be done to speed up calculations.
#' @param weight.by Should statistics for higher levels of the hierarchy be weighted by "cell" (default), whereby each statistic is recalculated on all the cells in a given group for each level, or by "cluster", whereby statistics are averaged, counting each item from the first level of the hierarchy (e.g., cluster) evenly.
#' @param seed The seed to use for reproducibility
#' @param python.exe What is the name of the executable for python (Default is 'python3')
#' 
#' @import Matrix
#' 
#' @examples
#' \dontrun{
#'   AIT.url  <- "https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SMART_seq_08082025.h5ad"
#'   AIT.file <- "Human_MTG_SMART_seq_08082025.h5ad"
#'   
#'   if (!file.exists(AIT.file)) {
#'     options(timeout = 9000)
#'     download.file(AIT.url, destfile = AIT.file, mode = "wb", timeout = 9000)
#'   }
#'   
#'   chargeTaxonomyHybrid(AIT.file)
#' }
#'
#' @export

chargeTaxonomyHybrid <- function(AIT.file,
                                 stats.dir = file.path(tempdir(),'charge_stats'),
                                 python.script = NULL,
                                 charge.file.name = 'CHARGE.RData',
                                 subsample = 100000000,
                                 weight.by = 'cell',
                                 seed = 42,
                                 python.exe = 'python3') {
  
  # 1) Resolve python script location
  if (is.null(python.script)) {
    # First try installed package location (inst/python -> python/)
    pkg_script <- system.file("python", "h5ad_to_charge_stats.py", package = "CHARGE_taxonomy")
    
    if (!is.null(pkg_script) && nzchar(pkg_script) && file.exists(pkg_script)) {
      python.script <- pkg_script
    } else {
      # Fallback for dev / source usage
      python.script <- file.path(getwd(), "h5ad_to_charge_stats.py")
    }
  }
  
  if (!file.exists(python.script)) {
    stop(
      "Could not find h5ad_to_charge_stats.py.\n",
      "Looked for:\n",
      "  1) system.file('python','h5ad_to_charge_stats.py', package='CHARGE_taxonomy')\n",
      "  2) ", python.script, "\n\n",
      "If developing locally, place h5ad_to_charge_stats.py in your working directory.\n",
      "If using the installed package, ensure it is located at inst/python/h5ad_to_charge_stats.py before install."
    )
  }
  
  # 2) Ensure output dir exists
  dir.create(stats.dir, recursive = TRUE, showWarnings = FALSE)
  
  # 3) Build command
  cmd <- sprintf(
    "%s %s --h5ad %s --outdir %s --subsample %s --weight_by %s --seed %s",
    shQuote(python.exe),
    shQuote(python.script),
    shQuote(AIT.file),
    shQuote(stats.dir),
    as.integer(subsample),
    shQuote(weight.by),
    as.integer(seed)
  )
  
  message("Running Python precompute: ", cmd)
  
  # 4) Run python
  status <- system(cmd)
  if (!identical(status, 0L)) stop("Python precompute failed with exit status ", status)
  
  # 5) Assemble CHARGE RData using backed h5ad + precomputed stats
  invisible(chargeTaxonomyFromStats(
    AIT.file = AIT.file,
    stats.dir = stats.dir,
    charge.file.name = charge.file.name,
    subsample = subsample,
    weight.by = weight.by,
    seed = seed
  ))
}
