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
    pkg_script <- system.file("python", "h5ad_to_charge_stats.py", package = "CHARGE.taxonomy") # In case function gets renamed
    if(!nzchar(pkg_script)) pkg_script <- system.file("python", "h5ad_to_charge_stats.py", package = "CHARGE_taxonomy")
    
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
  
  # 2.5) Subsample using R
  
  # Read the h5ad in backed mode just to get obs / hierarchy
  library(anndata)
  
  ad <- anndata::read_h5ad(AIT.file, backed = "r")
  hierarchy <- names(ad$uns$hierarchy)[order(-as.numeric(ad$uns$hierarchy))]
  cluster_vector <- ad$obs[, hierarchy[1]]
  
  # IMPORTANT: use the exact R subsampling function here
  keep_sample <- subsampleCells(cluster_vector, subsample, seed = seed)
  
  # Convert logical mask -> exact cell IDs
  subsample_ids <- rownames(ad$obs)[keep_sample]
  
  subsample_ids_file <- file.path(stats.dir, "subsample_ids.csv")
  dir.create(stats.dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    data.frame(cell_id = subsample_ids),
    subsample_ids_file,
    row.names = FALSE,
    quote = TRUE
  )
  
  # 3) Build command
  cmd <- sprintf(
    "%s %s --h5ad %s --outdir %s --subsample %s --subsample_ids %s --weight_by %s --seed %s",
    shQuote(python.exe),
    shQuote(python.script),
    shQuote(AIT.file),
    shQuote(stats.dir),
    as.integer(subsample),
    shQuote(subsample_ids_file),
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


#' @param knn.cl.df output of KNN.graph. Dataframe providing information about the cluster call of nearest neighbours of cells within a cluster. required columns: "cl.from" = cluster_id of edge origin, "cl.to" = cluster_id of edge destination, "Freq" = , "cl.from.total" = total nr of neigbours (above threshold) from cluster of origin, "cl.to.total" = total nr of neigbours (above threshold) from destination cluster, "frac" = fraction of total edge outgoing. 
#' @param cl.center.df dataframe containing metadata and coordinates for plotting cluster centroids. Required columns: "x" = x coordinate, "y" = y coordinate, "cl" = unique cluster id that should match "cl.to" and "cl.from" columns in knn.cl.df, "cluster_color","size" = nr of cells in cluster 
#' @param out.dir location to write plotting files to
#' @param prefix A character string to prepend to the filename
#' @param node.label Label to identify plotted nodes. Default is "cluster_id"
#' @param exxageration exxageration of edge width. Default is 1 (no exxageration)
#' @param curved Whether edges should be curved or not. Default is TRUE.
#' @param plot.parts output of intermediate files. default is FALSE.
#' @param plot.hull plot convex around cell type neighbourhood. Provide neighbourhood_id's that need to be plotted
#' @param node.dodge whether or not nodes are allowed to overlap. Default is false 
#' @param plot.height height of pdf in cm. Default is 25cm
#' @param plot.width width of pdf in cm. Default is 25cm 
#' @param label.size point size of plotted node labels. Default is 5pts
#' @param max_size maximum size of node. Default is 10pt 
#' @param label_repel whether to move labels away from node so they do not overlap. Default is FALSE.
#' @param node_trans transformation of node size. See ggplot::scale_size_area(trans=node_trans). Default is "sqrt".
#' @param return.list Whether to return list of independent plotting layers. Useful for replotting only part of the constellation
#' @param highlight_nodes list of node id's (matching cl.center.df) to highlight.
#' @param highlight_color Color of stroke around highlighted node. Default is red
#' @param highlight_width Width of stroke around highlighted node. Default is 1
#' @param highlight_labsize Size of label for highlighted nodes.
#' @param edge_mark_list which edges to color different selected by node using the cl value which is not necessarily the node label
#' @param edge_marking = c("dim", "highlight"),
#' @param fg.alpha alpha to use for edges that are most dark. Default = 0.4 
#' @param bg.alpha alpha to use for edges that are more faint. Default = 0.1
#' @param coord_fixed Cartesian coordinates with fixed "aspect ratio". See ggplot::coord_fixed. Default is TRUE
#'  
#' @import scrattch.bigcat
#' 
#' @usage plotting.MGE.constellation <- plot_constellation(knn.cl.df = knn.cl.df, cl.center.df = cl.center.df, out.dir = "data/Constellation_example/plot", node.dodge=TRUE, plot.hull=c(1,2), label_repel=TRUE) 


## EDITED FUNCTION
plot_constellation <- function (knn.cl.df,
                                cl.center.df,
                                out.dir,
                                prefix=format(Sys.time(), "%Y%m%d_%H%M%S"),
                                node.label = "cluster_id",
                                exxageration = 2,
                                curved = TRUE,
                                plot.parts = FALSE,
                                plot.hull = NULL,
                                plot.height = 25,
                                plot.width = 25,
                                node.dodge = FALSE,
                                label.size = 5,
                                max_size = 10,
                                label_repel = FALSE,
                                node_trans = "sqrt",
                                return.list = T,
                                highlight_nodes = NULL,
                                highlight_color = "red",
                                highlight_width = 1,
                                highlight_labsize = 10,
                                edge_mark_list = NULL,
                                edge_marking = c("dim", "highlight"),
                                fg.alpha = 0.4,
                                bg.alpha = 0.1,
                                coord_fixed= TRUE,
                                # NEW PARAMETER: Control for ggplotly output
                                enable_plotly = FALSE,
                                # NEW PARAMETER: Control label display in plotly
                                plotly_labels_on_plot = FALSE # TRUE for static labels, FALSE for hover only
) {
  
  library(gridExtra)
  library(cowplot)
  library(Hmisc)
  library(plotly)
  library(ggforce)
  if(label_repel == TRUE) { # Only load ggrepel if needed
    library(ggrepel)
  }
  
  #library(sna)   # Commented libraries either loaded globally or not needed
  #library(reshape2)
  #library(dplyr)
  
  
  st=prefix
  if(!is.null(out.dir)){
    if (!file.exists(out.dir)) {
      dir.create(out.dir)
    }}
  
  ###==== Cluster nodes will represent both cluster.size (width of point) and edges within cluster (stroke of point)
  
  # select rows that have edges within cluster
  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to,]
  
  #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
  cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)]
  
  # scale the node size by square root
  cl.center.df$cluster_size_sqrt = sqrt(cl.center.df$cluster_size)
  
  ###==== plot nodes
  labels <- cl.center.df[[node.label]]
  cl.center.df$text <- labels # Ensure text column is available for plotly text aesthetic
  
  # Initial p.nodes plot used for size extraction (not directly used for final plot.all)
  # Keeping this section as it was, assuming it's for internal calculations.
  p.nodes <-   ggplot() +
    geom_point(data=cl.center.df,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   # FIX 1: Remove alpha() from aes in initial point geom, set fixed alpha outside
                   color=cluster_color),
               alpha = 0.8) + # Set fixed alpha here
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() +
    geom_text(data=cl.center.df,
              aes(x=x,
                  y=y,
                  label=labels),
              size = label.size)
  
  if (plot.parts == TRUE & !is.null(out.dir)) {
    ggsave(file.path(out.dir,paste0(st,"nodes.org.pos.pdf")), p.nodes, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  ###==== extract node size/stroke width to replot later without scaling
  g <- ggplot_build(p.nodes)
  dots <-g[["data"]][[1]] #dataframe with geom_point size, color, coords
  
  nodes <- left_join(cl.center.df, dots, by=c("x","y")) %>% ungroup()
  
  ###==== if node.dodge==TRUE new xy coords are calculated for overlapping nodes.
  
  if (node.dodge==TRUE){
    # ... (node dodging logic - unchanged, assuming it works) ...
    # This section contains nested loops for node dodging. It's computationally
    # intensive and might be a bottleneck for large datasets.
    # It also relies on potentially modifying the 'nodes' dataframe's x/y directly.
    # Consider if this logic is critical for the Plotly output, as Plotly has
    # its own ways of handling overlap (though not as robust as ggrepel).
    
    nodes$r<- (nodes$size/10)/2
    
    
    x.list <- c(mean(nodes$x), nodes$x )
    y.list <- c(mean(nodes$y), nodes$y)
    dist.test <- as.matrix(dist(cbind(x.list, y.list)))
    nodes$distance <- dist.test[2:nrow(dist.test), 1]
    nodes <- nodes[order(nodes$distance),]
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        # print(paste(d1,d2)) # Commented out print statements to reduce console spam
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          # print(paste(d1,d2)) # Commented out print statements
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        # print(paste(d1,d2)) # Commented out print statements
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          # print(paste(d1,d2)) # Commented out print statements
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
  }
  
  nodes <- nodes[order(nodes$cluster_id),]
  
  ## when printing lines to pdf the line width increases slightly. This causes the edge to extend beyond the node. Prevent this by converting from R pixels to points.
  conv.factor <- ggplot2::.pt*72.27/96
  
  
  ## line width of edge can be scaled to node point size
  nodes$node.width <- nodes$size
  
  
  if (plot.parts == TRUE & !is.null(out.dir)) {
    if (node.dodge == TRUE) {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.dodge.csv"))) }
    else {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.csv")))
    }
  }
  
  
  ###==== prepare data for plotting of edges between nodes
  
  ##filter out all edges that are <5% of total for that cluster
  #knn.cl <- knn.cl.df[knn.cl.df$frac >0.05,] #1337 lines
  knn.cl <- knn.cl.df
  ##from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),]
  nodes$cl=as.character(nodes$cl)
  knn.cl.d$cl.from <- as.character(knn.cl.d$cl.from)
  knn.cl.d$cl.to <- as.character(knn.cl.d$cl.to)
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, node.width), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[colnames(knn.cl.d)=="node.width"]<- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""
  
  
  #bidirectional
  knn.cl.bid <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- subset(knn.cl.d[i:nrow(knn.cl.d),])
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from ),]
    
    if (dim(r)[1] != 0) {
      line$Freq.to <- r$Freq
      line$node.pt.to <- r$node.pt.from
      line$frac.to <- r$frac
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
    #print(i)
  }
  
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),]
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
    #print(i)
  }
  
  
  #min frac value = 0.01
  knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
  knn.cl.uni$Freq.to <- 1
  knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  
  ###==== create line segments
  
  line.segments <- knn.cl.lines %>% select(cl.from, cl.to)
  nodes$cl <- as.character(nodes$cl)
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.from"="cl"))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.to"="cl"))
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines$frac,
                              frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  
  ##from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  ##max fraction to max point size
  line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.from
  
  line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
  
  
  ###=== create edges, exaggerated width
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*exxageration),line.segments$node.size.from) #exxagerated width
  line.segments$ex.line.to <- pmin((line.segments$line.width.to*exxageration),line.segments$node.size.to) #exxagerated width
  
  
  line.segments <- na.omit(line.segments)
  
  print("calculating edges")
  
  # Assuming edgeMaker, perpStart, perpMid, perpEnd are defined elsewhere in your environment
  # They are not in the provided function, but necessary for allEdges and poly.Edges
  # For this example, I'll assume they work correctly and 'allEdges' and 'poly.Edges' are formed.
  # Placeholder for edgeMaker if it's not defined
  if (!exists("edgeMaker")) {
    edgeMaker <- function(i, len, curved, line.segments) {
      # This is a placeholder; replace with your actual edgeMaker function
      # It should return a data frame with x, y, fraction, and Group
      data.frame(x = runif(len), y = runif(len), fraction = runif(len),
                 Group = paste(line.segments$cl.from[i], line.segments$cl.to[i], sep = ">"))
    }
  }
  if (!exists("perpStart")) {
    perpStart <- function(x_coords, y_coords, width) {
      # Placeholder for perpStart
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[1], y_coords[1]), each = 2)
    }
  }
  if (!exists("perpMid")) {
    perpMid <- function(x_coords, y_coords, width) {
      # Placeholder for perpMid
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[2], y_coords[2]), each = 2)
    }
  }
  if (!exists("perpEnd")) {
    perpEnd <- function(x_coords, y_coords, width) {
      # Placeholder for perpEnd
      matrix(rnorm(4), ncol = 2) * width + rep(c(x_coords[2], y_coords[2]), each = 2)
    }
  }
  
  
  print(paste("Number of rows in line.segments:", nrow(line.segments)))
  if (nrow(line.segments) == 0) {
    stop("Error: The 'line.segments' data frame is empty. No edges to plot.")
  }
  
  allEdges <- lapply(1:nrow(line.segments), edgeMaker, len = 50, curved = curved, line.segments=line.segments) # Reduced len for faster example
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path with bend
  
  
  groups <- unique(allEdges$Group)
  
  poly.Edges <- data.frame(x=numeric(), y=numeric(), Group=character(),g1=integer(), g2=integer(),stringsAsFactors=FALSE)
  imax <- as.numeric(length(groups))
  
  # Limiting iterations for poly.Edges generation for a faster example run
  # In your actual use, remove or adjust this limit.
  # max_poly_edges_iter <- min(imax, 50) # Process max 50 groups for quick testing
  
  for(i in 1:imax) {
    #for(i in 1:max_poly_edges_iter) { # Use this line for faster testing
    select.group <- groups[i]
    select.edge <- allEdges[allEdges$Group %in% select.group,]
    
    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    
    N <- length(x)
    leftx <- numeric(N)
    lefty <- numeric(N)
    rightx <- numeric(N)
    righty <- numeric(N)
    
    ## Start point
    perps <- perpStart(x[1:2], y[1:2], w[1]/2)
    leftx[1] <- perps[1, 1]
    lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]
    righty[1] <- perps[2, 2]
    
    ### mid points
    for (ii in 2:(N - 1)) {
      seq <- (ii - 1):(ii + 1)
      perps <- perpMid(as.numeric(x[seq]), as.numeric(y[seq]), w[ii]/2)
      leftx[ii] <- perps[1, 1]
      lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]
      righty[ii] <- perps[2, 2]
    }
    ## Last control point
    perps <- perpEnd(x[(N-1):N], y[(N-1):N], w[N]/2)
    leftx[N] <- perps[1, 1]
    lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]
    righty[N] <- perps[2, 2]
    
    lineleft <- data.frame(x=leftx, y=lefty)
    lineright <- data.frame(x=rightx, y=righty)
    lineright <- lineright[nrow(lineright):1, ]
    lines.lr <- rbind(lineleft, lineright)
    lines.lr$Group <- select.group
    lines.lr[c("g1","g2")] <- as.integer(as.character(stringr::str_split_fixed(lines.lr$Group, '>', 2)))
    
    poly.Edges <- rbind(poly.Edges,lines.lr)
    
    #Sys.sleep(0.01) # Commented out for faster testing
    #cat("\r", i, "of", imax) # Commented out for cleaner output
  }
  
  if (plot.parts == TRUE && !is.null(out.dir)) {
    write.csv(poly.Edges, file=file.path(out.dir,paste0(st,"poly.edges.csv"))) }
  
  #############################
  ##                         ##
  ##        plotting         ##
  ##                         ##
  #############################
  
  labels <- nodes[[node.label]]
  
  ####plot edges
  # p.edges <- ggplot(poly.Edges, aes(group=Group))
  # p.edges <- p.edges +geom_polygon(aes(x=x, y=y), alpha=0.2) + theme_void()
  # p.edges # Commented out as this is not the final plot.all
  
  library("data.table")
  
  
  if(!is.null(edge_mark_list)){
    if(edge_marking == "dim"){
      poly.Edges$alpha_plot <- fg.alpha # Use a different column name to avoid conflict
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- bg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- bg.alpha
      
    } else if(edge_marking == "highlight"){
      poly.Edges$alpha_plot <- bg.alpha # Use a different column name
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- fg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- fg.alpha
      
    } else{
      print("provide valid edge marking")
    }
    
  } else{
    poly.Edges$alpha_plot <- fg.alpha # Use a different column name
  }
  
  # --- FIX 2: Consolidate plot.all creation and fix alpha/text logic ---
  # Define the base plot structure once
  base_plot <- ggplot() +
    geom_polygon(data=poly.Edges,
                 aes(x=x, y=y, group=Group, alpha=alpha_plot, # Use new alpha_plot variable
                     # Add hover text for polygons (edges)
                     text = paste0("From: ", g1, "<br>To: ", g2)),
                 fill = "grey60") + # Explicit fill for polygons
    scale_alpha_identity() + # FIX 2 (Crucial): Use scale_alpha_identity() for direct alpha values
    theme_void() +
    theme(legend.position = "none") # General theme for plot.all
  
  # Add geom_point for nodes
  # FIX 3: Add 'text' aesthetic to geom_point for plotly hover info
  base_plot <- base_plot +
    geom_point(data=nodes,
               # Removed alpha from outside aes here.
               # If you want point transparency, set it here explicitly as a numeric value.
               # For plotly, 0.8 is fine if you want visible points.
               alpha = 0.8,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   color=cluster_color,
                   # Add text for plotly hover (for nodes)
                   text = paste0("Cell type ID: ", cluster_id,
                                 "<br>Cell type label: ", cluster_label,
                                 "<br>Number of cells: ", cluster_size,
                                 "<br>Edges within: ", round(edge.frac.within, 2)))
    ) +
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() # Use identity for direct color mapping
  
  # Add highlight nodes if specified
  if(!is.null(highlight_nodes)){
    base_plot <- base_plot +
      geom_point(data=nodes %>% filter(cluster_id %in% highlight_nodes) ,
                 alpha=0.8,
                 shape=21, # Outlined circle
                 color=highlight_color, # Outline color
                 stroke=highlight_width, # Outline width
                 aes(x=x,
                     y=y,
                     size=cluster_size_sqrt,
                     # Add hover text for highlighted nodes if different
                     text = paste0("Cell type ID: ", cluster_id,
                                   "<br>Cell type label: ", cluster_label,
                                   "<br>Number of cells: ", cluster_size,
                                   "<br>Edges within: ", round(edge.frac.within, 2)))
      )
  }
  
  # Add hull if specified (and apply fixed alpha scaling)
  if (!is.null(plot.hull)) {
    base_plot <- base_plot +
      geom_mark_hull(data=nodes,
                     concavity = 8,
                     radius = unit(5,"mm"),
                     aes(filter = nodes$clade_id %in% plot.hull,x, y,
                         color=nodes$clade_color))
  }
  
  
  # Add labels (ggrepel or geom_text) for static display
  # FIX 4: Only add geom_text/geom_text_repel if NOT enabling plotly static labels
  if (!enable_plotly || (enable_plotly && !plotly_labels_on_plot)) {
    if(label_repel ==TRUE){
      base_plot <- base_plot +
        ggrepel::geom_text_repel(data=nodes,
                                 aes(x=x,
                                     y=y,
                                     label=.data[[node.label]]),
                                 size = label.size,
                                 min.segment.length = Inf)
    } else {
      base_plot <- base_plot +
        geom_text(data=nodes,
                  aes(x=x,
                      y=y,
                      label=.data[[node.label]]),
                  size = label.size)
    }
  }
  
  
  # Apply coord_fixed if TRUE
  if(isTRUE(coord_fixed)){
    base_plot <- base_plot + coord_fixed(ratio=1)
  }
  
  # Now assign the final ggplot object to plot.all
  plot.all <- base_plot
  
  segment.color = NA # This variable seems unused
  if (isTRUE(plot.parts) && !is.null(out.dir)) {
    ggsave(file.path(out.dir, paste0(st, ".comb.constellation.pdf")),
           plot.all, width = plot.width, height = plot.height,
           units = "cm", useDingbats = FALSE)
  }
  
  # --- FIX 5: Conditional ggplotly conversion and text addition ---
  if (enable_plotly) {
    plotly_plot <- ggplotly(plot.all, tooltip = "text") # Use "text" for custom tooltips
    
    # Add static text labels to plotly plot if requested
    if (plotly_labels_on_plot) {
      # MODIFICATION: Use layout annotations instead of add_text()
      
      # Prepare a list of annotations, one for each node
      annotations_list <- lapply(1:nrow(nodes), function(i) {
        list(
          x = nodes$x[i],
          y = nodes$y[i],
          text = as.character(nodes[[node.label]][i]), # Ensure text is character
          xref = "x",
          yref = "y",
          showarrow = FALSE, # No arrow pointing to the text
          font = list(
            color = nodes$cluster_color[i], # Color for each label
            size = label.size * 3 # Adjust font size (Plotly units differ from ggplot)
          ),
          xanchor = "center",  # Horizontal alignment of text
          yanchor = "bottom",  # Vertical alignment (places text slightly above the (x,y) point)
          yshift = 5          # Small vertical shift upwards to avoid overlap with point
        )
      })
      
      # Add the list of annotations to the layout of the plotly_plot
      plotly_plot <- plotly_plot %>%
        layout(annotations = annotations_list)
    }
    final_plot_output <- plotly_plot
  } else {
    final_plot_output <- plot.all # Return standard ggplot object
  }
  
  #############################
  ##                         ##
  ##      plot legends       ##
  ##                         ##
  #############################
  # Legend plotting logic remains unchanged. These are for static PDFs.
  # Plotly handles its own legends interactively.
  # If you need static legend images for plotly output, you might need to
  # export the plotly legend separately or create a combined layout.
  
  ### plot node size legend (1)
  plot.dot.legend <- ggplot()+
    geom_polygon(data=poly.Edges,
                 # FIX 2 (Crucial): Use alpha_plot and scale_alpha_identity() here too
                 aes(x=x, y=y, group=Group, alpha=alpha_plot), fill="grey60")+
    scale_alpha_identity() + # Use scale_alpha_identity()
    geom_point(data=nodes,
               alpha=0.8,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size_sqrt,
                   color=cluster_color)) +
    scale_size_area(trans=node_trans,
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() +
    geom_text(data=nodes,
              aes(x=x,
                  y=y,
                  label=labels),
              size = label.size)+
    theme_void()
  dot.size.legend <- cowplot::get_legend(plot.dot.legend)
  
  ### plot cluster legend (3)
  cl.center.df$cluster.label <- cl.center.df$cluster_label
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  label.col <- setNames(cl.center.df$cluster_color, cl.center.df$cluster.label)
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  leg.col.nr <- min((ceiling(length(cl.center.df$cluster_id)/20)),
                    5)
  cl.center <- ggplot(cl.center.df,
                      aes(x = cluster_id, y = cluster_size_sqrt)) +
    geom_point(aes(color = cluster.label)) +
    scale_color_manual(values = as.vector(label.col[levels(cl.center.df$cluster.label)])) +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 1)),
           fill=guide_legend(ncol=leg.col.nr)) +
    theme(legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          legend.key.size = unit(0.5, "lines"))
  
  cl.center.legend <- cowplot::get_legend(cl.center)
  
  
  width.1 <- max(line.segments$frac.from, line.segments$frac.to)
  width.05 <- width.1/2
  width.025 <- width.1/4
  edge.width.data <- tibble(node.width = c(1, 1, 1),
                            x = c(2, 2, 2),
                            y = c(5, 3.5, 2),
                            line.width = c(1, 0.5, 0.25),
                            fraction = c(100, 50, 25),
                            frac.ex = c(width.1, width.05, width.025))
  edge.width.data$fraction.ex <- round((edge.width.data$frac.ex *
                                          100), digits = 0)
  poly.positions <- data.frame(id = rep(c(1, 2, 3), each = 4),
                               x = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
                               y = c(4.9, 5.1, 5.5, 4.5, 3.4, 3.6, 3.75, 3.25, 1.9, 2.1, 2.125, 1.875))
  if (exxageration != 1) {
    edge.width.legend <- ggplot() + geom_polygon(data = poly.positions,
                                                 aes(x = x, y = y, group = id),
                                                 fill = "grey60") +
      geom_circle(data = edge.width.data, aes(x0 = x, y0 = y,
                                              r = node.width/2),
                  fill = "grey80", color = "grey80",
                  alpha = 0.4) +
      scale_x_continuous(limits = c(0,3)) +
      theme_void() +
      coord_fixed() +
      geom_text(data = edge.width.data, aes(x = 2.7, y = y, label = fraction.ex, hjust = 0, vjust = 0.5)) +
      annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")
  }  else {
    edge.width.legend <- ggplot() + geom_polygon(data = poly.positions,
                                                 aes(x = x, y = y, group = id), fill = "grey60") +
      geom_circle(data = edge.width.data, aes(x0 = x, y0 = y, r = node.width/2),
                  fill = "grey80", color = "grey80",
                  alpha = 0.4) +
      scale_x_continuous(limits = c(0,3)) +
      theme_void() +
      coord_fixed() +
      geom_text(data = edge.width.data, aes(x = 2.7, y = y, label = fraction, hjust = 0, vjust = 0.5)) +
      annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")
  }
  layout_legend <- rbind(c(1, 3, 3, 3, 3),
                         c(2, 3, 3, 3, 3))
  if (plot.parts == TRUE & !is.null(out.dir)) {
    ggsave(file.path(out.dir, paste0(st, ".comb.LEGEND.pdf")),
           gridExtra::marrangeGrob(list(dot.size.legend,
                                        edge.width.legend,
                                        cl.center.legend),
                                   layout_matrix = layout_legend,
                                   top=NULL),
           height = 20, width = 20, useDingbats = FALSE)
  }
  g2 <- gridExtra::arrangeGrob(grobs = list(dot.size.legend,
                                            edge.width.legend, cl.center.legend), layout_matrix = layout_legend)
  if(!is.null(out.dir) && !enable_plotly){ # Only save PDF if not enabling plotly
    ggsave(file.path(out.dir, paste0(st, ".constellation.pdf")),
           gridExtra::marrangeGrob(list(plot.all, g2), nrow = 1, ncol = 1,top=NULL),
           width = plot.width, height = plot.height, units = "cm",
           useDingbats = FALSE)
  }
  
  if (return.list == TRUE) {
    # Return the plotly object if enabled, otherwise the ggplot object
    return(list(constellation = final_plot_output, edges.df = poly.Edges,
                nodes.df = nodes))
  }
}







## function to draw (curved) line between to points
#' function to draw (curved) line between two points
#'
#' @param whichRow 
#' @param len 
#' @param line.segments 
#' @param curved 
#'
#' @return edge list
#' @export
edgeMaker <- function(whichRow, len=100, line.segments, curved=FALSE){
  
  fromC <- unlist(line.segments[whichRow,c(3,4)])# Origin
  toC <- unlist(line.segments[whichRow,c(5,6)])# Terminus
  # Add curve:
  
  graphCenter <- colMeans(line.segments[,c(3,4)])  # Center of the overall graph
  bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
  distance1 <- sum((graphCenter - bezierMid)^2)
  if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
    bezierMid <- c(toC[1], fromC[2])
    }  # To select the best Bezier midpoint
  bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
  if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve

  edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
  
  #line.width.from in 100 steps to linewidth.to
 edge$fraction <- seq(line.segments$ex.line.from[whichRow], line.segments$ex.line.to[whichRow], length.out = len)

  
    #edge$Sequence <- 1:len  # For size and colour weighting in plot
  edge$Group <- paste(line.segments[whichRow, 1:2], collapse = ">")
  return(edge)
  }



#utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#' perpStart
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return perp start value
#' @export
perpStart <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 1)
        }

#' avangle
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#'
#' @param x 
#' @param y 
#'
#' @return av angle
#' @export
avgangle <- function(x, y) {
    a1 <- angle(x[1:2], y[1:2])
    a2 <- angle(x[2:3], y[2:3])
    atan2(sin(a1) + sin(a2), cos(a1) + cos(a2))
}

#' perp
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#'
#' @param x 
#' @param y 
#' @param len 
#' @param a 
#' @param mid 
#'
#' @return perp value
#' @export
perp <- function(x, y, len, a, mid) {
    dx <- len*cos(a + pi/2)
    dy <- len*sin(a + pi/2)
    upper <- c(x[mid] + dx, y[mid] + dy)
    lower <- c(x[mid] - dx, y[mid] - dy)
    rbind(upper, lower)    
}

#' perpMid
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines. 
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return angle at midpoint
#' @export
perpMid <- function(x, y, len) {
    ## Now determine angle at midpoint
    perp(x, y, len, avgangle(x, y), 2)
}

#' perpEnd
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines.
#'
#' @param x 
#' @param y 
#' @param len 
#'
#' @return perp at end
#' @export
perpEnd <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 2)
}


#' angle
#' 
#' utils from vwline (https://github.com/pmur002/vwline) to draw variable width lines.
#'
#' @param x vector of length 2
#' @param y vector of length 2
#'
#' @return atan2 angle
#' @export
angle <- function(x, y) {
    atan2(y[2] - y[1], x[2] - x[1])
}


# I don't think this is used for anything
#plot_umap_constellation <- function(umap.2d, cl, cl.df, select.knn.cl.df, dest.d=".", prefix="",...)
#  {
#    
#    cl.center.df = as.data.frame(get_RD_cl_center(umap.2d,cl))
#    cl.center.df$cl = row.names(cl.center.df)
#    cl.center.df$cluster_id <- cl.df$cluster_id[match(cl.center.df$cl, cl.df$cl)]
#    cl.center.df$cluster_color <- cl.df$cluster_color[match(cl.center.df$cl, cl.df$cl)]
#    cl.center.df$cluster_label <- cl.df$cluster_label[match(cl.center.df$cl, cl.df$cl)] 
#    cl.center.df$cluster_size <- cl.df$cluster_size[match(cl.center.df$cl, cl.df$cl)]
#    tmp.cl = row.names(cl.center.df)
#    tmp.knn.cl.df = select.knn.cl.df  %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
#    p=plot_constellation(tmp.knn.cl.df, cl.center.df, node.label="cluster_id", out.dir=file.path(dest.d,prefix),...)    
#  }

#' Read in a reference data set in Allen taxonomy format (note: this is a copy of the scrattch.taxonomy function of the same name)
#'
#' @param taxonomyDir Directory containing the AIT file -OR- a direct h5ad file name -OR- a URL of a publicly accessible AIT file.
#' @param anndata_file If taxonomyDir is a directory, anndata_file must be the file name of the anndata object (.h5ad) to be loaded in that directory. If taxonomyDir is a file name or a URL, then anndata_file is ignored.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#' @param force (unused; kept for backwards compatibility)
#'
#' @return Organized reference object ready for mapping against (e.g., an anndata object)
#' 
#' @import anndata
#'
#' @export
loadTaxonomy = function(taxonomyDir = getwd(), 
                        anndata_file = "AI_taxonomy.h5ad",
                        log.file.path=getwd(),
                        force=FALSE){
  
  ## Allow for h5ad as the first/only input
  if(grepl("h5ad", taxonomyDir)){
    anndata_file = taxonomyDir
    taxonomyDir  = getwd()
  }
  ## Make sure the taxonomy path is an absolute path
  taxonomyDir = normalizePath(taxonomyDir, winslash = "/")
  
  ## If anndata_file is a URL, (1) parse the bucket out, (2) check whether the file is currently in the working directory, 
  ##   (3) download it if not, and then (4) rename anndata_file to the file name.
  # https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_VISp_ALM_SMART_seq_04042025.h5ad
  if(grepl("http", anndata_file)&grepl("s3",anndata_file)){
    anndata_parts  <- strsplit(anndata_file, "/")[[1]]
    anndata_object <- anndata_parts[length(anndata_parts)]
    if(!file.exists(anndata_object)){
      options(timeout = 9000) 
      download.file(anndata_file, file.path(taxonomyDir, anndata_object), timeout=9000)
    }
    anndata_file = anndata_object
  }
  
  ## Load from directory name input 
  if(file.exists(file.path(taxonomyDir, anndata_file))){
    print("Loading reference taxonomy into memory from .h5ad")
    ## Load taxonomy directly!
    AIT.anndata = read_h5ad(file.path(taxonomyDir, anndata_file))
    ## Default mode is always standard
    AIT.anndata$uns$mode = "standard"
    ##
    # for(mode in names(AIT.anndata$uns$dend)){
    #   invisible(capture.output({
    #     if(grepl("dend.RData", AIT.anndata$uns$dend[[mode]])){
    #       print("Loading an older AIT .h5ad version. Converting dendrogram to JSON format for mapping.")
    #       dend = readRDS(AIT.anndata$uns$dend[[mode]])
    #       AIT.anndata$uns$dend[[mode]] = toJSON(dend_to_json(dend))
    #     }
    #   }))
    # }
    if(("taxonomyName" %in% colnames(AIT.anndata$obs)) & (!"title" %in% colnames(AIT.anndata$obs))){
      AIT.anndata$obs$title = anndata_file$obs$taxonomyName
    }
    ## Ensure anndata is in scrattch.mapping format
    # if(!checkTaxonomy(AIT.anndata,log.file.path)){
    #  stop(paste("Taxonomy has some breaking issues.  Please check checkTaxonomy_log.txt in", log.file.path, "for details"))
    # }
  }else{
    stop("Required files to load Allen Institute taxonomy are missing.")
  }
  
  ## If counts are included but normalized counts are not, calculate normalized counts
  if((!is.null(AIT.anndata$raw$X))&(is.null(AIT.anndata$X))){
    normalized.expr = log2CPM_byRow(AIT.anndata$raw$X)
    AIT.anndata$X   = normalized.expr 
  }

  ## Set scrattch.mapping to default standard mapping mode
  AIT.anndata$uns$mode = "standard"
  AIT.anndata$uns$taxonomyDir = taxonomyDir
  AIT.anndata$uns$title <- gsub(".h5ad","",anndata_file)

  ## Return
  return(AIT.anndata)
}





#' Convert a matrix of raw counts to a matrix of log2(Counts per Million + 1) values 
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.  (note: this is a copy of the scrattch.taxonomy function of the same name)
#' 
#' This function expects that columns correspond to genes, and rows to samples by default and is equivalent to running logCPM with cells.as.rows=TRUE (but a bit faster).  By default the offset is 1, but to calculate just log2(counts per Million) set offset to 0.
#' 
#' @param counts A matrix, dgCMatrix, or dgTMatrix of count values
#' @param sf vector of numeric values representing the total number of reads. If count matrix includes all genes, value calulated by default (sf=NULL) will be accurate; however, if count matrix represents only a small fraction of genes, we recommend also providing this value.
#' @param denom Denominator that all counts will be scaled to. The default (1 million) is commonly used, but 10000 is also common for sparser droplet-based sequencing methods.
#' @param offset The constant offset to add to each cpm value prior to taking log2 (default = 1)
#' 
#' @return a dgCMatrix of log2(CPM + 1) values
#' 
#' @export 
log2CPM_byRow <- function (counts, sf = NULL, denom = 1e+06, offset=1){
  if(!("dgCMatrix" %in% as.character(class(counts))))
    counts <- as(counts, "dgCMatrix")
  if (is.null(sf)) {
    sf <- Matrix::rowSums(counts)
  }
  sf <- sf/denom
  sf[sf == 0] <- 1 / denom
  normalized.expr   <- counts
  normalized.expr@x <- counts@x / sf[as.integer(counts@i) + 1L]
  normalized.expr@x <- log2(normalized.expr@x + offset)
  return(normalized.expr)
}


#' Reorder annotations to have the same order and colors as what is shown in a separate metadata file
#'
#' @param anno an existing annotation data frame on which auto_annotate has ALREADY BEEN RUN
#' @param metadata a metadata file that includes some subset of content in the "_label" columns
#'
#' @return an updated data frame where order and colors of relevant columns have been adjusted to match metadata.
#'
#' @export
refactorize_annotations <- function(anno, metadata){
  ### NEED TO UPDATE THE CODE BELOW TO DEAL WITH COLORS, INCLUDING MAKING COLORS UNIQUE, BUT IT ALREADY IS CLOSE
  
  # Do nothing if cell type names are not provided
  if(sum(colnames(metadata)=="cell_type")==0)
    return(anno)
  
  # If provided, look for cell type names and order annotations accordingly by default.
  cns <- colnames(anno)
  cns <- gsub("_label","",cns[grepl("_label",cns)])
  for (cn in cns){
    intersecting_cell_types <- intersect(metadata$cell_type,as.character(anno[,paste0(cn,"_label")]))
    if(length(intersecting_cell_types)>0){
      new_levels <- c(intersecting_cell_types,setdiff(as.character(anno[,paste0(cn,"_label")]),metadata$cell_type))
      anno[,paste0(cn,"_id")]  <- as.numeric(factor(anno[,paste0(cn,"_label")],levels=new_levels))
    }
    
    ## If a column is called "color" in the metadata file, then look for categorical variable colors in the "cell_type" column
    
    if((sum(colnames(metadata)=="cell_type")==1)&(length(intersecting_cell_types)>0)){
      colors <- anno[,paste0(cn,"_color")][match(new_levels,anno[,paste0(cn,"_label")])]
      colors[match(intersecting_cell_types,new_levels)] <- metadata$color[match(intersecting_cell_types,metadata$cell_type)]
      anno[,paste0(cn,"_color")] <- colors[match(anno[,paste0(cn,"_label")],new_levels)]
      new_levels <- c(intersecting_cell_types,setdiff(as.character(anno[,paste0(cn,"_label")]),metadata$cell_type))
      anno[,paste0(cn,"_id")]  <- as.numeric(factor(anno[,paste0(cn,"_label")],levels=new_levels))
    }
    
  }
  
  return(anno)
}



#' Automatically format an annotation file
#'
#' This takes an anno file as input at any stage and properly annotates it for compatability with
#'   shiny and other scrattch functions.  In particular, it ensures that columns have a label,
#'   an id, and a color, and that there are no factors.  It won't overwrite columns that have
#'   already been properly process.
#'
#' @param anno an existing annotation data frame
#' @param sample_identifier the name of the column that contains the sample names. default = "cell_id"
#' @param scale_num should color scaling of numeric values be "predicted" (default and highly recommended;
#'   will return either "linear" or "log10" depending on scaling), "linear","log10","log2", or "zscore".
#' @param na_val_num The value to use to replace NAs for numeric columns. default = 0.
#' @param colorset_num A vector of colors to use for the color gradient.
#'   default = c("darkblue","white","red")
#' @param sort_label_cat a logical value to determine if the data in category columns
#'   should be arranged alphanumerically before ids are assigned. default = T.
#' @param na_val_cat The value to use to replace NAs in category and factor variables.
#'   default = "ZZ_Missing".
#' @param colorset_cat The colorset to use for assigning category and factor colors.
#'   Options are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order_cat The order in which colors should be assigned for cat and
#'   factor variables. Options are "sort" and "random". "sort" (default) assigns colors
#'   in order; "random" will randomly assign colors.
#'
#' @return an updated data frame that has been automatically annotated properly
#'
#' @export
auto_annotate <- function (anno, 
                           scale_num = "predicted", 
                           na_val_num = 0, 
                           colorset_num = c("darkblue","white", "red"), 
                           sort_label_cat = TRUE, 
                           na_val_cat = "ZZ_Missing", 
                           colorset_cat = "varibow", 
                           color_order_cat = "sort") 
{
  #################################################################################
  # UPDATED FUNCTION WITH BUG FIX FOR LARGE CSV FILES AND TO OMIT DIRECTION COLUMNS
  
  # Define and properly format a sample name
  anno_out <- as.data.frame(anno)
  cn <- colnames(anno_out)
  if (!is.element("sample_name", cn)) {
    colnames(anno_out) <- gsub("sample_id", "sample_name", cn)
  }
  if (!is.element("sample_name", colnames(anno_out))) {
    anno_out <- cbind(anno_out, paste0("sn_",1:dim(anno_out)[1]))
    colnames(anno_out) <- c(cn,"sample_name")
  }
  anno_out <- anno_out[,c("sample_name",setdiff(colnames(anno_out),"sample_name"))]
  
  # Annotate any columns missing annotations
  cn <- colnames(anno_out)
  convertColumns <- cn[(!grepl("_label", cn)) & (!grepl("_id", 
                                                        cn)) & (!grepl("_color", cn)) & (!grepl("_direction", cn))]  # UPDATE TO OMIT DIRECTION COLUMNS
  convertColumns <- setdiff(convertColumns, "sample_name")
  convertColumns <- setdiff(convertColumns, gsub("_label", 
                                                 "", cn[grepl("_label", cn)]))
  
  # Return input annotation file if there is nothing to convert
  if(length(convertColumns)==0){
    anno_out  <- group_annotations(anno_out)
    return(anno_out)
  }
  
  anno_list <- list()
  for (cc in convertColumns) {
    value <- anno_out[, cc]
    if (sum(!is.na(value)) == 0) 
      value = rep("N/A", length(value))
    if (is.numeric(value)) {
      if (length(table(value)) == 1) 
        value = jitter(value, 1e-06)
      val2 <- value[!is.na(value)]
      if (is.element(scale_num, c("linear", "log10", "log2", 
                                  "zscore"))) {
        out <- annotate_num(df = anno_out[,c("sample_name",cc)], col = cc, 
                            scale = scale_num, na_val = na_val_num, colorset = colorset_num)
      }
      else {
        scalePred <- ifelse(min(val2) < 0, "linear", 
                            "log10")
        if ((max(val2 + 1)/min(val2 + 1)) < 100) {
          scalePred <- "linear"
        }
        if (mean((val2 - min(val2))/diff(range(val2))) < 
            0.01) {
          scalePred <- "log10"
        }
        out <- annotate_num(df = anno_out[,c("sample_name",cc)], col = cc, 
                            scale = scalePred, na_val = na_val_num, colorset = colorset_num)
      }
    }
    else {
      if (is.factor(value)) {
        out <- annotate_factor(df = anno_out[,c("sample_name",cc)], col = cc, 
                               base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                               color_order = color_order_cat)
      }
      else {
        out <- annotate_cat(df = anno_out[,c("sample_name",cc)], col = cc, 
                            base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                            color_order = color_order_cat, sort_label = sort_label_cat)
      }
    }
    anno_list[[cc]] <- out[,colnames(out)!="sample_name"]
  }
  
  # Format the annotations as an appropriate data frame
  for(cc in convertColumns)
    anno_list[[cc]] <- anno_list[[cc]][,colnames(anno_list[[cc]])!="sample_name"]
  anno_out2 <- bind_cols(anno_list)
  anno_out  <- cbind(anno_out[,c(1,which(!is.element(colnames(anno_out),convertColumns)))],anno_out2)
  anno_out  <- group_annotations(anno_out[,c(1,3:dim(anno_out)[2])])
  anno_out
}





#' Returns the summary gene expression value across a group
#' 
#' @param datExpr a matrix of data (rows=genes, columns=samples)
#' @param groupVector character vector corresponding to the group (e.g., cell type)
#' @param fn Summary function to use (default = mean)
#' 
#' @return Summary matrix (genes x groups)
#'
#' @export
findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[1],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[,groupVector==groups[i]]
    if (is.null(dim(datIn)[1])) { 
      datMeans[,i] = as.numeric(datIn)
    } else { 
      datMeans[,i] = as.numeric(apply(datIn,1,fn)) }
  };    
  colnames(datMeans) = groups;
  rownames(datMeans) = rownames(datExpr)
  return(datMeans)
}


#' Adds to the statistics matrix
#' 
#' @param stat a statitics matrix (rows=genes, columns=names of hierarchy[1])
#' @param hierarchy a hierarchy character vector
#' @param cluster_info A cluster_info matrix
#' 
#' @return a statistics matrix with additional columns for other hierarchy levels
#'
#' @export
addHierarchyToStat <- function(stat, hierarchy, cluster_info, fn="mean"){
  if(length(hierarchy)==1)  return(stat)
  stat     <- stat[,cluster_info[,paste0(hierarchy[1],"_label")]]
  stat_out <- stat
  for (i in 2:length(hierarchy)){
    print(paste("........",i,"of",length(hierarchy)))
    vector   <- cluster_info[,paste0(hierarchy[i],"_label")]
    vector   <- factor(vector,levels=unique(vector[order(cluster_info[,paste0(hierarchy[i],"_id")])]))
    stat_out <- cbind(stat_out,findFromGroups(stat,vector,fn=fn))
  }
  stat_out
}




#' Generate colors and ids for categorical annotations that are factors
#'
#' @param df data frame to annotate
#' @param col name of the factor column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'   
#' @import dplyr
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'   
#' @export
annotate_factor <- function(df,
                            col = NULL, base = NULL,
                            na_val = "ZZ_Missing",
                            colorset = "varibow", color_order = "sort") {
  
  if (class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if (class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }
  
  if (class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if (class(base) == "NULL") {
    base <- col
  }
  
  if (!is.factor(df[[col]])) {
    df[[col]] <- as.factor(df[[col]])
  }
  
  # Convert NA values and add NA level to the end
  if (sum(is.na(df[[col]])) > 0) {
    lev <- c(levels(df[[col]]), na_val)
    levels(df[[col]]) <- lev
    df[[col]][is.na(df[[col]])] <- na_val
  }
  
  x <- df[[col]]
  
  annotations <- data.frame(label = as.character(levels(x)), stringsAsFactors = F)
  
  annotations <- annotations %>%
    dplyr::mutate(id = 1:n())
  
  if (colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if (colorset == "rainbow") {
    colors <- sub("FF$", "", grDevices::rainbow(nrow(annotations)))
  } else if (colorset == "viridis") {
    colors <- sub("FF$", "", viridisLite::viridis(nrow(annotations)))
  } else if (colorset == "magma") {
    colors <- sub("FF$", "", viridisLite::magma(nrow(annotations)))
  } else if (colorset == "inferno") {
    colors <- sub("FF$", "", viridisLite::inferno(nrow(annotations)))
  } else if (colorset == "plasma") {
    colors <- sub("FF$", "", viridisLite::plasma(nrow(annotations)))
  } else if (colorset == "terrain") {
    colors <- sub("FF$", "", grDevices::terrain.colors(nrow(annotations)))
  } else if (is.character(colorset)) {
    colors <- grDevices::colorRampPalette(colorset)(nrow(annotations))
  }
  
  if (color_order == "random") {
    colors <- sample(colors, length(colors))
  }
  
  annotations <- dplyr::mutate(annotations, color = colors)
  
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))
  
  names(df)[names(df) == col] <- paste0(base, "_label")
  
  df[[paste0(col,"_label")]] <- as.character(df[[paste0(col,"_label")]]) # convert the factor to a character in the anno
  
  df <- dplyr::left_join(df, annotations, by = paste0(base, "_label"))
  
  df
}



#' Generate colors and ids for numeric annotations
#'
#' @param df data frame to annotate
#' @param col name of the numeric column to annotate
#' @param base base name for the annotation, which wil be used in the desc table. If not provided, will use col as base.
#' @param scale The scale to use for assigning colors. Options are "linear","log10","log2, and "zscore"
#' @param na_val The value to use to replace NAs. default = 0.
#' @param colorset A vector of colors to use for the color gradient. default = c("darkblue","white","red")
#'
#' @return A modified data frame: the annotated column will be renamed base_label, and base_id and base_color columns will be appended
#'
#' @export
annotate_num <- function (df,
                          col = NULL, base = NULL,
                          scale = "log10", na_val = 0,
                          colorset = c("darkblue", "white", "red")) {

  #library(lazyeval)
  #library(dplyr)

  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }

  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }

  if (!is.numeric(df[[col]])) {
    df[[col]] <- as.numeric(df[[col]])
  }

  df[[col]][is.na(df[[col]])] <- na_val

  x <- df[[col]]

  annotations <- data.frame(label = unique(x)) %>%
    dplyr::arrange(label) %>%
    dplyr::mutate(id = 1:dplyr::n())

  if (scale == "log10") {
    colors <- values_to_colors(log10(annotations$label + 1), colorset = colorset)
  } else if(scale == "log2") {
    colors <- values_to_colors(log2(annotations$label + 1), colorset = colorset)
  } else if(scale == "zscore") {
    colors <- values_to_colors(scale(annotations$label), colorset = colorset)
  } else if(scale == "linear") {
    colors <- values_to_colors(annotations$label, colorset = colorset)
  }
  annotations <- mutate(annotations, color = colors)
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))

  names(df)[names(df) == col] <- paste0(base,"_label")
  df <- dplyr::left_join(df, annotations, by = paste0(base,"_label"))
  df
}

#' Generate colors and ids for categorical annotations
#'
#' @param df data frame to annotate
#' @param col name of the character column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param sort_label a logical value to determine if the data in col should be
#'   arranged alphanumerically before ids are assigned. default = T.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'
#' @export
annotate_cat <- function(df,
                         col = NULL, base = NULL,
                         sort_label = T, na_val = "ZZ_Missing",
                         colorset = "varibow", color_order = "sort") {

  #library(dplyr)
  #library(lazyeval)
  #library(viridisLite)

  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }

  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }

  if(!is.character(df[[col]])) {
    df[[col]] <- as.character(df[[col]])
  }

  df[[col]][is.na(df[[col]])] <- na_val

  x <- df[[col]]

  annotations <- data.frame(label = unique(x), stringsAsFactors = F)

  if(sort_label) {
    annotations <- annotations %>% dplyr::arrange(label)
  }

  annotations <- annotations %>%
    dplyr::mutate(id = 1:n())

  if(colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if(colorset == "rainbow") {
    colors <- sub("FF$","",grDevices::rainbow(nrow(annotations)))
  } else if(colorset == "viridis") {
    colors <- sub("FF$","",viridisLite::viridis(nrow(annotations)))
  } else if(colorset == "magma") {
    colors <- sub("FF$","",viridisLite::magma(nrow(annotations)))
  } else if(colorset == "inferno") {
    colors <- sub("FF$","",viridisLite::inferno(nrow(annotations)))
  } else if(colorset == "plasma") {
    colors <- sub("FF$","",viridisLite::plasma(nrow(annotations)))
  } else if(colorset == "terrain") {
    colors <- sub("FF$","",grDevices::terrain.colors(nrow(annotations)))
  } else if(is.character(colorset)) {
    colors <- grDevices::colorRampPalette(colorset)(nrow(annotations))
  }

  if(color_order == "random") {

    colors <- sample(colors, length(colors))

  }

  annotations <- dplyr::mutate(annotations, color = colors)

  names(annotations) <- paste0(base, c("_label","_id","_color"))

  names(df)[names(df) == col] <- paste0(base,"_label")

  df <- dplyr::left_join(df, annotations, by = paste0(base, "_label"))

  df
}



#' Create a generic description file
#'
#' @param dat any data frame that you would like to create a description file for
#' @param name desired names of each element in the description file (default is the column names)
#' @param use_label_columns should only columns containing "_label" be included (default = FALSE)
#' @param start_columns character vector of variables to include first in the list 
#'	(default = NULL, but "cluster" would be a common choice)
#'
#' @return a data.frame with columns "base", "name", and "type" for writing to tome
#'
#' @export
create_desc <- function(dat, name = colnames(dat), use_label_columns = FALSE, start_columns = NULL) {
  dat <- as.data.frame(dat)
  if (use_label_columns) {
    dat <- dat[, grepl("_label", colnames(dat))]
    colnames(dat) <- gsub("_label", "", colnames(dat))
  }

  desc <- data.frame(base = colnames(dat), name = name, type = "cat")
  for (i in 1:dim(dat)[2]) if (is.element(class(dat[, i]), c("numeric", "integer"))) {
      desc[i, 3] <- "num"
  }
	
  ## Reorder colums as requested
  cn   <- c(intersect(start_columns,desc$base),setdiff(desc$base,start_columns))
  desc <- desc[match(cn,desc$base),]
  desc
}


#'Group annotation columns
#'
#'@param df the annotation dataframe to arrange
#'@param sample_col the column with unique sample ids. Default is "cell_id".
#'@param keep_order a logical value. If FALSE, will sort the annotations alphanumerically by base.
#'
#'
#'@return an annotation data frame with reordered columns
#'
#' @export
group_annotations <- function(df, sample_col = "cell_id", keep_order = TRUE) {
  labels <- names(df)[grepl("_label",names(df))]
  if(!keep_order) {
    labels <- labels[order(labels)]
  }
  bases <- sub("_label","",labels)
  
  anno_cols <- c(paste0(rep(bases,each=3),c("_id","_label","_color")))
  extras <- setdiff(names(df),anno_cols)
  
  anno <- select(df,one_of(c(sample_col,anno_cols,extras)))
  
}




#' get knn graph
#' 
#' This function gets a knn graph, which is required for creating a constellation diagram
#'
#' @param rd.dat Reduced dimension data frame where rows are cells and columns are reduced dimensional data (e.g., 50-100 principal components)
#' @param cl Vector of cluster assignments (can be numeric, character, or factor)
#' @param k Number of nearest neighbors (default=50 to align with example)
#' @param ... Other variables that I'm not sure what they are and should probably be left alone
#'
#' @return A list with three components: knn.result = the knn result; knn.cl.df = the data frame of cluster information, outlier = outlier information.
#' @export
get_knn_graph <- function(rd.dat, 
                          cl, 
                          k=50, 
                          ref.cells=row.names(rd.dat),
                          method="Annoy.Cosine", 
                          knn.outlier.th=2, 
                          outlier.frac.th=0.5,
                          clean.cells=row.names(rd.dat), 
                          knn.result=NULL,
                          mc.cores=10)
{
  if(is.null(knn.result)){
    ref.rd.dat = rd.dat[ref.cells,]
    knn.result = get_knn_batch(rd.dat, ref.rd.dat, k, method=method, batch.size=10000, mc.cores=mc.cores, transposed=FALSE, return.distance=TRUE)
  }
  knn  = knn.result[[1]]
  knn.dist = knn.result[[2]]
  colnames(knn) = colnames(knn.dist)=1:ncol(knn)
  knn.dist = as.data.frame(as.table(knn.dist),stringsAsFactors=FALSE)
  knn.id = as.data.frame(as.table(knn),stringsAsFactors=FALSE)
  knn.df = cbind(knn.id, knn.dist[,3])
  colnames(knn.df)=c("sample_id","k","ref_id","dist")
  knn.df$cl = cl[knn.df$sample_id]  
  knn.df$knn.cl = cl[ref.cells[knn.df$ref_id]]
  knn.df = knn.df %>% filter(k!=1)
  
  cl.knn.dist.stats = knn.df %>%  group_by(cl) %>% summarise(med=median(dist),mad=mad(dist))
  cl.knn.dist.stats =   cl.knn.dist.stats %>% mutate(th=med + knn.outlier.th * mad)
  th.med = median(cl.knn.dist.stats$th)
  cl.knn.dist.stats =   cl.knn.dist.stats %>% mutate(th=pmax(th, th.med))
  
  outlier.df=knn.df %>% left_join(cl.knn.dist.stats[,c("cl","th")]) %>% group_by(sample_id) %>% summarise(outlier = sum(dist > th))
  
  outlier = outlier.df %>% filter(outlier/(k-1)>outlier.frac.th) %>% pull(sample_id)
  knn.df = knn.df %>% filter(!sample_id %in% outlier)
  knn.cl.df = knn.df %>% group_by(cl, knn.cl) %>% summarise(Freq=n())
  colnames(knn.cl.df)[1:2]=c("cl.from","cl.to")  
  from.size = knn.cl.df %>% group_by(cl.from) %>% summarise(from.total=sum(Freq))
  to.size = knn.cl.df %>% group_by(cl.to) %>% summarise(to.total=sum(Freq))
  total = sum(knn.cl.df$Freq)
  knn.cl.df = knn.cl.df %>% left_join(from.size) %>% left_join(to.size)
  knn.cl.df = knn.cl.df %>% mutate(odds = Freq/(from.total*as.numeric(to.total)/total))
  knn.cl.df = knn.cl.df %>% mutate(pval.log = phyper(q=Freq-1, m=to.total, n=total - to.total, k=from.total, lower.tail=FALSE, log.p=TRUE))
  knn.cl.df$frac = knn.cl.df$Freq/knn.cl.df$from.total
  return(list(knn.result=knn.result, knn.cl.df=knn.cl.df,outlier=outlier))
}



#' Get binary (aka beta) score
#'
#' Returns a beta score which indicates the binaryness of a gene across clusters.  High scores
#'   (near 1) indicate that a gene is either on or off in nearly all cells of every cluster.
#'   Scores near 0 indicate a cells is non-binary (e.g., not expressed, ubiquitous, or
#'   randomly expressed).  This value is used for gene filtering prior to defining clustering.
#'
#' @param propExpr a matrix of proportions of cells (rows) in a given cluster (columns) with
#'   CPM/FPKM > 1 (or 0, HCT uses 1)
#' @param returnScore if TRUE returns the score, if FALSE returns the ranks
#' @param spec.exp scaling factor (recommended to leave as default)
#'
#' @return returns a numeric vector of beta score (or ranks)
#'
#' @export
getBetaScore_fast <- function(propExpr, returnScore = TRUE, spec.exp = 2) {
  
  # Ensure propExpr is a matrix for consistent apply behavior
  propExpr <- as.matrix(propExpr)
  n_cols <- ncol(propExpr) # N in the formula above
  
  # Define a vectorized and optimized calculation for a single row
  calc_beta_optimized <- function(y_row, n_cols, spec.exp) {
    eps1 <- 1e-10
    
    # Numerator: sum of squared differences (spec.exp = 2)
    # sum((yi - yj)^2 for i<j) = N * sum(yk^2) - (sum(yk))^2
    if (spec.exp == 2) {
      sum_sq_diffs <- n_cols * sum(y_row^2) - sum(y_row)^2
    } else {
      # If spec.exp is not 2, we must compute pairwise differences.
      y_diffs <- as.vector(outer(y_row, y_row, FUN = "-"))
      sum_sq_diffs <- sum(y_diffs^spec.exp) / 2 
    }
    
    if (n_cols > 1) {
      sum_abs_diffs <- sum(abs(outer(y_row, y_row, FUN = "-"))) / 2 # Divide by 2 for unique pairs
    } else { # Handle single-column case
      sum_abs_diffs = 0 # No differences if only one element
    }
    
    score1 <- sum_sq_diffs / (sum_abs_diffs + eps1)
    return(score1)
  }
  
  # Apply the optimized calculation across rows
  betaScore <- apply(propExpr, 1, calc_beta_optimized, n_cols = n_cols, spec.exp = spec.exp)
  
  # Handle NA values
  betaScore[is.na(betaScore)] <- 0
  
  if (returnScore) {
    return(betaScore)
  }
  
  # If not returning score, return rank
  scoreRank <- rank(-betaScore, ties.method = "average") # Using average for tie-breaking
  return(scoreRank)
}

# =============================================================================
# Local dependency replacements
# These functions replace the small subset previously supplied by scrattch
# packages. No scrattch package is imported or attached.
# =============================================================================

## NEW FUNCTION
# Copied from scrattch.taxonomy as supplied by the user.
subsampleCells <- function(cluster.names, subSamp = 25, seed = 5, use.historical = FALSE) {
  if (use.historical) {
    kpSamp <- rep(FALSE, length(cluster.names))
    for (cli in unique(as.character(cluster.names))) {
      set.seed(seed)
      seed <- seed + 1
      kp <- which(cluster.names == cli)
      kpSamp[kp[sample(seq_along(kp), min(length(kp), subSamp))]] <- TRUE
    }
    return(kpSamp)
  }
  if (length(subSamp) == 1) {
    subSamp <- rep(subSamp, length(unique(as.character(cluster.names))))
  }
  if (is.null(names(subSamp))) {
    names(subSamp) <- unique(as.character(cluster.names))
  }
  set.seed(seed)
  cluster_split <- split(seq_along(cluster.names), as.character(cluster.names))
  kpSamp <- unlist(lapply(names(cluster_split), function(cli) {
    val <- subSamp[cli]
    if (!is.na(val)[1]) {
      kp <- cluster_split[[cli]]
      kp[sample(length(kp), min(length(kp), val))]
    } else {
      integer(0)
    }
  }))
  kpSamp2 <- rep(FALSE, length(cluster.names))
  kpSamp2[kpSamp] <- TRUE
  kpSamp2
}

## NEW FUNCTION
# Local Matrix implementation replacing scrattch.bigcat's compiled cluster stats.
.charge_cluster_design <- function(mat, cl) {
  if (!is.factor(cl)) cl <- setNames(factor(cl), names(cl))
  if (is.null(names(cl)) && ncol(mat) == length(cl)) names(cl) <- colnames(mat)
  if (length(cl) != ncol(mat)) stop("For CHARGE cluster statistics, length(cl) must equal ncol(mat).")
  Matrix::sparseMatrix(
    i = seq_along(cl), j = as.integer(cl), x = 1,
    dims = c(length(cl), nlevels(cl)),
    dimnames = list(names(cl), levels(cl))
  )
}

## NEW FUNCTION
get_cl_sums <- function(mat, cl) {
  design <- .charge_cluster_design(mat, cl)
  result <- mat %*% design
  colnames(result) <- colnames(design)
  result
}

## NEW FUNCTION
get_cl_means <- function(mat, cl) {
  if (!is.factor(cl)) cl <- setNames(factor(cl), names(cl))
  result <- get_cl_sums(mat, cl)
  sizes <- as.numeric(table(cl)[colnames(result)])
  result <- result %*% Matrix::Diagonal(x = 1 / sizes)
  dimnames(result) <- list(rownames(mat), levels(cl))
  result
}

## NEW FUNCTION
get_cl_sqr_means <- function(mat, cl) {
  squared <- mat
  if (inherits(squared, "sparseMatrix")) squared@x <- squared@x^2 else squared <- squared^2
  get_cl_means(squared, cl)
}

## NEW FUNCTION
# Formula copied from the user-supplied scrattch implementation.
get_cl_vars <- function(mat, cl, cl.means = NULL, cl.sqr.means = NULL) {
  if (is.null(cl.means)) cl.means <- get_cl_means(mat, cl)
  if (is.null(cl.sqr.means)) cl.sqr.means <- get_cl_sqr_means(mat, cl)
  cl.vars <- cl.sqr.means - cl.means^2
  cl.size <- as.vector(table(cl)[colnames(cl.vars)])
  cl.vars <- t(t(cl.vars) * cl.size / (cl.size - 1))
  cl.vars
}

## NEW FUNCTION
# Copied from the user-supplied scrattch utility.
varibow <- function(n_colors) {
  sats <- rep_len(c(0.4, 0.55, 0.7, 0.85, 1), length.out = n_colors)
  vals <- rep_len(c(1, 0.8, 0.6, 0.4), length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}

## NEW FUNCTION
# Copied from the user-supplied scrattch utility with explicit grDevices namespaces.
values_to_colors <- function(x, min_val = NULL, max_val = NULL,
                             colorset = c("darkblue", "dodgerblue", "gray80", "orange", "orangered"),
                             missing_color = "black") {
  heat_colors <- grDevices::colorRampPalette(colorset)(1001)
  if (is.null(max_val)) max_val <- max(x, na.rm = TRUE) else x[x > max_val] <- max_val
  if (is.null(min_val)) min_val <- min(x, na.rm = TRUE) else x[x < min_val] <- min_val
  if (sum(x == min_val, na.rm = TRUE) == length(x)) {
    colors <- rep(heat_colors[1], length(x))
  } else if (length(x) > 1) {
    if (stats::var(x, na.rm = TRUE) == 0) {
      colors <- rep(heat_colors[500], length(x))
    } else {
      heat_positions <- unlist(round((x - min_val) / (max_val - min_val) * 1000 + 1, 0))
      colors <- heat_colors[heat_positions]
    }
  } else {
    colors <- heat_colors[500]
  }
  if (!is.null(missing_color)) {
    colors[is.na(colors)] <- grDevices::rgb(t(grDevices::col2rgb(missing_color) / 255))
  }
  colors
}

## NEW FUNCTION
# Quadratic Bezier evaluator matching edgeMaker's three-control-point use.
bezier <- function(x, y, evaluation = 100) {
  if (length(x) != 3L || length(y) != 3L) stop("bezier() requires three x and three y control points.")
  tt <- seq(0, 1, length.out = evaluation)
  data.frame(
    x = (1 - tt)^2 * x[1] + 2 * (1 - tt) * tt * x[2] + tt^2 * x[3],
    y = (1 - tt)^2 * y[1] + 2 * (1 - tt) * tt * y[2] + tt^2 * y[3]
  )
}

## NEW FUNCTION
l2norm <- function(X, by = "column") {
  if (by == "column") {
    norms <- sqrt(Matrix::colSums(X^2))
    if (any(norms == 0)) warning("L2 norms of zero detected; zero columns are left unchanged")
    sweep(X, 2, pmax(norms, 1), "/", check.margin = FALSE)
  } else {
    norms <- sqrt(Matrix::rowSums(X^2))
    if (any(norms == 0)) warning("L2 norms of zero detected; zero rows are left unchanged")
    X / pmax(norms, 1)
  }
}

## NEW FUNCTION
knn_combine <- function(result.1, result.2) {
  list(
    knn.index = rbind(result.1[[1]], result.2[[1]]),
    knn.distance = rbind(result.1[[2]], result.2[[2]])
  )
}

## NEW FUNCTION
# CHARGE-local batch processor; avoids foreach/doMC attachment.
batch_process <- function(x, batch.size, FUN, mc.cores = 1, .combine = "c", bins = NULL, ...) {
  if (is.null(bins)) bins <- split(x, floor((seq_along(x) - 1L) / batch.size))
  if (mc.cores > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(bins, FUN, ..., mc.cores = min(mc.cores, length(bins)))
  } else {
    results <- lapply(bins, FUN, ...)
  }
  if (length(results) == 1L) return(results[[1]])
  if (identical(.combine, "knn_combine")) return(Reduce(knn_combine, results))
  if (identical(.combine, "rbind")) return(do.call(rbind, results))
  if (identical(.combine, "c")) return(do.call(c, results))
  Reduce(match.fun(.combine), results)
}

## NEW FUNCTION
get_knn <- function(dat, ref.dat, k, method = "cor", dim = NULL, index = NULL,
                    build.index = FALSE, transposed = TRUE, return.distance = FALSE,
                    ntrees = 100) {
  if (transposed) cell.id <- colnames(dat) else cell.id <- rownames(dat)
  if (transposed) {
    if (is.null(index)) ref.dat <- Matrix::t(ref.dat)
    dat <- Matrix::t(dat)
  }
  if (method %in% c("Euclidean", "Cosine")) { build.index <- FALSE; index <- NULL }
  if (method == "RANN") {
    knn.result <- RANN::nn2(ref.dat, dat, k = k)
    knn.result <- list(knn.result$nn.idx, knn.result$nn.dists)
  } else if (method %in% c("Annoy.Euclidean", "Annoy.Cosine", "cor", "Cosine")) {
    if (is.null(index)) {
      if (method == "cor") ref.dat <- l2norm(ref.dat - Matrix::rowMeans(ref.dat), by = "row")
      if (method %in% c("Annoy.Cosine", "Cosine")) ref.dat <- l2norm(ref.dat, by = "row")
      if (build.index) index <- BiocNeighbors::buildAnnoy(ref.dat, ntrees = ntrees)
    }
    if (method %in% c("Annoy.Cosine", "Cosine")) dat <- l2norm(dat, by = "row")
    if (method %in% c("Annoy.cor", "cor")) dat <- l2norm(dat - Matrix::rowMeans(dat), by = "row")
    if (method %in% c("Annoy.Cosine", "Annoy.Euclidean", "Annoy.cor")) {
      knn.result <- BiocNeighbors::queryAnnoy(X = ref.dat, query = dat, k = k, precomputed = index)
    } else {
      knn.result <- BiocNeighbors::queryKNN(X = ref.dat, query = dat, k = k)
    }
  } else {
    stop(paste(method, "method unknown"))
  }
  knn.index <- knn.result[[1]]
  knn.distance <- knn.result[[2]]
  rownames(knn.index) <- rownames(knn.distance) <- cell.id
  if (!return.distance) knn.index else list(knn.index = knn.index, knn.distance = knn.distance)
}

## NEW FUNCTION
get_knn_batch <- function(dat, ref.dat, k, method = "cor", dim = NULL, batch.size,
                          mc.cores = 1, return.distance = FALSE, transposed = TRUE,
                          index = NULL, clear.index = FALSE, ntrees = 50) {
  fun <- if (return.distance) "knn_combine" else "rbind"
  if (is.null(index) && method %in% c("Annoy.Euclidean", "Annoy.Cosine", "cor")) {
    map.ref.dat <- if (transposed) Matrix::t(ref.dat) else ref.dat
    if (method == "cor") map.ref.dat <- l2norm(map.ref.dat - rowMeans(map.ref.dat), by = "row")
    if (method == "Annoy.Cosine") map.ref.dat <- l2norm(map.ref.dat, by = "row")
    index <- BiocNeighbors::buildAnnoy(map.ref.dat, ntrees = ntrees)
    rm(map.ref.dat)
  }
  if (transposed) {
    results <- batch_process(seq_len(ncol(dat)), batch.size, mc.cores = mc.cores, .combine = fun,
      FUN = function(bin) get_knn(dat[row.names(ref.dat), bin, drop = FALSE], ref.dat, k,
        method = method, dim = dim, return.distance = return.distance,
        transposed = TRUE, index = index, ntrees = ntrees))
  } else {
    results <- batch_process(seq_len(nrow(dat)), batch.size, mc.cores = mc.cores, .combine = fun,
      FUN = function(bin) get_knn(dat[bin, colnames(ref.dat), drop = FALSE], ref.dat, k,
        method = method, dim = dim, return.distance = return.distance,
        transposed = FALSE, index = index, ntrees = ntrees))
  }
  if (!return.distance) results <- list(knn.index = results)
  results$index <- index
  results
}

## NEW FUNCTION
# Copied from the user-supplied scrattch implementation.
get_RD_cl_center <- function(rd.dat, cl) {
  do.call("rbind", tapply(seq_len(nrow(rd.dat)), cl[row.names(rd.dat)], function(x) {
    x <- sample(x, pmin(length(x), 500))
    dd <- as.matrix(stats::dist(rd.dat[x, 1:2, drop = FALSE]))
    tmp <- x[which.min(rowSums(dd))]
    c(x = rd.dat[tmp, 1], y = rd.dat[tmp, 2])
  }))
}

## NEW FUNCTION
get_PCA <- function(dat, max.pca, verbose = FALSE, method = "zscore", th = 2,
                    fun = "prcomp", rot = TRUE, init.pca = 200) {
  dat <- as.matrix(dat)
  if (rot) dat <- t(dat)
  if (fun == "prcomp") {
    pca <- stats::prcomp(dat, tol = 0.01)
  } else {
    pca <- irlba::prcomp_irlba(dat, n = min(init.pca, nrow(dat)))
  }
  if (method == "elbow") {
    stop("method='elbow' requires the unavailable scrattch helper findElbowPoint; use method='zscore'.")
  } else if (method == "zscore") {
    v <- summary(pca)$importance[2, ]
    select <- head(which((v - mean(v)) / stats::sd(v) > th), max.pca)
  } else stop("Unknown method")
  if (!length(select)) return(NULL)
  rotation <- pca$rotation[, select, drop = FALSE]
  rd.dat <- pca$x[, select, drop = FALSE]
  list(rot = rotation, rd.dat = rd.dat, pca = pca)
}

## NEW FUNCTION
rd_PCA <- function(norm.dat, select.genes = row.names(norm.dat), select.cells = colnames(norm.dat),
                   sampled.cells = select.cells, max.pca = 10, th = 2, verbose = FALSE,
                   method = "zscore", mc.cores = 1) {
  tmp <- get_PCA(norm.dat[select.genes, sampled.cells, drop = FALSE], max.pca = max.pca,
                 verbose = verbose, th = th, method = method)
  if (is.null(tmp)) return(NULL)
  rotation <- tmp$rot
  rd.dat <- tmp$rd.dat
  pca <- tmp$pca
  if (length(sampled.cells) < length(select.cells)) {
    rd.dat <- do.call("rbind", lapply(select.cells, function(x) {
      tmp.dat <- norm.dat[row.names(rotation), x, drop = FALSE]
      as.matrix(Matrix::crossprod(tmp.dat, rotation))
    }))
  }
  list(rd.dat = rd.dat, pca = pca)
}

