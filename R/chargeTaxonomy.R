# This function is tested on scrattch docker 1.2

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
#'  - cluster_info
#'  NOTE: proportions are not provided because they can be calculated trivaially from the above variables
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param check.taxonomy Should the function check if the input variable is a valid AIT taxonomy (default = TRUE)
#' @param subsample How much subsampling should be done before running statistics. Default (recommended!) is none, but subsampling can be done to speed up calculations.
#' @param charge.file.name File name (and path) to write CHARGE file to
#' @param seed The seed to use for reproducibility
#' 
#' @import scrattch.taxonomy
#' @import rsvd
#' @import Matrix
#' 
#' @examples
#' \dontrun{
#'   AIT.file    <- "HMBA_BG_Macaque_AIT.h5ad" #"Macaque_HMBA_basalganglia_AIT_pre-print.h5ad"
#'   AIT.anndata <- loadTaxonomy(AIT.file)
#'   #raw = AIT.anndata$raw$X
#'   #raw@x <- round(2^(raw@x)-1)
#'   #AIT.anndata$raw$X@x <- raw@x  # This does not work
#'   chargeTaxonomy(AIT.anndata)
#' }
#'
#' @export
chargeTaxonomy <- function(AIT.anndata,
                           check.taxonomy = TRUE,
                           subsample      = 100000000,
                           charge.file.name = paste0(AIT.anndata$title,"_CHARGE.RData",
						   seed = 42)
){
  #############################################################
  print("===== Checking taxonomy, if requested =====")
  if(check.taxonomy){
    AIT.anndata = checkTaxonomy(AIT.anndata)
    if(!AIT.anndata$uns$valid) warning("Taxonomy does not pass taxonomy checks, potentially causing chargeTaxonomy to fail. Check the log file if needed.") 
  } else {
    warning("Taxonomy is not being checked. If chargeTaxonomy fails, we recommend running checkTaxonomy to ensure valid AIT file.")
  }
  
  cell_counts <- AIT.anndata$raw$X
  if(is.null(cell_counts)){
    cell_counts <- AIT.anndata$X
    if(is.null(cell_counts)){
      error("chargeTaxonomy required raw counts in the raw$X slot or the X to run.")
    }
  }
  
  norm_counts <- AIT.anndata$X
  if(is.null(norm_counts)){
    print("Creating normalized count matrix from cell_counts.")
    norm_counts = log2CPM_byRow(cell_counts)
  }
  if(max(norm_counts)>100){
    print("Counts do not appear to be log-normalized; assuming X holds counts.")
    norm_counts = log2CPM_byRow(cell_counts)
  }

  #############################################################
  print("===== Collecting required cluster information =====")
  # hierarchy
  hierarchy = names(AIT.anndata$uns$hierarchy)[order(-as.numeric(AIT.anndata$uns$hierarchy))]
  
  # cluster_info
  cluster_vector = AIT.anndata$obs[,hierarchy[1]]
  all_clusters   = unique(cluster_vector)
  if(is.factor(all_clusters)) 
    all_clusters = levels(all_clusters)
  if(!is.null(AIT.anndata$uns$cluster_info)){
    cluster_info = AIT.anndata$uns$cluster_info
  } else {
    cluster_info = AIT.anndata$obs
  }
  cluster_info = cluster_info[match(all_clusters,cluster_info[,hierarchy[1]]),hierarchy]

  # Rename anything called "_id" or "_label" to avoid breaking shiny
  hierarchy_old <- hierarchy
  hierarchy <- gsub("_id","",hierarchy)
  hierarchy <- gsub("_label","",hierarchy)
  colnames(cluster_info) <- hierarchy
  names(hierarchy_old) <- hierarchy

  # Subsample vector
  set.seed(seed)
  keep_sample <- subsampleCells(cluster_vector,subsample)
  if(mean(keep_sample)<1 ){
    cell_counts    <- cell_counts[,keep_sample]
	norm_counts    <- norm_counts[,keep_sample]
    cluster_vector <- cluster_vector[keep_sample]
  }
  
  print("===== Building statistics =====")

  # Building sums, counts, and props from the COUNT MATRIX
  print("... transpose count matrix")
  t_cell_counts  <- as(Matrix::t(cell_counts),"dgCMatrix")
  cluster_factor <- factor(cluster_vector,levels=all_clusters)
  rownames(t_cell_counts) <- colnames(AIT.anndata)
  colnames(t_cell_counts) <- names(cluster_factor) <- rownames(AIT.anndata)[keep_sample]

  # Building means and sds from LOG2CPM
  print("... transpose logCPM matrix")
  t_norm_counts  <- as(Matrix::t(norm_counts),"dgCMatrix")
  rownames(t_norm_counts) <- colnames(AIT.anndata)
  colnames(t_norm_counts) <- names(cluster_factor) <- rownames(AIT.anndata)[keep_sample]
  

  count_n <- as.numeric(table(cluster_factor))
  names(count_n) <- all_clusters
  print("... sums")
  sums    <- round(get_cl_sums(t_cell_counts,cluster_factor))
  print("... counts")
  counts  <- round(get_cl_sums(t_cell_counts>0,cluster_factor))
  print("... props")
  props   <- t(t(counts)/as.numeric(count_n))
  print("... means")
  means   <- get_cl_means(t_norm_counts,cluster_factor)
  print("... sds")
  sds     <- sqrt(get_cl_vars(t_norm_counts,cluster_factor,means))


  #############################################################
  print("===== Add cluster ids and colors =====")
  ## Add ids and colors to the cluster information
  ## NOTE: I'LL NEED TO UPDATE THIS PART ONCE I SORT OUT HOW TO EMBED CLUSTER COLORS
  ## (AND CLUSTER ORDER?) WITHIN THE CLUSTER_INFO DATA FRAME
  cluster_info$sample_name = paste0("i",1:dim(cluster_info)[1])
  for (i in 1:dim(cluster_info)[2]){
    cluster_info[,i] <- factor(cluster_info[,i],levels = unique(cluster_info[,i]))
  }    
  cluster_info <- auto_annotate(cluster_info)
  
  
  #############################################################
  print("===== Building statistics for the rest of the hierarchy =====")
  # NOTE: THIS TREATS EACH CELL WITH THE SAME WEIGHT (e.g., BIGGER CLUSTER GET WEIGHTED HIGHER)
  for (i in 2:length(hierarchy)){
    print(paste(hierarchy[i],"........",i,"of",length(hierarchy)))
    cluster_vector2 = AIT.anndata$obs[,hierarchy[i]][keep_sample]
	all_clusters2   = unique(cluster_vector2)
    if(is.factor(all_clusters2)) 
      all_clusters2 = levels(all_clusters2)
	cluster_factor2 <- factor(cluster_vector2,levels=all_clusters2)
	names(cluster_factor2) <- names(cluster_factor)
	
    nm      <- names(count_n)
    count_n <- c(count_n,as.numeric(table(cluster_factor2)))
    names(count_n) <- c(nm,all_clusters2)
    print("... sums")
    sums    <- cbind(sums,round(get_cl_sums(t_cell_counts,cluster_factor2)))
    print("... counts")
    counts  <- cbind(counts,round(get_cl_sums(t_cell_counts>0,cluster_factor2)))
    print("... props")
    props   <- cbind(props,t(t(counts)/as.numeric(count_n)))
    print("... means")
	means2  <- get_cl_means(t_norm_counts,cluster_factor2)
    means   <- cbind(means,means2)
    print("... sds")
    sds     <- cbind(sds,sqrt(get_cl_vars(t_norm_counts,cluster_factor2,means2)))
  }
  
  # NOTE: THIS SUMMARIZES EVERYTHING BY CLUSTER, TREATING EACH CLUSTER WITH THE SAME WEIGHT
  # Uncomment this code and comment the code above if we want to do it this way
  #print("Building sums:")
  #sums   <- addHierarchyToStat(sums,hierarchy,cluster_info,"sum")
  #print("Building counts:")
  #counts <- addHierarchyToStat(counts,hierarchy,cluster_info,"sum")
  #print("Building means:")
  #means  <- addHierarchyToStat(means,hierarchy,cluster_info)
  #print("Building sds:")
  #sds    <- addHierarchyToStat(sds,hierarchy,cluster_info)
  #print("Building props:")
  #props  <- addHierarchyToStat(props,hierarchy,cluster_info)
  #print("Building counts:")
  #count_n2 <- addHierarchyToStat(rbind(count_n,count_n),hierarchy,cluster_info,"sum")
  #count_n  <- setNames(as.numeric(count_n2[1,]),colnames(count_n2))
  
 
  #############################################################
  print("===== Create input files for constellation plots =====")
  
  # Create a reduced dimensionality matrix (essentially creating principal components)
  # ---- We are starting from the normalized data here
  # ---- Finally, we take the top 1000 highly variable genes
  
  ## Get the top variable genes
  if(is.null(AIT.anndata$var$highly_variable_genes_standard)){
    betaScore <- getBetaScore_fast(props[rowMaxs(props)>0.5,1:length(all_clusters)],returnScore=FALSE)
    betaScore <- sort(betaScore)
    top.genes <- names(betaScore)[1:min(1000,length(betaScore))]
  } else {
    top.genes <- AIT.anndata$var$highly_variable_genes_standard
  }
  
  
  ## GET PCS (e.g., reduced dimension data frame)
  print("... calculate principal components")
  
  # OLD WAY TO DO IT
  #max_k_to_estimate   <- min(c(nrow(norm_counts),ncol(norm_counts),200))
  #pca_result          <- rpca(norm_counts[,top.genes], k = max_k_to_estimate, center = TRUE, scale = TRUE)
  #exp_variance_ratio  <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  #cumulative_variance <- cumsum(exp_variance_ratio)
  #variance_threshold  <- 0.9
  #n_components_chosen <- which.max(cumulative_variance >= variance_threshold)
  #rd.dat              <- pca_result$x[,1:n_components_chosen]
  #colnames(rd.dat)    <- paste0("PC",1:dim(rd.dat)[2])  
  
  
  ###############################
  # dimension reduction USING CINDY'S CODE

  select.markers = colnames(norm_counts)[top.genes]  # Replacing the super slow de.genes code for now.

  rd.dat = rd_PCA(t_norm_counts,
                select.genes=select.markers, 
                select.cells=rownames(norm_counts), 
                max.pca = 50,
                method = "elbow",
                sampled.cells=rownames(norm_counts), 
                th=0.5)

  rd.dat <- rd.dat$rd.dat
  # rd.dat is used as input for knn graph
  
  ## GET THE UMAP
  print("... read in or calculate UMAP")
  embedding <- NULL
  print(1)
  if(length(AIT.anndata$obsm)>0){
    embedding <- AIT.anndata$uns$default_embedding[[1]]
    if(is.null(embedding)) embedding = names(AIT.anndata$obsm)[1]
  }
  print(2)
  if(length(embedding)==1){
    umap.df <- AIT.anndata$obsm[[embedding]]
  } else {
    umap.df <- umap(rd.dat)$layout
    umap.df <- as.data.frame(umap.df)
  }
  print(3)
  rownames(umap.df) <- rownames(rd.dat)
  print(4)
  
  ##  NOW BUILD THE CONSTELLATION
  constellation <- list()
  for (level in hierarchy){
    print(paste("... creating constellation for",level))
    cl.cl  = AIT.anndata$obs[,hierarchy_old[level]]
    names(cl.cl) = rownames(rd.dat)
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
    cl.df$cluster_size <- count_n[types] #as.numeric(table(AIT.anndata$obs[,hierarchy[1]])[cl.df[,2]]) #count_n[types]
    
    ## Define cl.center.df (centroids for clusters)
    cl.center.df$cluster_id    <- cl.df[,1]  # id column
    cl.center.df$cluster_label <- rownames(cl.df)
    cl.center.df$cluster_color <- cl.df[,3]  # color column
    cl.center.df$cluster_size  <- cl.df$cluster_size
    rownames(cl.center.df)     <- rownames(cl.df)
    
    # set cl as cluster_id since that was used to summarise the edges
    #cl.center.df$cl = as.integer(as.character(row.names(cl.center.df) ))
    cl.center.df$cl <- cl.df[,1]  # The _id column
    
    # Convert labels to ids 
    select.knn.cl.df$cl.from <- cl.df[,1][match(select.knn.cl.df$cl.from,cl.df[,2])]
    select.knn.cl.df$cl.to   <- cl.df[,1][match(select.knn.cl.df$cl.to,cl.df[,2])]
    
    # Define the edges (I think that is what tmp.knn.cl.df does???)
    tmp.cl = cl.center.df$cluster_id
    tmp.knn.cl.df = select.knn.cl.df %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
    
    # Create the plot
	#cl.center.df$cl <- rownames(cl.center.df)
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
  print("===== Save file =====")
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

