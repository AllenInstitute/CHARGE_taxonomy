#' Read in a reference data set in Allen taxonomy format (note: this is a copy of the scrattch.taxonomy function of the same name)
#'
#' @param taxonomyDir Directory containing the AIT file -OR- a direct h5ad file name -OR- a URL of a publicly accessible AIT file.
#' @param anndata_file If taxonomyDir is a directory, anndata_file must be the file name of the anndata object (.h5ad) to be loaded in that directory. If taxonomyDir is a file name or a URL, then anndata_file is ignored.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#'
#' @return Organized reference object ready for mapping against.
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
  normalized.expr   <- counts
  normalized.expr@x <- counts@x/sf[as.integer(counts@i + offset)]
  normalized.expr@x <- log2(normalized.expr@x+1)
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
  
  # library(dplyr)
  # library(lazyeval)
  # library(viridisLite)
  
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
