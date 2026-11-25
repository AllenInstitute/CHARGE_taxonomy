
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
  
  # The point layer has one row per input node. Joining on coordinates would
  # duplicate nodes when two types have exactly the same centroid.
  stopifnot(nrow(dots) == nrow(cl.center.df),
            isTRUE(all.equal(dots$x, cl.center.df$x, check.attributes = FALSE)),
            isTRUE(all.equal(dots$y, cl.center.df$y, check.attributes = FALSE)))
  nodes <- cl.center.df
  nodes$size <- dots$size
  
  ###==== if node.dodge==TRUE new xy coords are calculated for overlapping nodes.
  
  if (node.dodge == TRUE && nrow(nodes) > 1L){
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
  
  
  ###==== Prepare the retained between-type edges (zero is valid).
  nodes$cl <- as.character(nodes$cl)
  if (anyNA(nodes$cl) || anyDuplicated(nodes$cl)) stop('Node IDs must be unique.')
  from <- match(as.character(knn.cl.df$cl.from), nodes$cl)
  to <- match(as.character(knn.cl.df$cl.to), nodes$cl)
  frac <- knn.cl.df$frac
  keep <- !is.na(from) & !is.na(to) & from != to & is.finite(frac) & frac > 0
  from <- from[keep]; to <- to[keep]; frac <- frac[keep]
  key <- paste(from, to, sep = ':')
  rev.idx <- match(paste(to, from, sep = ':'), key)
  # Keep one ribbon for a reciprocal pair, with its two directional widths.
  take <- which(is.na(rev.idx) | seq_along(rev.idx) < rev.idx)
  f <- from[take]; t <- to[take]
  frac.to <- frac[rev.idx[take]]
  frac.to[is.na(frac.to)] <- 0.01  # original narrow taper for one-way edges
  line.segments <- data.frame(
    cl.from = nodes$cl[f], cl.to = nodes$cl[t],
    x.from = nodes$x[f], y.from = nodes$y[f], x.to = nodes$x[t], y.to = nodes$y[t],
    frac.from = frac[take], frac.to = frac.to,
    node.size.from = nodes$node.width[f] / 10,
    node.size.to = nodes$node.width[t] / 10
  )
  # A ribbon between coincident centroids has no visible path.
  line.segments <- line.segments[with(line.segments, x.from != x.to | y.from != y.to), , drop = FALSE]
  poly.Edges <- data.frame(x = numeric(), y = numeric(), Group = character(),
                           g1 = character(), g2 = character())
  if (nrow(line.segments) > 0L) {
    max.frac <- max(line.segments$frac.from, line.segments$frac.to)
    line.segments$ex.line.from <- with(line.segments,
      pmin(frac.from / max.frac * node.size.from * exxageration, node.size.from))
    line.segments$ex.line.to <- with(line.segments,
      pmin(frac.to / max.frac * node.size.to * exxageration, node.size.to))
    polygons <- lapply(seq_len(nrow(line.segments)), function(i) {
      path <- edgeMaker(i, len = 50L, curved = curved, line.segments = line.segments)
      n <- nrow(path)
      xy <- vapply(seq_len(n), function(j) {
        ij <- seq.int(max(1L, j - 1L), min(n, j + 1L))
        half.width <- path$fraction[j] / 2
        if (j == 1L) p <- perpStart(path$x[ij], path$y[ij], half.width)
        else if (j == n) p <- perpEnd(path$x[ij], path$y[ij], half.width)
        else p <- perpMid(path$x[ij], path$y[ij], half.width)
        c(p[1L, ], p[2L, ])
      }, numeric(4L))
      data.frame(x = c(xy[1L, ], rev(xy[3L, ])),
                 y = c(xy[2L, ], rev(xy[4L, ])), Group = as.character(i),
                 g1 = line.segments$cl.from[i], g2 = line.segments$cl.to[i])
    })
    poly.Edges <- do.call(rbind, polygons)
  }
  # Preserve integer endpoint IDs for callers that use the original numeric IDs.
  if (is.numeric(nodes$cluster_id) && all(nodes$cl == as.character(nodes$cluster_id)) &&
      all(nodes$cluster_id == as.integer(nodes$cluster_id))) {
    poly.Edges$g1 <- as.integer(poly.Edges$g1)
    poly.Edges$g2 <- as.integer(poly.Edges$g2)
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
    edge_marking <- match.arg(edge_marking)
    if(edge_marking == "dim"){
      poly.Edges$alpha_plot <- rep(fg.alpha, nrow(poly.Edges)) # Use a different column name to avoid conflict
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- bg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- bg.alpha
      
    } else if(edge_marking == "highlight"){
      poly.Edges$alpha_plot <- rep(bg.alpha, nrow(poly.Edges)) # Use a different column name
      poly.Edges$alpha_plot[poly.Edges$g1 %in% edge_mark_list] <- fg.alpha
      poly.Edges$alpha_plot[poly.Edges$g2 %in% edge_mark_list] <- fg.alpha
      
    } else{
      print("provide valid edge marking")
    }
    
  } else{
    poly.Edges$alpha_plot <- rep(fg.alpha, nrow(poly.Edges)) # Use a different column name
  }
  
  # --- FIX 2: Consolidate plot.all creation and fix alpha/text logic ---
  # Define the base plot structure once
  base_plot <- ggplot() + theme_void() + theme(legend.position = "none")
  if (nrow(poly.Edges) > 0L) {
    base_plot <- base_plot + geom_polygon(data=poly.Edges,
                 aes(x=x, y=y, group=Group, alpha=alpha_plot, # Use new alpha_plot variable
                     # Add hover text for polygons (edges)
                     text = paste0("From: ", g1, "<br>To: ", g2)),
                 fill = "grey60") + # Explicit fill for polygons
    scale_alpha_identity()
  }
  
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
  # The CHARGE app only needs the interactive plot, not PDF legend grobs.
  if (is.null(out.dir) || (enable_plotly && !plot.parts)) {
    if (return.list) return(list(constellation = final_plot_output,
                                 edges.df = poly.Edges, nodes.df = nodes))
    return(invisible(NULL))
  }

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
  
  
  width.1 <- max(c(0, line.segments$frac.from, line.segments$frac.to))
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
  if (!nrow(line.segments)) edge.width.legend <- ggplot() + theme_void()
  if (is.null(dot.size.legend)) dot.size.legend <- grid::nullGrob()
  if (is.null(cl.center.legend)) cl.center.legend <- grid::nullGrob()
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

  t <- seq(0, 1, length.out = len)
  edge <- data.frame(
    x = (1 - t)^2 * fromC[1] + 2 * (1 - t) * t * bezierMid[1] + t^2 * toC[1],
    y = (1 - t)^2 * fromC[2] + 2 * (1 - t) * t * bezierMid[2] + t^2 * toC[2]
  )

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
