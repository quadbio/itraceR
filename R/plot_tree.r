#' Plot the lineage tree
#' 
#' Plot the lineage tree in the wheel shape
#' 
#' @import igraph
#' @import dplyr
#' @import ggplot2
#' @import ggraph
#' @param metadata_df Metadata data frame
#' @param col_lab The metadata column name of cell annotations
#' @param col_barcode The metadata column name of iTracer barcodes
#' @param col_scar The metadata column name of iTracer scars
#' @param edge_alpha Edge alpha of the tree
#' @param edge_width Edge width of the tree
#' @rdname plot_tree
#' @export plot_tree
plot_tree <- function(metadata_df,
                      col_lab = NA,
                      col_barcode = 'GeneBarcodeMerged',
                      col_scar = 'ScarMerged',
                      edge_alpha = 0.25,
                      edge_width = 0.1)
{
  if(is.na(col_lab)){
    metadata_df[,'_mockLab'] <- 'mock'
    col_lab <- '_mockLab'
  }
  tree_df1 <- data.frame("cells" = rownames(metadata_df), 
                         "barcode"= metadata_df[,col_barcode], 
                         "scar" = metadata_df[,col_scar])
  
  tree_df1 <- filter(tree_df1, 
                     barcode %in% unique(tree_df1$barcode)[!is.na(unique(tree_df1$barcode))])
  
  tree_df2 <- data.frame(root = "root", 
                         barcode = tree_df1$barcode, 
                         scar = tree_df1$scar, 
                         cells = tree_df1$cells) %>%
    mutate(barcode_scar = paste(barcode, scar, sep='_')) %>%
    arrange(barcode)
  
  tree_df2_root_bar <- data.frame(from = tree_df2$root, 
                                  to = tree_df2$barcode, 
                                  label = rep("barcode", length(tree_df2$root)))
  tree_df2_bar_scar <- data.frame(from = tree_df2$barcode, 
                                  to = tree_df2$barcode_scar, 
                                  label = rep("scar", length(tree_df2$barcode)))
  tree_df2_scar_cell <- data.frame(from = tree_df2$barcode_scar,
                                   to = tree_df2$cells,
                                   label = as.character(metadata_df[tree_df2$cells, col_lab]))
  
  tree_df3 <- rbind(tree_df2_root_bar,
                    tree_df2_bar_scar,
                    tree_df2_scar_cell)
  
  root_label <- data.frame(to = "root", label = "root")
  
  node_labels <- unique(tree_df3[,2:3])
  node_labels <- rbind(node_labels, root_label)
  rownames(node_labels) = node_labels[,1]
  
  g <- graph_from_data_frame(tree_df3)
  V(g)$label <- as.character(node_labels[V(g)$name,"label"])
  tree <- ggraph(g, 'dendrogram', circular = TRUE) + 
    geom_edge_diagonal(alpha = edge_alpha, width = edge_width) + 
    geom_node_point(aes(colour = label)) +
    coord_fixed() +
    theme_void() +
    theme(aspect.ratio = 1)
  
  return(tree)
}
