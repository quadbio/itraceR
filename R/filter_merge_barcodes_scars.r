#' Filter the iTracer barcode data frame
#' 
#' Filter the iTracer barcode data frame output by the extraction, based on number of supporting
#' reads per UMI, barcode sequences, UMI conflicts between barcodes, and UMI-barcode conflicts
#' between cells.
#' 
#' 
#' @import stringdist
#' @import igraph
#' @import reshape2
#' @import mixtools
#' @import dplyr
#' @import ggplot2
#' @param df_barcode Barcode data frame
#' @param coverage_cutoff_mode Algorithm to determine read coverage threshold. 'linear' is the default method described in the paper. 'log_gmm' applies a bimodal Gaussian mixture model. 'log_turn' is similar to 'linear' but with read number log-transformed.
#' @param cutoff_prob Cutoff of the bimodal Gaussian mixture model posterior probability. Only applicable when coverage_cutoff_mode is 'log_gmm'.
#' @param num_log_breaks Number of breaks to count frequency on the log-transformed read coverage per UMI. Only applicable when coverage_cutoff_mode is 'log_turn'.
#' @param do_plot Plot read coverage histogram and the threshold
#' @rdname do_filter_barcodes
#' @export do_filter_barcodes
do_filter_barcodes <- function(df_barcode, # table output by extract_barcodes()
                               coverage_cutoff_mode = c("linear","log_gmm","log_turn"),
                               cutoff_prob = 1 - 1E-10,
                               num_log_breaks = 20,
                               do_plot = FALSE) 
{
  coverage_cutoff_mode <- match.arg(coverage_cutoff_mode)
  
  # filter 1: read coverage
  if (coverage_cutoff_mode == "linear"){
    count <- melt(table(df_barcode$Reads))
    count <- merge(data.frame(count=3:max(count$Var1)),count, by.y = "Var1",by.x = "count",all = TRUE)
    count[is.na(count)] <- 0
    count$fit <- predict(loess(value~log(count), count, span = 0.1), count)
    count$fit[count$fit<1] <- 1
    turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
    idx_high_coverage <- which(df_barcode$Reads >= turn)
  } else if (coverage_cutoff_mode == "log_gmm"){
    logReads <- log(df_barcode$Reads)
    model <- normalmixEM(logReads, lambda = c(0.5,0.5), mu = c(min(logReads), max(logReads)))
    idx_high_coverage <- which(model$posterior[,which.max(model$mu)] >= cutoff_prob)
  } else{
    logReads <- log(df_barcode$Reads)
    breaks <- as.numeric(cut(logReads, breaks = num_log_breaks, right = F, include.lowest = T))
    count <- melt(table(breaks))
    count <- merge(data.frame(count=1:max(count$breaks)),count, by.y = "breaks",by.x = "count",all = TRUE)
    count[is.na(count)] <- 0
    count$fit <- predict(loess(value~count, count, span = 1), count)
    count$fit[count$fit<1] <- 1
    turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
    idx_high_coverage <- which(breaks >= turn)
  }
  if (do_plot){
    plt <- ggplot(df_barcode, aes(x= Reads)) +
      geom_histogram(bins=100) +
      scale_y_log10() + 
      xlab("Barcode reads (n.)") +
      ylab("Occurrance (n.)")+
      geom_vline(xintercept = min(df_barcode$Reads[idx_high_coverage]), linetype = "dashed", color = "red")+
      theme_minimal()
    print(plt)
  }
  df_barcode <- df_barcode[idx_high_coverage,]
  
  # filter 2: start and end with either A or T
  idx_reasonable_seq <- which(substring(df_barcode$GeneBarcode, 1, 1) %in% c("A","T") & substring(df_barcode$GeneBarcode, nchar(df_barcode$GeneBarcode), nchar(df_barcode$GeneBarcode)) %in% c("A","T"))
  df_barcode <- df_barcode[idx_reasonable_seq,]
  
  # filter 3: if there are molecules where the same UMI for different barcodes is detected in a cell, only remain the one with higher coverage
  df_barcode_cells <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), df_barcode$Cell)
  idx_max_coverage_per_umi <- sort(unlist(lapply(df_barcode_cells, function(df){
    df_each_umi <- split(df, df$UMI)
    return(unlist(lapply(df_each_umi, function(df) df$idx[which.max(df$Reads)])))
  })))
  df_barcode <- df_barcode[idx_max_coverage_per_umi,]
  
  # filter 4: filter out molecules where the same UMI and barcode is detected across multiple cells
  df_barcode_barcode_umi <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), paste0(df_barcode$UMI, "_", df_barcode$GeneBarcode))
  df_barcode_barcode_umi <- df_barcode_barcode_umi[sapply(df_barcode_barcode_umi, nrow) == 1]
  idx_once_barcode_umi <- sort(unlist(lapply(df_barcode_barcode_umi, function(x) x$idx)))
  df_barcode <- df_barcode[idx_once_barcode_umi,]
  
  # filter 5: filter out molecules where the hamming distance between detected barcodes within a cell is 1, to pick the barcode that are covered by more UMI or that have a higher coverage
  df_barcode_cells <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), df_barcode$Cell)
  idx_diff_enough <- sort(unlist(lapply(df_barcode_cells, function(df){
    barcodes <- sort(unique(df$GeneBarcode))
    umi_barcodes <- table(df$GeneBarcode)[barcodes]
    read_barcodes <- sapply(split(df[,c("Reads","GeneBarcode")], df$GeneBarcode), function(x) sum(x$Reads))[barcodes]
    
    hamming_dist <- stringdist::stringdistmatrix(barcodes, barcodes, method="hamming")
    hamming_adj <- hamming_dist <= 1
    rownames(hamming_adj) <- barcodes
    colnames(hamming_adj) <- barcodes
    graph <- igraph::graph_from_adjacency_matrix(hamming_adj, mode = "undirected")
    components <- setNames(igraph::components(graph)$membership, names(igraph::V(graph)))
    valid_barcodes <- sapply(tapply(barcodes, components[barcodes], list), function(x) x[order(umi_barcodes[x], read_barcodes[x], decreasing=T)[1]] )
    
    return(df$idx[df$GeneBarcode %in% valid_barcodes])
  })))
  df_barcode <- df_barcode[idx_diff_enough,]
  
  return(df_barcode)
}

#' Filter the iTracer scar data frame
#' 
#' Filter the iTracer scar data frame output by the extraction, based on number of supporting
#' reads per UMI, bUMI conflicts between barcodes, and UMI-barcode conflicts between cells.
#' 
#' 
#' @import stringdist
#' @import igraph
#' @import reshape2
#' @import mixtools
#' @import dplyr
#' @import ggplot2
#' @param df_scar Scar data frame
#' @param coverage_cutoff_mode Algorithm to determine read coverage threshold. 'linear' is the default method described in the paper. 'log_gmm' applies a bimodal Gaussian mixture model. 'log_turn' is similar to 'linear' but with read number log-transformed.
#' @param cutoff_prob Cutoff of the bimodal Gaussian mixture model posterior probability. Only applicable when coverage_cutoff_mode is 'log_gmm'.
#' @param num_log_breaks Number of breaks to count frequency on the log-transformed read coverage per UMI. Only applicable when coverage_cutoff_mode is 'log_turn'.
#' @param do_plot Plot read coverage histogram and the threshold
#' @rdname do_filter_scars
#' @export do_filter_scars
do_filter_scars <- function(df_scar,
                            coverage_cutoff_mode = c("linear","log_gmm","log_turn"),
                            cutoff_prob = 1 - 1E-10,
                            num_log_breaks = 20,
                            do_plot = FALSE) # table output by ExtractScarsFromBam.pl
{
  coverage_cutoff_mode <- match.arg(coverage_cutoff_mode)
  
  # filter 1: read coverage
  if (max(df_scar$Reads) <= 5){
    idx_high_coverage <- 1:nrow(df_scar)
  } else if (coverage_cutoff_mode == "linear"){
    count <- melt(table(df_scar$Reads))
    count <- merge(data.frame(count=3:max(count$Var1)),count, by.y = "Var1",by.x = "count",all = TRUE)
    count[is.na(count)] <- 0
    count$fit <- predict(loess(value~log(count), count, span = 0.1), count)
    count$fit[count$fit<1] <- 1
    turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
    idx_high_coverage <- which(df_scar$Reads >= turn)
  } else if (coverage_cutoff_mode == "log_gmm"){
    logReads <- log(df_scar$Reads)
    model <- normalmixEM(logReads, lambda = c(0.5,0.5), mu = c(min(logReads), max(logReads)))
    idx_high_coverage <- which(model$posterior[,which.max(model$mu)] >= cutoff_prob)
  } else{
    logReads <- log(df_scar$Reads)
    breaks <- as.numeric(cut(logReads, breaks = num_log_breaks, right = F, include.lowest = T))
    count <- melt(table(breaks))
    count <- merge(data.frame(count=1:max(count$breaks)),count, by.y = "breaks",by.x = "count",all = TRUE)
    count[is.na(count)] <- 0
    count$fit <- predict(loess(value~count, count, span = 1), count)
    count$fit[count$fit<1] <- 1
    turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
    idx_high_coverage <- which(breaks >= turn)
  }
  if (do_plot){
    plt <- ggplot(df_scar, aes(x= Reads)) +
      geom_histogram(bins=100) +
      scale_y_log10() + 
      xlab("Barcode reads (n.)") +
      ylab("Occurrance (n.)")+
      geom_vline(xintercept = min(df_scar$Reads[idx_high_coverage]), linetype = "dashed", color = "red")+
      theme_minimal()
    print(plt)
  }
  df_scar <- df_scar[idx_high_coverage,]
  
  # filter 2: if there are molecules where the same UMI for different scars is detected in a cell, only remain the one with higher coverage
  df_scar_cells <- split(data.frame(idx = 1:nrow(df_scar), df_scar, stringsAsFactors=F), df_scar$Cell)
  idx_max_coverage_per_umi <- sort(unlist(lapply(df_scar_cells, function(df){
    df_each_umi <- split(df, df$UMI)
    return(unlist(lapply(df_each_umi, function(df) df$idx[which.max(df$Reads)])))
  })))
  df_scar <- df_scar[idx_max_coverage_per_umi,]
  
  # filter 3: filter out molecules where the same UMI and scar is detected across multiple cells
  df_scar_scar_umi <- split(data.frame(idx = 1:nrow(df_scar), df_scar, stringsAsFactors=F), paste0(df_scar$UMI, "_", df_scar$Scar))
  df_scar_scar_umi <- df_scar_scar_umi[sapply(df_scar_scar_umi, nrow) == 1]
  idx_once_barcode_umi <- sort(unlist(lapply(df_scar_scar_umi, function(x) x$idx)))
  df_scar <- df_scar[idx_once_barcode_umi,]
  
  return(df_scar)
}


#' Merge the filtered iTracer barcode and scar data frames
#' 
#' Merge the filtered iTracer barcode and scar data frames, by considering only those presented
#' at the same UMIs in the same cells
#' 
#' @param df_barcode Barcode data frame
#' @param df_scar Scar data frame
#' @rdname merge_barcode_scars
#' @export merge_barcode_scars
merge_barcode_scars <- function(df_barcode, # filtered barcode table
                                df_scar) # filtered scar table
{
  # only shared cells are considered
  cells <- intersect(df_barcode$Cell, df_scar$Cell)
  df_barcode <- df_barcode[df_barcode$Cell %in% cells,]
  df_scar <- df_scar[df_scar$Cell %in% cells,]
  
  # only shared UMIs in the same cells are considered
  df_barcode_cells <- split(df_barcode, df_barcode$Cell)
  df_scar_cells <- split(df_scar, df_scar$Cell)
  df_barcode_scar <- do.call(rbind, lapply(cells, function(cell){
    df_barcode_cell <- df_barcode_cells[[cell]]
    df_scar_cell <- df_scar_cells[[cell]]
    umis <- intersect(df_barcode_cell$UMI, df_scar_cell$UMI)
    
    barcode_family <- setNames(df_barcode_cell$GeneBarcode, df_barcode_cell$UMI)[umis]
    scar_family <- setNames(df_scar_cell$Scar, df_scar_cell$UMI)[umis]
    orders <- order(barcode_family)
    df <- unique(data.frame(barcode = barcode_family[orders], scar = scar_family[orders], stringsAsFactors = F))
    return(data.frame(cell = cell,
                      GeneBarcodeMerged = paste(df$barcode, collapse="_"),
                      ScarMerged = paste(df$scar, collapse="_"), stringsAsFactors=F))
  }))
  df_barcode_scar <- df_barcode_scar[which(df_barcode_scar$GeneBarcodeMerged!=""),]
  return(df_barcode_scar)
}
