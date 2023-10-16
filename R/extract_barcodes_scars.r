#' Extract iTracer barcodes
#' 
#' Extract iTracer barcode from the BAM file of iTracer barcode library to the barcode reference with CellRanger
#' 
#' @import Rsamtools
#' @import dplyr
#' @param file BAM file
#' @param yield_size Chunk size when reading the BAM file. Default is to read in the whole file
#' @param min_reads Minimal read number for a UMI to be considered as valid
#' @param library Label of the library
#' @rdname extract_barcode
#' @export extract_barcodes
extract_barcodes <- function(file,
                             yield_size = NA,
                             min_reads = 3,
                             library = NA)
{
  barcode_pats <- list(pSBbi_1 = '.*AGTGATCC([A-Z]{8,})CCAACC.*',
                       pSBbi_2 = '.*AGCTGATCC([A-Z]{8,})CCAAGC*',
                       pSMAL = '.*AAACCGGT([A-Z]{6,})GAATTCGA.*')
  
  bamFile <- Rsamtools::BamFile(file)
  Rsamtools::yieldSize(bamFile) <- yield_size
  param <- Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(), tag = c('CR','CB','UR','UB'))
  open(bamFile)
  
  df_barcodes <- setNames(data.frame(matrix(nrow=0, ncol=8)),c('Sequence','Barcode','Library','Cell','Reads','UMI','GeneBarcode','Plasmid'))
  while(TRUE){
    cont <- Rsamtools::scanBam(bamFile, param = param)[[1]]
    umis <- cont$tag$UR
    umis[which(!is.na(cont$tag$UB))] <- cont$tag$UB[which(!is.na(cont$tag$UB))]
    cells <- cont$tag$CR
    cells[which(!is.na(cont$tag$CB))] <- cont$tag$CB[which(!is.na(cont$tag$CB))]
    seq <- as.character(cont$seq)
    
    if(length(seq)==0)
      break
    
    idx <- lapply(barcode_pats, function(pat) grep(pat, seq, perl = T))
    if(sum(lengths(idx)>0)==0)
      next
    
    df_block <- Reduce(rbind, lapply(names(barcode_pats)[lengths(idx)>0], function(x){
      pat <- barcode_pats[[x]]
      idx <- idx[[x]]
      bc <- gsub(pat, '\\1', seq[idx], perl=TRUE)
      data.frame(Sequence = seq[idx],
                 Barcode = gsub('\\-[0-9]+$','',cells[idx]),
                 Library = library,
                 Cell = cells[idx],
                 Reads = 1,
                 UMI = umis[idx],
                 GeneBarcode = bc,
                 Plasmid = gsub('_[0-9]+','',x))
    }))
    df_barcodes <- rbind(df_barcodes[,colnames(df_block)], df_block) %>%
      group_by(Sequence,Barcode,Library,Cell,UMI,GeneBarcode,Plasmid) %>%
      summarise(Reads = sum(Reads), .groups='keep')
  }
  df_barcodes <- df_barcodes %>%
    filter(Reads >= min_reads)
  
  return(df_barcodes)
}

#' Extract iTracer scars
#' 
#' Extract iTracer scars from the BAM file of iTracer scar library to the scar reference with CellRanger
#' 
#' @import Rsamtools
#' @import dplyr
#' @param file BAM file
#' @param yield_size Chunk size when reading the BAM file. Default is to read in the whole file
#' @param target Names of the valid target sequences in the reference
#' @param min_reads Minimal read number for a UMI to be considered as valid
#' @param library Label of the library
#' @rdname extract_scars
#' @export extract_scars
extract_scars <- function(file,
                          yield_size = NA,
                          target = c('Scar','Scar_Tomato'),
                          min_reads = 3,
                          library = NA)
{
  illegal_cigar_pat <- c('^[0-9]+S', 'S$', 'N$', '^[0-9]+N')
  
  bamFile <- Rsamtools::BamFile(file)
  Rsamtools::yieldSize(bamFile) <- yield_size
  param <- Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(), tag = c('CR','CB','UR','UB'))
  open(bamFile)
  
  df_scars <- setNames(data.frame(matrix(nrow=0, ncol=9)),c('Sequence','Barcode','Library','Cell','UMI','Reads','CIGAR','Scar','Target'))
  while(TRUE){
    cont <- Rsamtools::scanBam(bamFile, param = param)[[1]]
    umis <- cont$tag$UR
    umis[which(!is.na(cont$tag$UB))] <- cont$tag$UB[which(!is.na(cont$tag$UB))]
    cells <- cont$tag$CR
    cells[which(!is.na(cont$tag$CB))] <- cont$tag$CB[which(!is.na(cont$tag$CB))]
    seq <- as.character(cont$seq)
    cigar <- cont$cigar
    start <- cont$pos
    targeted <- cont$rname
    
    if(length(seq)==0)
      break
    
    idx <- which(targeted %in% target)
    idx_legal <- which(rowSums(sapply(illegal_cigar_pat, function(pat) grepl(pat, cigar)))==0)
    idx <- intersect(idx, idx_legal)
    df_block <- data.frame(Sequence = seq[idx],
                           Barcode = gsub('\\-[0-9]+$','',cells[idx]),
                           Library = library,
                           Cell = cells[idx],
                           UMI = umis[idx],
                           Reads = 1,
                           CIGAR = cigar[idx],
                           Scar = paste(start[idx], cigar[idx], sep=':'),
                           Target = targeted[idx])
    
    df_scars <- rbind(df_scars[,colnames(df_block)], df_block) %>%
      group_by(Sequence,Barcode,Library,Cell,UMI,CIGAR,Scar,Target) %>%
      summarise(Reads = sum(Reads), .groups='keep')
  }
  df_scars <- df_scars %>%
    filter(Reads >= min_reads) %>%
    mutate(Target = factor(Target, levels=target))
  
  return(df_scars)
}
