# itraceR
The R package to preprocess barcode and scar data of iTracer. iTracer is the dual-channel lineage recorder that combines reporter barcodes with inducible CRISPR–Cas9 scarring and is compatible with single-cell and spatial transcriptomics. The iTracer system has been successfully applied to brain organoids and introduced in the paper '[Lineage recording in human cerebral organoids](https://www.nature.com/articles/s41592-021-01344-8)' published in Nature Methods in 2021. The originally released codes and notebooks can be found in the [iTracer_analysis](https://github.com/quadbio/iTracer_analysis) repository.

## Quick start
The package can be installed with the *devtools* package in R
```R
devtools::install_github('quadbio/itraceR')
```

Once it is done, the *itraceR* package can then be loaded and used for the iTracer readout preprocessing
```R
library(itraceR)

# Extract barcode UMIs
# Note: file_barcode is the path to the BAM file output by CellRanger mapping of the barcode library (e.g. barcode/outs/possorted_genome_bam.bam)
df_barcode <- extract_barcodes(file_barcode)

# Extract scar UMIs
# Note: file_scar is the path to the BAM file output by CellRanger mapping of the scar library (e.g. scar/outs/possorted_genome_bam.bam)
df_scar <- extract_scars(file_scar)

# Filter barcode UMIs based on number of reads, barcode sequences, etc.
# Note: do_plot=T to plot the read number per UMI distribution; remove it if such plot is not needed
df_barcode <- do_filter_barcodes(df_barcode, do_plot=T)

# Filter scar UMIs based on number of reads, etc.
# Note: do_plot=T to plot the read number per UMI distribution; remove it if such plot is not needed
df_scar <- do_filter_scars(df_scar , do_plot=T)

# Merge barcode UMIs and scar UMIs data frames into the final iTracer readout data frame, by taking only barcode UMIs and scar UMIs sharing the same cell barcode and UMI barcode
df_barcode_scar <- merge_barcode_scars(df_barcode, df_scar)
```

Next, the resulted iTracer readout data frame can be incorporated into the corresponding Seurat object of the transcriptomic data
```R
# Incorporate the iTracer readout data frame obtained above to the metadata of the corresponding Seurat object
# Note: in this example, seurat is the Seurat object for the corresponding transcriptomic data
# Note: by default the barcode information is added to the new column ‘GeneBarcodeMerged’, and the scar is added to ‘ScarMerged’
seurat <- incorporate_to_seurat(seurat, df_barcode_scar)
```

Assuming that the transcriptomic data has been analyzed with cell type label stored as 'annot' in the object, the wheel lineage tree similar to those presented in the paper can be made as following
```R
# Plot the lineage tree
# Note: if col_lab is not specified, the function assumes no cell type annotation, and all cells will be colored the same
plot_tree(seurat@meta.data, col_lab = 'annot')
```
