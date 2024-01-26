#!/usr/bin/env Rscript

# Load necessary libraries and install them if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("ggplot2", "cowplot")
bioconductor_packages <- c("rGADEM", "motifStack", "Biostrings", "DiffBind", "ChIPseeker", 
                           "TxDb.Hsapiens.UCSC.hg19.knownGene", "ComplexHeatmap")

# Install CRAN packages if they are not installed
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install Bioconductor packages if they are not installed
for (pkg in bioconductor_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load necessary libraries
library(rGADEM)
library(motifStack)
library(Biostrings)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ComplexHeatmap)

# Quality Control of Raw Reads
run_fastqc <- function(fastq_files, output_dir) {
  fastqc_cmd <- paste("fastqc", paste(fastq_files, collapse=" "), "-o", output_dir)
  system(fastqc_cmd)
}

# Read Alignment with Bowtie2
run_bowtie2 <- function(index_base, read_files, output_sam) {
  bowtie2_cmd <- paste("bowtie2 -x", index_base, "-1", read_files[1], "-2", read_files[2], "-S", output_sam)
  system(bowtie2_cmd)
}

# Peak Calling with MACS2
run_macs2 <- function(treatment_bam, control_bam=NULL, output_dir, genome_size="hs", is_pe=FALSE) {
  macs2_cmd <- paste("macs2 callpeak -t", treatment_bam, "-f", ifelse(is_pe, "BAMPE", "BAM"), "-g", genome_size, "--outdir", output_dir)
  if (!is.null(control_bam)) {
    macs2_cmd <- paste(macs2_cmd, "-c", control_bam)
  }
  system(macs2_cmd)
}

# Annotation of Peaks
annotate_peaks <- function(peak_file) {
  peaks <- readPeakFile(peak_file)
  annotated_peaks <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
  return(annotated_peaks)
}

# Motif Analysis with GADEM
run_gadem <- function(peak_sequences) {
  gadem <- GADEM(peak_sequences)
  return(gadem)
}

# Differential Binding Analysis with DiffBind
run_diffbind <- function(sampleSheet) {
  dba <- dba(sampleSheet = sampleSheet)
  dba <- dba.count(dba)
  dba <- dba.contrast(dba, categories=DBA_CONDITION)
  dba <- dba.analyze(dba)
  results <- dba.report(dba, th=1, method=DBA_ALL_METHODS)
  return(results)
}

# MA Plot for Differential Binding Analysis
plot_ma <- function(results) {
  # Check if 'AveExpr' column exists; if not, calculate it or use an alternative
  if (!"AveExpr" %in% colnames(results)) {
    # Placeholder: AveExpr needs to be calculated or set to a default column
    results$AveExpr <- rowMeans(results[, c("ColumnName1", "ColumnName2")])  # Adjust with actual column names
  }
  
  ggplot_object <- ggplot(results, aes(x = AveExpr, y = log2FoldChange)) +
    geom_point(aes(color = FDR < 0.05)) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = "MA Plot", x = "Average Expression", y = "Log2 Fold Change") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  print(ggplot_object)  # Print the plot
}

# Volcano Plot for Differential Binding Analysis
plot_volcano <- function(results) {
  ggplot_object <- ggplot(results, aes(x = log2FoldChange, y = -log10(FDR))) +
    geom_point(aes(color = FDR < 0.05)) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 FDR") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    theme(legend.position = "none")  # Hide the legend
  print(ggplot_object)  # Print the plot
}

# Heatmap for Differential Binding Analysis
plot_heatmap <- function(normalized_data, results) {
  # Assuming 'results' is a DataFrame with a 'gene' column and 'normalized_data' is a matrix with rownames as genes
  # Select top differentially expressed genes for the heatmap
  top_genes <- head(order(results$FDR), 50)  # Get the top 50 genes by FDR
  data_subset <- normalized_data[top_genes, ]
  
  # Create the heatmap
  heatmap_object <- Heatmap(data_subset, name = "Binding Intensity", show_row_names = TRUE, show_column_names = TRUE)
  draw(heatmap_object)  # Draw the heatmap
}

# Main function to run the workflow
run_chipseq_workflow <- function() {
  # Example usage:
  fastq_files <- c("sample1.fastq", "sample2.fastq")  # List your FASTQ files here
  output_dir <- "fastqc_results"
  run_fastqc(fastq_files, output_dir)
  
  index_base <- "path/to/genome_index"
  read_files <- c("sample1_R1.fastq", "sample1_R2.fastq")  # Replace with your read files
  output_sam <- "aligned_reads.sam"
  run_bowtie2(index_base, read_files, output_sam)
  
  treatment_bam <- "path/to/treatment.bam"
  control_bam <- "path/to/control.bam"
  output_dir <- "macs2_output"
  run_macs2(treatment_bam, control_bam, output_dir)
  
  peak_file <- "macs2_output/peaks.narrowPeak"
  annotated_peaks <- annotate_peaks(peak_file)
  
  # Assuming the peak sequences are extracted and stored in 'peak_sequences'
  peak_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, names(annotated_peaks))
  gadem <- run_gadem(peak_sequences)
  
  sampleSheet <- "samples.csv"
  results <- run_diffbind(sampleSheet)
  
  # Assuming 'results' and 'normalized_data' are already defined
  plot_ma(results)
  plot_volcano(results)
  plot_heatmap(normalized_data, results)
}

# Run the workflow
run_chipseq_workflow()
