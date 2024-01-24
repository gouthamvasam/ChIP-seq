#!/usr/bin/env Rscript

# Import necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rGADEM")
BiocManager::install("motifStack")
BiocManager::install("Biostrings")
BiocManager::install("DiffBind")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("ggplot2")
BiocManager::install("ComplexHeatmap")

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

# Read Alignment
run_bowtie2 <- function(index_base, read_files, output_sam) {
  bowtie2_cmd <- paste("bowtie2 -x", index_base, "-1", read_files[1], "-2", read_files[2], "-S", output_sam)
  system(bowtie2_cmd)
}

# Peak Calling
run_macs2 <- function(treatment_bam, control_bam=NULL, output_dir, genome_size="hs", is_pe=FALSE) {
  macs2_cmd <- paste("macs2 callpeak -t", treatment_bam, "-f", ifelse(is_pe, "BAMPE", "BAM"), "-g", genome_size, "--outdir", output_dir)
  if (!is.null(control_bam)) {
    macs2_cmd <- paste(macs2_cmd, "-c", control_bam)
  }
  system(macs2_cmd)
}

# Annotation
annotate_peaks <- function(peak_file) {
  peaks <- readPeakFile(peak_file)
  annotated_peaks <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
  return(annotated_peaks)
}

# Motif Analysis
run_gadem <- function(peak_sequences) {
  gadem <- GADEM(peak_sequences)
  return(gadem)
}

# Differential Binding Analysis
run_diffbind <- function(sampleSheet) {
  dba <- dba(sampleSheet = sampleSheet)
  dba <- dba.count(dba)
  dba <- dba.contrast(dba, categories=DBA_CONDITION)
  dba <- dba.analyze(dba)
  results <- dba.report(dba, th=1, method=DBA_ALL_METHODS)
  return(results)
}

# MA Plot
plot_ma <- function(results) {
  ggplot(results, aes(x = AveExpr, y = log2FoldChange)) +
    geom_point(aes(color = pvalue < 0.05)) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = "MA Plot", x = "Average Expression", y = "Log2 Fold Change") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
}

# Volcano Plot
plot_volcano <- function(results) {
  ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = pvalue < 0.05)) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")
}

# Heatmap
plot_heatmap <- function(normalized_data, de_results) {
  top_genes <- de_results[order(de_results$pvalue), ][1:50, ]$gene
  data_subset <- normalized_data[top_genes, ]
  heatmap <- Heatmap(data_subset, name = "Binding Intensity", show_row_names = FALSE, show_column_names = TRUE)
  draw(heatmap, heatmap_legend_side = "bot")
}

# Main function to run the workflow
run_chipseq_workflow <- function() {
  # Example usage:
  fastq_files <- c("sample1.fastq", "sample2.fastq")
  output_dir <- "fastqc_results"
  run_fastqc(fastq_files, output_dir)
  
  index_base <- "path/to/genome_index"
  read_files <- c("sample1_R1.fastq", "sample1_R2.fastq")
  output_sam <- "aligned_reads.sam"
  run_bowtie2(index_base, read_files, output_sam)
  
  treatment_bam <- "path/to/treatment.bam"
  control_bam <- "path/to/control.bam"
  output_dir <- "macs2_output"
  run_macs2(treatment_bam, control_bam, output_dir)
  
  peak_file <- "macs2_output/peaks.narrowPeak"
  annotated_peaks <- annotate_peaks(peak_file)
  
  peak_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, names(annotated_peaks))
  gadem <- run_gadem(peak_sequences)
  
  sampleSheet <- "samples.csv"
  results <- run_diffbind(sampleSheet)
  
  plot_ma(results)
  plot_volcano(results)
  
  normalized_data <- read.csv("normalized_expression.csv", row.names = 1)
  plot_heatmap(normalized_data, results)
}

# Run the workflow
run_chipseq_workflow()