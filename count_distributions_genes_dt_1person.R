

# -------------------------------------------------------------------------
# TPM Distributions  ------------------------------------------------------
# -------------------------------------------------------------------------

  setwd('/Users/meco9877/Documents/Projects/ifnb_regression/')
  library(tidyverse)
  savedir_figs <- 'proseq/figures/counts/'

  # Load in the list of genes DE in only 1 individual
  gene_in1 <- read_tsv('../proseq/data/deseq_res/intrahuman_ifnb_proseq_genes_dein1.tsv')
  colnames(gene_in1) <- c('Feature', 'ID_DE')

  # Load df with counts related to transcript isoform
  counts_raw <- read_csv('../proseq/data/counts/pop_proseq_cleaned_counts.csv')
  
  # Load information on which isoform is longest per gene
  longest_isoform <- read_tsv("~/Documents/Reference/resources/annotations/hg38_diff53prime_refseq_longisof_noputatives_transcripts.bed",
                              col_names = FALSE)
  colnames(longest_isoform) <- c('Chr', 'Start', 'End', 'Isoforms', 'Blank', 'Strand', 'Longest')
  longest_isoform$isoform <- sub('.*:(.*)\\..*$', '\\1', longest_isoform$Longest)
  
  # Get the list of longest isoforms and associated length
  counts_longest <- counts_raw %>%
    filter(Transcript %in% longest_isoform$isoform) 
  
  counts_longest[1:3,]

  

# Calculate TPMs ----------------------------------------------------------

  counts_tpm <- counts_longest[, -c(2:7)]
  
  # 1. Divide by length of each gene in kilobases
  counts_tpm[, 3:dim(counts_tpm)[2]] <- counts_tpm[, 3:dim(counts_tpm)[2]] / (counts_tpm[[2]]/1000)
  
  # 2. Count up RPK values per sample and divide by 10^6
  rpk_per_sample <- colSums(counts_tpm[, 3:dim(counts_tpm)[2]]) / 10^6
  
  # 3. Divide RPK by RPK per sample
  counts_tpm[, 3:dim(counts_tpm)[2]] <- counts_tpm[, 3:dim(counts_tpm)[2]] / rpk_per_sample
  write_csv(counts_tpm, '../proseq/data/counts/counts_tpm.csv')
  
  

# Filter by Avg TPM = 1 ---------------------------------------------------
  
  counts_tpm$avg_tpm <- rowMeans(counts_tpm[, -c(1:2)])
  tpm_expressed <- counts_tpm[counts_tpm$avg_tpm > 1, ]
  

# Split into two pools, DE in 1 and other ---------------------------------

  de_in1_genes <- gene_in1$Feature
  tpm_de_in1 <- tpm_expressed[tpm_expressed$Gene %in% de_in1_genes, ]
  tpm_other <- tpm_expressed[!tpm_expressed$Gene %in% de_in1_genes, ]
  dim(tpm_de_in1)
  dim(tpm_other)
  dim(tpm_expressed)
  

# Get quintiles of lengths  -----------------------------------------------

  set.seed(42)
  n_bins = 10
  
  bin_breaks <- quantile(tpm_de_in1$Length, probs = seq(0, 1, length.out = n_bins + 1))
  bin_breaks[1] <- -Inf
  bin_breaks[length(bin_breaks)] <- Inf
  
  tpm_de_in1$length_bin <- cut(tpm_de_in1$Length, breaks = bin_breaks, 
                               include.lowest = TRUE)
  tpm_other$length_bin <- cut(tpm_other$Length, breaks = bin_breaks, 
                              include.lowest = TRUE)
  
# Get # genes per bin to sample  ------------------------------------------

  bin_targets <- table(tpm_de_in1$length_bin)
  bin_targets
  

# Sample from other  ------------------------------------------------------

  matched_list <- vector("list", length(bin_targets))
  
  for (i in seq_along(bin_targets)){
    bin_name <- names(bin_targets)[i]
    n_needed <- bin_targets[i]
    pool <- tpm_other[tpm_other$length_bin == bin_name, ]
    
    matched_list[[i]] <- pool[sample(nrow(pool), n_needed, replace = FALSE), ]
  }
  
  tpm_other_matched <- do.call(rbind, matched_list)
  

# Double check selection  -------------------------------------------------

  plot(density(tpm_de_in1$Length), col = "blue", lwd = 2)
  lines(density(tpm_other_matched$Length), col = "red", lwd = 1)
  
  
# Combine df and plot  ----------------------------------------------------

  tpm_de_in1$cat <- 'DE in 1 Individual'
  tpm_other_matched$cat <- 'Other'
  
  tpm_plot <- rbind(tpm_de_in1, tpm_other_matched) 
  tpm_plot <- tpm_plot %>% 
    pivot_longer(cols = -c('Gene', 'Length', 'avg_tpm', 'length_bin', 'cat'), 
                 names_to = 'sample', 
                 values_to = 'TPM')
  
  
  ggplot(tpm_plot, aes(x = log10(TPM), 
                       fill = cat, 
                       color = cat)) + 
    geom_density(alpha = 0.4, linewidth = 1) + 
    geom_vline(xintercept = 0)+ 
    theme_light() + 
    theme(
      plot.title = element_text(size = 26, hjust = 0.5), 
      axis.text = element_text(size = 20), 
      axis.title = element_text(size = 23), 
      legend.text = element_text(size = 20), 
      legend.title = element_blank()
    ) +
    scale_fill_manual(values = c('#4A2040', '#F9E0D9')) + 
    scale_color_manual(values = c('#4A2040', '#F9E0D9'))+ 
    xlab(expression(log[10](TPM))) + 
    ggtitle('TPM Distribution for Genes') 

  ggsave(paste0(savedir_figs,
                'tpm_distributions.png'),
         height = 5, width = 10)
  