

# -------------------------------------------------------------------------
# PURPOSE -----------------------------------------------------------------
# -------------------------------------------------------------------------

  # The purpose of this script is to graph results from DESeq2.
  # The current experiment performed PRO-seq and RNA-seq on LCLs from 
  # 9 individuals from Asia, Africa, North America, and South America. 
  # The experimental conditions are 1 hr of interferon-beta treatment 
  # (100ng/uL) and a BSA control. 


# -------------------------------------------------------------------------
# User variables ----------------------------------------------------------
# -------------------------------------------------------------------------

  seq_mode <- 'PRO' # options = c('PRO', 'RNA')
  

# -------------------------------------------------------------------------
# Set environment ---------------------------------------------------------
# -------------------------------------------------------------------------

  
  library(tidyverse)
  library(ComplexHeatmap)
  
  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  path
  source(paste0('deseq2_helper_functions.R'))
  setwd(paste0(path, '/..'))
  getwd()


# -------------------------------------------------------------------------
# Establish directories for loading and saving data -----------------------
# -------------------------------------------------------------------------

  
  if (seq_mode == 'PRO') {
    savedir_figs <- 'proseq/figures/deseq_figs/'
    savedir_data <- 'proseq/data/deseq_res/'
    data_load <- 'proseq/data/counts/pop_pro_dt_res.csv'
    
  } else {
    savedir_figs <- 'rnaseq/figures/deseq_figs/'
    savedir_data <- 'rnaseq/data/deseq_res/'
    data_load <- 'rnaseq/data/deseq_res/pop_rna_dt_res.csv'
  }
  
  print(savedir_figs)
  print(savedir_data)
  print(data_load)

  
# -------------------------------------------------------------------------
# Read in data ------------------------------------------------------------
# -------------------------------------------------------------------------

  deseq_res <- read_csv(data_load)
  

# -------------------------------------------------------------------------
# Upset Plots -------------------------------------------------------------
# -------------------------------------------------------------------------

  # I want an upset plot looking at what is common between people. I will 
  # split upregulated and downregulated genes and look at them separately. 
  
  
# `create_upset_list` takes 
# c(deseq results, direction of change, yes/no include population category)
  upset_mat_up_y <- create_upset_list(deseq_res, 'up', 'yes')
  upset_mat_up_n <- create_upset_list(deseq_res, 'up', 'no')
  
  upset_mat_down_y <- create_upset_list(deseq_res, 'down', 'yes')
  upset_mat_down_n <- create_upset_list(deseq_res, 'down', 'no')
  
#-- Optional step - Eliminate rows related to only a single sample 
  # upset_filtered <- upset_mat[rowSums(upset_mat)>1, ]
  
# Create upset plots that include the population category, save as png
  # upregulated 
  png(filename = paste0(savedir_figs, 'intrahuman_upset_pop_up.png'),
      width = 300, height = 300, units='px')
  make_upset_plot(upset_mat_up_y, 12, 'Genes Upregulated, Including Population')
  dev.off()
  
  # downregulated
  png(filename = paste0(savedir_figs, 'intrahuman_upset_pop_down.png'),
      width = 300, height = 300, units='px')
  make_upset_plot(upset_mat_down_y, 12, 
                  'Genes Downregulated, Including Population')
  dev.off()
  
# Create upset plots that will compare across all individuals, without the
# population category to focus on the inter-individual comparisons
  
  # Upregulated 
  png(filename = paste0(savedir_figs, 'intrahuman_upset_up.png'),
      width = 350, height = 250, units='px')
  make_upset_plot(upset_mat_up_n, 12, 
                  'Genes Upregulated')
  dev.off()
  
  # Downregulated 
  png(filename = paste0(savedir_figs, 'intrahuman_upset_down.png'),
      width = 300, height = 300, units='px')
  make_upset_plot(upset_mat_down_n, 12, 'Genes Downregulated')
  dev.off()
  

# -------------------------------------------------------------------------
# Make barplot ------------------------------------------------------------
# -------------------------------------------------------------------------

  # The purpose of these barplots is to compare the number of genes up and 
  # downregulated by person, using the DESeq2 population comparisons, and 
  # taking the intersection of the genes up/down regulated in all indivdiual
  # results. 
  
  
# Get the number of genes up/down per DESeq comparison
  table_up <- table_sig_genes(deseq_res, 'up')
  table_down <- table_sig_genes(deseq_res, 'down')
  
# This uses the lists created when making the upset plots to pull the number
# of genes up/down regulated in all individuals' datasets. 
  up_genes_9 <- num_genes_in_9(upset_mat_up_n, 'Up')
  down_genes_9 <- num_genes_in_9(upset_mat_down_n, 'Down')
  
# Bind the two above tables into one for graphing
  table_de <- rbind(table_up, table_down, up_genes_9, down_genes_9)
  
  title <- bquote(.(seq_mode) * '-seq after IFN' *beta* 'Treatment' * .(seq_mode))
  p <- bargraph_deres(table_de, title)
  p  
  ggsave(filename = paste0(savedir_figs, seq_mode, 'seq_barplot_de.png'),
         p,
         height = 6, width = 8)

  