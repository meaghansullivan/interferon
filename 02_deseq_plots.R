

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

  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  path
  setwd(paste0(path, '/..'))
  getwd()

  library(tidyverse)
  library(ComplexHeatmap)
  source(paste0('/Users/meco9877/Documents/Projects/ifnb_regression/intrahuman',
                '_ifnb/deseq2_helper_functions.R'))



# -------------------------------------------------------------------------
# Establish directories for loading and saving data -----------------------
# -------------------------------------------------------------------------

  
  if (seq_mode == 'PRO') {
    savedir_figs <- 'proseq/figures/deseq_figs/'
    savedir_data <- 'proseq/data/deseq_res/'
    data_load <- 'proseq/data/deseq_res/pop_pro_dt_gene_res.csv'
    
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
  
  # We had an expectation that the DESeq2 population comparisons would be 
  # very inclusive in its gene list. With the comparison, a gene does not 
  # have to be differentially expressed in all people for it to be called 
  # as differentially expressed.
  
  # In contrast, looking for the intersection of all individual results is 
  # very conservative. 
  
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

  

# -------------------------------------------------------------------------
# Export genes/features sig in 8/9 people ---------------------------------
# -------------------------------------------------------------------------

  # For up- and down-regulated genes, change to upset format
  upset_combo_up <- make_comb_mat(upset_mat_up_n)
  upset_combo_down <- make_comb_mat(upset_mat_down_n)
  
  # For each of the possible combinations of 8 people with one left out, I 
  # want to get that individual and all the features associated with them that 
  # are not differential in that individual, but are in the population. 
  
  # Get all the possible combinations of 8-'1' and 1-'0'
  codes_for_8_up <- comb_name(upset_combo_up[colSums(upset_combo_up) == 8])
  codes_for_8_down <- comb_name(upset_combo_down[colSums(upset_combo_down) == 8])
  
  # df of each gene sig in 8/9 people, and the one in whom it wasnt
  for_regression_up <- get_features_regressions(upset_mat_up_n, 
                                                upset_combo_up, 
                                                codes_for_8_up)
  
  for_regression_down <- get_features_regressions(upset_mat_down_n, 
                                                  upset_combo_down, 
                                                  codes_for_8_down)
  
  all_genes_for_regression <- rbind(for_regression_up, for_regression_down)
  colnames(all_genes_for_regression) <- c('Feature', 'Target_individual')
  
  write_csv(all_genes_for_regression,
            file = paste0(savedir_data, str_to_lower(seq_mode), 'seq_genes_sig_in8_ppl.csv'))


# -------------------------------------------------------------------------
# Export genes/features sig in 1/9 people ---------------------------------
# -------------------------------------------------------------------------
  
  # For each person who has some # of genes DT in only them, I 
  # want to get that individual and all the features associated with them that 
  # are only differential in that individual
  
  # Get all the possible combinations of 8-'1' and 1-'0'
  codes_for_1_up <- comb_name(upset_combo_up[colSums(upset_combo_up) == 1])
  codes_for_1_down <- comb_name(upset_combo_down[colSums(upset_combo_down) == 1])
  
  # df of each gene sig in 8/9 people, and the one in whom it wasnt
  for_regression_up <- get_features_regressions(upset_mat_up_n, 
                                                upset_combo_up, 
                                                codes_for_1_up, 
                                                num_de = '1')
  
  for_regression_down <- get_features_regressions(upset_mat_down_n, 
                                                  upset_combo_down, 
                                                  codes_for_1_down, 
                                                  num_de = '1')
  
  all_genes_for_regression <- rbind(for_regression_up, for_regression_down)
  colnames(all_genes_for_regression) <- c('Feature', 'Target_individual')
  
  write_csv(all_genes_for_regression,
            file = paste0(savedir_data, str_to_lower(seq_mode), 'seq_genes_sig_in1_ppl.csv'))
  
  
  
  
  
# -------------------------------------------------------------------------
# Num genes significant in given num of people ----------------------------
# -------------------------------------------------------------------------

  
  setsize_summary <- genes_per_setsize(upset_combo_up, upset_combo_down, seq_mode)
  write_csv(setsize_summary, 
            file = paste0(savedir_data, 'pop_', seq_mode, 'seq_de_summary.csv'))
  
  plot_genes_per_setsize(setsize_summary, seq_mode, 'count')
  ggsave(filename = paste0(savedir_figs, seq_mode, 
                           'seq_barplot_count_per_setsize.png'),
         height = 6, width = 10)
  
  plot_genes_per_setsize(setsize_summary, seq_mode, 'percent')
  ggsave(filename = paste0(savedir_figs, seq_mode, 
                           'seq_barplot_percent_per_setsize.png'),
         height = 6, width = 10)
  
  


# -------------------------------------------------------------------------
# New upset plot function  ------------------------------------------------
# -------------------------------------------------------------------------

  ifn_resp <- read_csv('~/Documents/Reference/resources/gene_sets/hallmark_ifnalpha_responsive.v2024.1.Hs.csv')
  deseq_res$Class <- ifelse(deseq_res$gene_id %in% ifn_resp$genes, 
                            'IFN-responsive', 
                            'Other')

  upset_up <- annotated_upset_deseq_res(deseq_res, "up", min_num_genes = 16, 
                                       title = 'Upregulated with 1hr IFNb Treatment')
  upset_up
  ggsave(filename = paste0(savedir_figs, str_to_lower(seq_mode), 
                           'seq_annot_upset_up_narrow.png'),
         upset_up,
         height = 7, width = 9)
  
  upset_down <- annotated_upset_deseq_res(deseq_res, "down", min_num_genes = 12, 
                                          title = 'Downregulated with 1hr IFNb Treatment')
  upset_down
  ggsave(filename = paste0(savedir_figs, str_to_lower(seq_mode), 
                           'seq_annot_upset_down.png'),
         height = 7, width = 13)
  
  

  