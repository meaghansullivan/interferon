
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Purpose -----------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

  # This script cleans the RNA-seq data from the counting step and saves
  # the clean data for use in downstream scripts. 


# -------------------------------------------------------------------------
# Set environment  --------------------------------------------------------
# -------------------------------------------------------------------------

  library(tidyverse)
  
  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  path
  setwd(paste0(path, '/..'))
  getwd()
  source('intrahuman_ifnb/utility_functions.R')


# -------------------------------------------------------------------------
# Import data  ------------------------------------------------------------
# -------------------------------------------------------------------------

  # In this section, I will also select the longest isoform. 
  counts_import <- read_tsv(paste0('proseq/data/counts/fixed_genes_IFNB_Co',
                                   'unts_9_29_25_counts.txt'), 
                            skip = 1)
  
  # alt_alias <- read_csv('population_identifiers.csv')
  longest_isoform <- load_longest_isoform_annotation()


# -------------------------------------------------------------------------
# Simplify counts to only the longest isoform -----------------------------
# -------------------------------------------------------------------------

  counts <- get_longest_isoform(counts_import, longest_isoform)
  counts <- clean_col_names(counts, 'population')
  # in the next section, will need only rows with no NAs
  counts <- counts[complete.cases(counts),]
  

# -------------------------------------------------------------------------
# Remove three low-quality samples ----------------------------------------
# -------------------------------------------------------------------------

  drop <- c('BSA_USA-3-DS_2', 'BSA_USA-2_2', 'BSA_USA-1-DS_2',
            'IFN_USA-3-DS_2', 'IFN_USA-2_2', 'IFN_USA-1-DS_2')
  counts <- counts[, !colnames(counts) %in% drop]
  

# -------------------------------------------------------------------------
# Write out clean data ----------------------------------------------------
# -------------------------------------------------------------------------

  write_csv(counts, 
            file = 'proseq/data/counts/pop_proseq_cleaned_counts.csv')
