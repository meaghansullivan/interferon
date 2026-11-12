
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
  counts_import <- read_tsv(paste0('rnaseq/data/counts/interhuman_ifnb_',
                                   'rnaseq_featureCounts_transcript_rna',
                                   'seq.txt'))
  
  longest_isoform <- load_longest_isoform_annotation()


# -------------------------------------------------------------------------
# Simplify counts to only the longest isoform -----------------------------
# -------------------------------------------------------------------------

  counts <- get_longest_isoform(counts_import, longest_isoform)
  counts <- clean_col_names(counts, 'population')
  # in the next section, will need only rows with no NAs
  counts <- counts[complete.cases(counts),]
  
  print(colnames(counts))
  print(counts[1:3, ])

# -------------------------------------------------------------------------
# Write out clean data ----------------------------------------------------
# -------------------------------------------------------------------------

  write_csv(counts, 
            file = 'rnaseq/data/counts/pop_rnaseq_cleaned_counts.csv')
