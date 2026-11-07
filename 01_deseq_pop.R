

# -------------------------------------------------------------------------
# PURPOSE -----------------------------------------------------------------
# -------------------------------------------------------------------------

  # The purpose of this script is to perform differential expression
  # analysis using DESeq. This script can be used for analyzing RNA-seq and 
  # nascent-RNA seq (PRO-seq). 

  # Here it is being applied to PRO-seq and RNA-seq data generated from LCLs 
  # originating from a small population of individuals. 
  # This data set consists of LCLs from 9 individuals located across Africa, 
  # Asian, North America, and South America. 

  # The experiment consists of a 1 hour interferon-beta stimulation (100ng/uL)
  # and a BSA control. 


# -------------------------------------------------------------------------
# User variables ----------------------------------------------------------
# -------------------------------------------------------------------------

  seq_mode <- 'RNA' # options = c('PRO', 'RNA')
  
  deseq_design <- ~ Treatment + Person + Treatment:Person



# -------------------------------------------------------------------------
# Set environment ---------------------------------------------------------
# -------------------------------------------------------------------------

  # Import libraries and source custom functions
  library(tidyverse)
  library(DESeq2)
  source(paste0('/Users/meco9877/Documents/Projects/ifnb_regression/intrahuman',
                '_ifnb/deseq2_helper_functions.R'))

  
  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  path
  setwd(paste0(path, '/..'))
  getwd()



# -------------------------------------------------------------------------
# Establish paths for loading and saving  ---------------------------------
# -------------------------------------------------------------------------
    
  if (seq_mode == 'PRO') {
    savedir_figs <- 'proseq/figures/deseq_figs/'
    savedir_data <- 'proseq/data/deseq_res/'
    data_load <- 'proseq/data/counts/pop_proseq_cleaned_counts.csv'
    
  } else {
    savedir_figs <- 'rnaseq/figures/deseq_figs/'
    savedir_data <- 'rnaseq/data/deseq_res/'
    data_load <- 'rnaseq/data/counts/pop_rnaseq_cleaned_counts.csv'
  }
  
  print(savedir_figs)
  print(savedir_data)
  print(data_load)
  

  
# -------------------------------------------------------------------------
# Read in data ------------------------------------------------------------
# -------------------------------------------------------------------------

  
  counts_import <- read_csv(data_load)
  print(counts_import[1:3, 1:9])
  
# -------------------------------------------------------------------------
# Create data and metadata for deseq --------------------------------------
# -------------------------------------------------------------------------
  
# Drop non-count columns (chr#, gene length, isoform names, etc)
  if(seq_mode == 'PRO') {
    counts_df <- data.frame(counts_import[-c(1:8)])  
  } else {
    counts_df <- data.frame(counts_import[-c(1:4)])
  }
  
  rownames(counts_df) <- counts_import$Gene
  print(counts_df[1:4, 1:6])
  
# Create metadata 
  metadata <- data.frame(SampleIDs = colnames(counts_df))
  separated_info <- strsplit(metadata$SampleIDs, split = '_')
  new_cols <- do.call(rbind, separated_info)
  colnames(new_cols) <- c('Treatment', 'Person', 'Rep')
  metadata <- cbind(metadata, new_cols)
  
# Set factor levels for metadata
  metadata$Treatment <- factor(metadata$Treatment, levels = c("BSA", "IFN"))
  metadata$Person <- factor(metadata$Person, levels = unique(metadata$Person))
  metadata$Genotype <- ifelse(metadata$Person %in% c('USA.1.DS', 'USA.3.DS'), 
                              'T21', 'D21')
  metadata$Genotype <- factor(metadata$Genotype, levels = unique(metadata$Genotype))
  # Group is used for custom colors when graphing data 
  metadata$Group <- factor(paste0(metadata$Genotype, ' - ', metadata$Treatment))
  # Sibling is used for custom colors when graphing data 
  metadata$Sibling <- factor(ifelse(metadata$Person %in% c('USA.2', 'USA.3.DS'), 
                                    'Related', 'Unrelated'))
  metadata$Sibling <- paste0(metadata$Sibling, ' - ', metadata$Treatment)

  # Look at metadata, just double check that everything is correct   
  print(metadata[1:4, ])

  
# -------------------------------------------------------------------------
# Create DESeq2 Object ----------------------------------------------------
# -------------------------------------------------------------------------
  
# Create DESeq object
  dds <-DESeqDataSetFromMatrix(countData = counts_df, 
                               colData = metadata, 
                               design = deseq_design) 
  
  
# Filter out lowly expressed genes. PRO-seq only have two replicates and RNA-seq
# has 3 replicates, so adjust read count cut-off to reflect this. 
  
  if (seq_mode == 'pro'){
    smallestGroupSize <- 2
    keep <- rowSums(counts(dds) >= 15) >= smallestGroupSize 
    dds <- dds[keep,] 
    
  } else {
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds) >= 40) >= smallestGroupSize 
    dds <- dds[keep,] # go from 27,350 to 13,475
    
  }


  
# -------------------------------------------------------------------------
# Get chr21-aware norm factors --------------------------------------------
# -------------------------------------------------------------------------

  # This population includes two individuals with Down syndrome, a genetic
  # condition that occurs when chr21 is triplicated. Therefore, we know that
  # chr21 genes will have higher abundance in the dataset, and chr21 genes 
  # should not be used to calculate size factors for normalizing library size
  
  chr21_genes <- counts_import[counts_import$Chr == 'chr21', ]$Gene
  control_genes <- getControlGenes(counts_import, dds)
  normFactors <- normChr21(dds, metadata, chr21_genes)
  
  
# -------------------------------------------------------------------------
# Run DESeq2 --------------------------------------------------------------
# -------------------------------------------------------------------------

  # DESeq has 4 steps:
  #   1. Estimate the size factor
  #   2. Estimate dispersion
  #   3. Fit linear model
  #   4. Hypothesis testing with Wald Test    
  
  dds <- runDESeq2(dds, control_genes, normFactors, savedir_figs)
  
  
  
# -------------------------------------------------------------------------
# Extract VST counts ------------------------------------------------------
# -------------------------------------------------------------------------

# Transform counts to be used in downstream visualizations
  vs_counts <- vst(dds, blind=FALSE)
  
# Get a matrix of counts using assay 
  vs_matrix <- assay(vs_counts)
  
# Save counts 
  vs_df <- as.data.frame(vs_matrix) %>% rownames_to_column('Gene')
  write_csv(vs_df, 
            file = paste0(savedir_data, 'pop_', str_to_lower(seq_mode), 'seq_counts_vst.csv'))  

  
# -------------------------------------------------------------------------
# Get DESeq results: POPULATION -------------------------------------------
# -------------------------------------------------------------------------
  
# Get contrast matrix - calls the custom function "make_contrasts()"
  ifnb_bsa_population_design <- make_contrasts(dds, 
                                               person = 'population', 
                                               weighted = TRUE)
  # Get results: 
  ifnb_bsa_pop_res <- lfcShrink(dds, 
                                contrast = ifnb_bsa_population_design, 
                                type='ashr')
  df_res <- clean_res(ifnb_bsa_pop_res, 'Population')
  
  
# -------------------------------------------------------------------------
# Get DESeq results: INDIVIDUALS ------------------------------------------
# -------------------------------------------------------------------------

# List of all the individuals, I will loop through this list to get the 
# DESeq2 results for each person
  people_list <- unique(metadata$Person)
  
  for (person in people_list) {
    print(person)
    # make the appropriate contrast matrix
    contrast_temp <- make_contrasts(dds, person = person)
    
    # Get results:
    res_temp <- lfcShrink(dds,
                          contrast = contrast_temp,
                          type='ashr')
    res_df <- clean_res(res_temp, person)

    # Add results to master df
    df_res <- rbind(df_res, res_df) 
  }  

  
# -------------------------------------------------------------------------
# Write DESeq2 Results ----------------------------------------------------
# -------------------------------------------------------------------------
  
  write_csv(df_res,
            file = paste0(savedir_data, 'pop_', str_to_lower(seq_mode),
                          ifelse(seq_mode == 'PRO', '_dt_gene_res.csv', '_dt_res.csv')))

  
# -------------------------------------------------------------------------
# Make MA Plots -----------------------------------------------------------
# -------------------------------------------------------------------------
  
  for (res_name in unique(df_res$Comparison)){
    print(res_name)
    p <- ma_plot(df_res, res_name, seq_mode)
    print(p)
    ggsave(filename = paste0(savedir_figs, res_name, '_', str_to_lower(seq_mode), 
                             '_maplot.png'),
           p,
           height = 5, width = 6)
  }
  
  
  
# -------------------------------------------------------------------------
# PCA plots ---------------------------------------------------------------
# -------------------------------------------------------------------------
  
  plotPCA(vs_counts, intgroup=c("Group"), 
          ntop = 500) + 
    theme_bw() + 
    ggtitle(ifelse(seq_mode == "PRO", 
                   'PCA of Nascent RNA sequencing',
                   'PCA of RNA sequencing'), 
            subtitle = 'Colored by Treatment and Genotype') + 
    theme(
      axis.text = element_text(size = 14), 
      axis.title = element_text(size = 16), 
      plot.title = element_text(size = 18, hjust = 0.5), 
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    ) + 
    scale_color_manual(values = c('#66c2a4', '#238b45', '#fecc5c', '#fd8d3c'))
  ggsave(filename = paste0(savedir_figs, str_to_lower(seq_mode),
                           '_genes_PCA_treatment.png'),
         height = 4.5, width = 6)
  
  
  plotPCA(vs_counts, intgroup=c("Sibling"), 
          ntop = 500) + 
    theme_bw() + 
    ggtitle(ifelse(seq_mode == "PRO", 
                   'PCA of Nascent RNA sequencing',
                   'PCA of RNA sequencing'), 
            subtitle = 'Siblings Highlighted') + 
    theme(
      axis.text = element_text(size = 14), 
      axis.title = element_text(size = 16), 
      plot.title = element_text(size = 18, hjust = 0.5), 
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    ) + 
    scale_color_manual(values = c('#7B435B', '#BE89A0', '#767656', '#A3A380'))
  ggsave(filename = paste0(savedir_figs, str_to_lower(seq_mode), 
                           '_genes_PCA_siblings.png'),
         height = 4.5, width = 6)
  