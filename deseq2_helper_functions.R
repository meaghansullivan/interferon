

# -------------------------------------------------------------------------
# Load the longest isoform ------------------------------------------------
# -------------------------------------------------------------------------

  # This function loads an annotation file with information about which
  # transcript isoform is the longest. It then cleans the names in the 
  # Gene and Transcript columns 

load_longest_isoform_annotation <- function() {
  longest_isoform <- read_tsv(
    paste0('~/Documents/Reference/resources/annota',
           'tions/hg38_diff53prime_refseq_longisof',
           '_noputatives_transcripts.bed'),
    col_names = FALSE
  )
  
  colnames(longest_isoform) <- c(
    'Chr', 'Start', 'End', 'Isoforms', 'Blank', 'Strand', 'Longest'
  )

  longest_isoform$Transcript <- sub('.+:', '', longest_isoform$Longest)
  longest_isoform$Transcript <- sub('\\.\\d', '', longest_isoform$Transcript)
  longest_isoform$Gene <- sub(':.+', '', longest_isoform$Longest)
  
  key <- longest_isoform[c('Gene', 'Transcript', 'Chr')]
  return(key)
}



# Select longest isoform --------------------------------------------------

  # This function joins the "key" made in `get_longest_isoform_annotation()`
  # to the input data such that the key is at the front of the dataset. 

get_longest_isoform <- function (data, longest_isoform) {
  
  if (grepl(':', data[1])){
    data$Transcript <- sub('.*:(.*)\\.\\d*', '\\1', data$Geneid)
    data$Geneid <- sub(':.*', '', data$Geneid)
  }else {
    colnames(data)[1] <- 'Transcript'
  }
  
  counts_longest_isoform <- left_join(longest_isoform, data)
  
  return(counts_longest_isoform)
}



# -------------------------------------------------------------------------
# Clean column names ------------------------------------------------------
# -------------------------------------------------------------------------

  # This takes the input data and removes text in the column names that 
  # comes from the file annotations and sequencing information. 
  # It also fixes inconsistencies in column names (i.e. is the time series
  # time point in minutes or hours). 

  # Data was given pseudonyms by the person who prepared the samples. I 
  # change those back to country of origin of the cell line used.

clean_col_names <- function (data, dataset) {
  
  if(dataset == 'ts'){
    
    colnames(data) <- sub('.sorted.bam', '', colnames(data))
    colnames(data) <- sub('_combo', '_batch1', colnames(data))
    colnames(data) <- ifelse(grepl('0m|90m|150m|180m', colnames(data)) &
                               grepl('D21|T21', colnames(data)) &
                               grepl('_S\\d+_L00\\d', colnames(data)),
                             paste0(colnames(data), '_batch2'),
                             colnames(data))
    colnames(data) <- ifelse(!grepl('batch|Gene|Transcript|Chr|Length', 
                                    colnames(data)),
                             paste0(colnames(data), '_batch1'),
                             colnames(data))
    colnames(data) <- sub('_S\\d+_L00\\d', '', colnames(data))
    colnames(data) <- sub('_RNA', '', colnames(data))
    colnames(data) <- sub('_IFN', '', colnames(data))
    colnames(data) <- sub('Eric', 'D21', colnames(data))
    colnames(data) <- sub('Ethan', 'T21', colnames(data))
    
  } else if (dataset == 'population'){
    
    colnames(data) <- sub('.mmfilt.sorted.bam', '', colnames(data)) # PRO
    colnames(data) <- sub('.sorted.bam', '', colnames(data)) # RNA
    colnames(data) <- sub('RNA.', '', colnames(data)) # RNA
    colnames(data) <- sub('./PRO-', '', colnames(data)) # PRO
    colnames(data) <- gsub('\\.', '_', colnames(data)) # RNA
    colnames(data) <- gsub('-', '_', colnames(data)) # PRO
    colnames(data) <- gsub('IFNB', 'IFN', colnames(data)) # PRO
    
    colnames(data) <- sub('Sengbe', 'Sierra Leone', colnames(data))
    colnames(data) <- sub('Khaondo', 'Kenya', colnames(data))
    colnames(data) <- sub('Niyilolawa', 'Nigeria', colnames(data))
    colnames(data) <- sub('Pedro', 'Peru', colnames(data))
    colnames(data) <- sub('Srivathani', 'Sri Lanka', colnames(data))
    colnames(data) <- sub('ChenChao', 'China', colnames(data))
    colnames(data) <- sub('Dave', 'USA-1-DS', colnames(data))
    colnames(data) <- sub('Eric', 'USA-2', colnames(data))
    colnames(data) <- sub('Ethan', 'USA-3-DS', colnames(data))
    
    data <- data[, !grepl('DR|Ursula', colnames(data))]
    
  }else {
    print('need dataset type')
  }

  return(data)
}



# -------------------------------------------------------------------------
# Get control (nonChr21) genes  -------------------------------------------
# -------------------------------------------------------------------------

  # This function is used when getting size-factors for DESeq2 when 
  # comparing sample from people with or without Down syndrome

getControlGenes <- function(counts){
  
  chr21genes <- counts[counts$Chr == 'chr21', ]$Gene
  
  # This creates a vector of T/F, whether every gene (i.e. line) in dfChr21 is 
  # a normal (control) gene or not-normal (on chr 21)
  control_genes <-ifelse(rownames(dds) %in% chr21genes, FALSE, TRUE)
}


# -------------------------------------------------------------------------
# Normalize by chr21 ------------------------------------------------------
# -------------------------------------------------------------------------

  # This function take in a list of genes that are on chr21 and uses the 
  # information to create normalization factors that reflect the expected 
  # ploidy of each cell type 

normChr21 <- function(dds, metadata, chr21genes){
  
  # Create a vectory with ploidy #, one entry for each gene in deseq_format. 
  # Need one for T21 samples and one for D21 samples 
  T21 <- ifelse(rownames(dds) %in% chr21genes, 3, 2)
  D21 <- rep(2, nrow(dds))
  
  # Now I can get normalization factors 
  normalizationfactors <- matrix(do.call(cbind, mget(paste0(metadata$Genotype))),
                                 ncol=ncol(dds), nrow=nrow(dds),
                                 dimnames=list(
                                   rownames(dds),colnames(dds)
                                 ))
  
  #divide by 2 = baseploidy, because typically there are 2 copies of gene per cell
  normalizationfactors <- normalizationfactors/2 
  return(normalizationfactors)
}


# -------------------------------------------------------------------------
# Run DESeq2 step-by-step -------------------------------------------------
# -------------------------------------------------------------------------

runDESeq2 <- function(dds, control_genes, normalizationfactors, savedir){
  
  # This function takes in the DESeq2 object, a list of control genes that 
  # can be used for estimateing size factors and the normalization factors 
  # from `normChr21()` and performs DESeq2, but doing each step independently. 
  # It is necessary to do each step independently so that the custom 
  # normalization factors and control genes can be used. 
  
  # Intermediate files are automatically saved to the `save_dir`
  
  #1. Estimate the Size Factor
  dds_normfactors <-estimateSizeFactors(
    dds,
    normMatrix=normalizationfactors,
    controlGenes=control_genes) 
  
  #2. Estimate Dispersion
  dds_dispersion <- estimateDispersionsGeneEst(dds_normfactors) 
  png(filename = paste0(savedir, 'dispersion_estimates.png'), 
      height = 480, width = 672, units = 'px')
  par(
    cex = 1.2, 
    cex.axis = 1.2,
    cex.lab = 1.4, 
    cex.main = 1.4 
  )
  plotDispEsts(dds_dispersion)
  dev.off()
  
  #3. Fit linear model
  dds_fit <- estimateDispersionsFit(dds_dispersion) 
  png(filename = paste0(savedir, 'dispersion_fit.png'), 
      height = 480, width = 672, units = 'px')
  par(
    cex = 1.2, 
    cex.axis = 1.2, 
    cex.lab = 1.4, 
    cex.main = 1.4 
  )
  plotDispEsts(dds_fit)
  dev.off()
  
  # Get shrunken dispersion estimates by shrinking gene-wise dispersion
  # estimates towards the trend 
  dds_map <- estimateDispersionsMAP(dds_fit) 
  png(filename = paste0(savedir, 'dispersion_shrunk.png'), 
      height = 480, width = 672, units = 'px')
  par(
    cex = 1.2,      
    cex.axis = 1.2, 
    cex.lab = 1.4,   
    cex.main = 1.4   
  )
  plotDispEsts(dds_map)
  dev.off()
  
  #4. Hypothesis testing with Wald Test
  dds_final <- nbinomWaldTest(dds_map,
                              betaPrior = FALSE) 
  
  return(dds_final)
}



# -------------------------------------------------------------------------
# Create complex contrasts, used to extract results from DESeq obj --------
# -------------------------------------------------------------------------


make_contrasts <- function(dds, person = NULL, weighted = FALSE){
  #1. Get coefficient names 
  coef_names <- resultsNames(dds)
  
  #2. Identify main treatment effect and interaction terms 
  treat_main <- grep("^Treatment_IFN_vs_BSA$", coef_names, value=TRUE)
  interactions <- grep("^TreatmentIFN.Person", coef_names, value=TRUE)
  
  #3. Start contrast vector (one value per coefficient)
  con_vector <- rep(0, length(coef_names))
  names(con_vector) <- coef_names
  
  #4. How to get effects for a specific person 
  if (!person == 'population'){
    # Baseline person just use main term
    if (person == "China"){
      con_vector[treat_main] <- 1
    } else {
      interaction_term <- paste0("TreatmentIFN.Person", person)
      con_vector[treat_main] <- 1
      con_vector[interaction_term] <- 1
    }
  } else {
    #5. Get contrasts (weighted or not) for the population 
    # For each person, get the # of samples 
    people <- levels(colData(dds)$Person)
    person_counts <- table(colData(dds)$Person)
    
    # Weight - 1/ # of people 
    weights <- rep(1/length(people), length(people))
    names(weights) <- people
    
    # The main term is always 1  
    con_vector[treat_main] <- sum(weights)
    
    # add interaction term for everyone but the baseline 
    for (p in people[-1]) {
      term <- paste0("TreatmentIFN.Person", p)
      con_vector[term] <- weights[[p]]
    }
  }
  return(con_vector)
}


# -------------------------------------------------------------------------
# Clean results from DESeq2 object ----------------------------------------
# -------------------------------------------------------------------------

clean_res <- function(results_df, res_name){
  temp <- as.data.frame(results_df) %>% 
    rownames_to_column('gene_id') %>% 
    tibble() %>% 
    mutate(Sig = ifelse(padj <= 0.05, 'padj <= 0.05', 'ns')) %>% 
    drop_na() %>% 
    mutate(Comparison = res_name)
}



# -------------------------------------------------------------------------
# Get number differential genes per comparison ----------------------------
# -------------------------------------------------------------------------


table_sig_genes <- function(deseq_res, fc_direction) {
  
  if(fc_direction == 'up'){
    sig <- deseq_res %>% 
      filter(Sig == 'padj <= 0.05' & log2FoldChange > 0) #& !grepl('Population', Comparison))
    table <- sig %>% 
      group_by(Comparison) %>% 
      summarise(sig_genes = n()) %>% 
      mutate(Category = 'Up')
  } else (fc_direction == 'down') {
    sig <- deseq_res %>% 
      filter(Sig == 'padj <= 0.05' & log2FoldChange < 0) #& !grepl('Population', Comparison))
    table <- sig %>% 
      group_by(Comparison) %>% 
      summarise(sig_genes = n()) %>% 
      mutate(Category = 'Down')
  }
  
  return(table)
}


# -------------------------------------------------------------------------
# Creates a bar graph representing # of genes up/down by comparison -------
# -------------------------------------------------------------------------

bargraph_deres <- function(table_de_res, user_subtitle){
  
  # Create factor levels so that everything plots in the correct order
  table_de_res$Comparison <- factor(table_de_res$Comparison, 
                                levels = c('Population',
                                           'All Individuals',
                                           'China', 
                                           'Kenya', 
                                           'Nigeria', 
                                           'Peru', 
                                           'Sierra.Leone', 
                                           'Sri.Lanka', 
                                           'USA.1.DS', 
                                           'USA.2', 
                                           'USA.3.DS'))
  table_de_res$Category <- factor(table_de_res$Category, 
                                  levels = c('Down', 'Up', 
                                             'Down in All', 'Up in All'))
  
  p <- ggplot(table_de_res, aes(x = Comparison, y = sig_genes, fill = Category)) + 
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.9) + 
    theme_bw() + 
    labs(x = " ", y = "Genes Significantly Changed\n", 
         title = expression("Genes Changing During IFN" * beta * " Treatment"), 
         subtitle = user_subtitle) +
    theme(plot.title = element_text(face="bold", size = 24, hjust = 0.5), 
          plot.subtitle = element_text( size = 20, hjust = 0.5),
          axis.title.y = element_text(face="bold", size = 20),
          legend.title = element_text(face="bold", size = 18),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
          axis.text.y = element_text(size = 18)) + 
    scale_fill_manual(values = c( '#E1E5F2', '#7E90C8', '#4EA43D', '#093824'))
  
  return(p)
}



# -------------------------------------------------------------------------
# Creates a MA plot and highlights several IFN-responsive genes -----------
# -------------------------------------------------------------------------

ma_plot <- function(deseq_res, comparisons, assay){
  
  # Get results for specific comparison
  df_temp <- deseq_res %>%
    filter(Comparison == comparisons)
  
  # These are a small selection of IFN-responsive genes
  df_label <- df_temp %>%
    filter(gene_id %in% c('IFIT2', 'IFI6', 'IFIT1',
                          'IFIT3', 'RTP4', 'OASL', 'USP18', 'IRF7', 'OAS1'))
  
  # Make df's of what's up and down
  df_temp_up <- df_temp %>% filter(log2FoldChange > 0 & padj <= 0.05)
  df_temp_down <- df_temp %>% filter(log2FoldChange < 0 & padj <= 0.05)
  
  # Get sizes of what's up and down
  num_up <- dim(df_temp_up)[1]
  print(num_up)
  num_down <- dim(df_temp_down)[1]
  print(num_down)
  
  # Plot
  plot_res <- ggplot(df_temp, aes(y=log2FoldChange, x=log10(baseMean),
                                  color = Sig)) +
    geom_hline(yintercept = 0, color = 'gray20') + 
    geom_point(alpha = 0.4) +
    scale_color_manual(values = c('gray70', 'red4')) +
    theme_light() +
    ggtitle(paste0('DE of Genes in ', assay, '-seq: ', res_name)) + #,
            # subtitle = res_name)+
    # ylim(-9, 9) +
    # xlim(0, 5) +
    geom_point(data = df_label, aes(y=log2FoldChange, x=log10(baseMean)),
               shape = 1, color = 'black', size = 3) +
    geom_text(data = df_label, aes(y=log2FoldChange, x=log10(baseMean), label = gene_id),
              color = 'black', check_overlap = TRUE, vjust = 1.5) +
    annotate("text", 
             x = -Inf, 
             y = Inf, 
             hjust = -0.02, 
             vjust = 1.25,
             label = paste("Up:", num_up), hjust = 1, vjust = 1, 
             fontface = "bold") +
    annotate("text", 
             x = -Inf, 
             y = -Inf, 
             hjust = -0.02, 
             vjust = -1.25,
             label = paste("Down:", num_down), hjust = 1, vjust = 0, 
             fontface = "bold")

  return(plot_res)
}


# -------------------------------------------------------------------------
# Create a list ready for the upset function ------------------------------
# -------------------------------------------------------------------------

create_upset_list <- function(deseq_results, direction, population){
  
  if (population == 'no'){
    deseq_results <- deseq_results[ !deseq_results$Comparison == 'Population', ]
  }
  
  people_list <- unique(deseq_results$Comparison)
  upset_list <- list()
  
  if (direction == 'up'){
    
    for (person in people_list) {
      print(person)
      upset_list[[person]] <- deseq_results[grepl(person, deseq_results$Comparison) & 
                                             deseq_results$Sig == 'padj <= 0.05' & 
                                             deseq_results$log2FoldChange > 0, ]$gene_id
      
    }
    
  } else if (direction == 'down') {

    for (person in people_list) {
      print(person)
      upset_list[[person]] <- deseq_results[grepl(person, deseq_results$Comparison) & 
                                             deseq_results$Sig == 'padj <= 0.05' & 
                                             deseq_results$log2FoldChange < 0, ]$gene_id
    
    }
  }
  
  upset_matrix <- list_to_matrix(upset_list)
  
  return(upset_matrix)
}



# -------------------------------------------------------------------------
# Make the upset plot, including categories that are above the min --------
# -------------------------------------------------------------------------

make_upset_plot <- function (upset_matrix, min_num_genes, user_title) {
  
  upset_combo <- make_comb_mat(upset_matrix)
  
  upset_p <- UpSet(
    upset_combo [comb_size(upset_combo ) > min_num_genes],
    comb_order = order(-comb_size(upset_combo [comb_size(upset_combo ) > min_num_genes])),
    column_title = user_title
  )
  
  draw(upset_p)
  
}
