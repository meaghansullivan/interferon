

# -------------------------------------------------------------------------
# PURPOSE -------------------------------------------------------------
# -------------------------------------------------------------------------

  # Helper functions for the population (9-individual) DESeq2 pipeline:
  # run DESeq2 with the chr21 ploidy correction (needed because 2 of
  # the 9 donors are trisomy 21), pull per-person/population contrasts
  # and results out of the DESeq2 object. Plot results (MA, bar graph, UpSet).

  # Sourced by: 01_deseq_pop.R, 02_deseq_plots.R 

  # Functions below are grouped as:
  #   1. DESeq2 run with chr21 ploidy correction
  #   2. Contrasts & extracting results from the DESeq2 object
  #   3. Summarizing significant genes
  #   4. Single-comparison plots (MA, bar graph)
  #   5. Multi-individual overlap (UpSet) plots and gene extraction
  #   6. Results vs. library depth


# -------------------------------------------------------------------------
# Set environment -------------------------------------------------------
# -------------------------------------------------------------------------

  library(ComplexUpset)
  library(tidyverse)


# =========================================================================
# 1. DESeq2 run with chr21 ploidy correction ==============================
# =========================================================================

# Get control (nonChr21) genes  --------------------------------------------

  # Returns a logical vector (one entry per gene in `dds`) that is FALSE 
  # for chr21 genes and TRUE for everything else.
  # Passed to `estimateSizeFactors(controlGenes = ...)` inside `runDESeq2()`
  # so that chr21 genes don't bias library-size normalization.

  getControlGenes <- function(counts, dds){
    chr21genes <- counts[counts$Chr == 'chr21', ]$Gene
  
    # This creates a vector of T/F, whether every gene (i.e. line) is
    # a normal (control) gene or not-normal (on chr 21)
    control_genes <-ifelse(rownames(dds) %in% chr21genes, FALSE, TRUE)
  }


# Normalize by chr21 ------------------------------------------------------

  # Builds a per-gene, per-sample normalization-factor matrix that encodes
  # expected ploidy (3 copies of chr21 genes for T21 samples, 2 copies
  # everywhere else). Given to DESeq2, size factors calcualtions account for 
  # the extra chr21 dosage in the two trisomy-21 individuals 

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
  
    #divide by 2 = normal human ploidy
    normalizationfactors <- normalizationfactors/2
    return(normalizationfactors)
  }


# Run DESeq2 step-by-step --------------------------------------------------

  # Runs the DESeq2 steps (size factors, dispersion estimation/fit/shrink,
  # Wald test) individually, which I am doing to account for trisomy 21.
  # Uses `getControlGenes()` and `normChr21()`, defined above. Saves
  # dispersion plots from each step.

  # Used in: 01_script.R, 05_script.R

  runDESeq2 <- function(dds, control_genes, normalizationfactors, savedir){
  
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
      cex = 1.2,       # overall text scaling
      cex.axis = 1.2,  # axis tick labels
      cex.lab = 1.4,   # x/y axis labels
      cex.main = 1.4   # title
    )
    plotDispEsts(dds_map)
    dev.off()
  
    #4. Hypothesis testing with Wald Test
  
    dds_final <- nbinomWaldTest(dds_map,
                                betaPrior = FALSE)
  
    return(dds_final)
  
  }


# =========================================================================
# 2. Contrasts & extracting results from the DESeq2 object ================
# =========================================================================

# Create complex contrasts -------------------------------------------------

  # Builds the numeric contrast vector passed to `results(dds, contrast =
  # ...)` for a given individual (their specific IFN-vs-BSA effect) or for
  # the population, which is the average effect across all 9 individuals.
  # There is an option for weighted sample count, although I don't end up 
  # using that here. 
  
  # The main effect is "China", since I let the factors be set alphabetically

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
  
      # Design is different if weighted or not
      if (weighted){
        # Weight = samples for that person / total # of samples
        weights <- person_counts / sum(person_counts)
      } else {
        # Weight - 1/ # of people
        weights <- rep(1/length(people), length(people))
      }
  
      # Name vector so can pull out weights easily (not necessary if weighted,
      # but necessary if unweighted)
      names(weights) <- people
  
      # Everyone contributes to main term, this is always 1
      con_vector[treat_main] <- sum(weights)
  
      # add interaction term for everyone but the baseline
      for (p in people[-1]) {
        term <- paste0("TreatmentIFN.Person", p)
        con_vector[term] <- weights[[p]]
      }
    }
    return(con_vector)
  }


# Clean results from DESeq2 object -----------------------------------------

  # Converts a `results()`/`lfcShrink()` output object into a tibble.
  # Add a Sig column (Significance) to make graphing easier later.
  # Adds a Comparison column, which is the person under question.
  
clean_res <- function(results_df, res_name){
  temp <- as.data.frame(results_df) %>%
    rownames_to_column('gene_id') %>%
    tibble() %>%
    mutate(Sig = ifelse(padj <= 0.05, 'padj <= 0.05', 'ns')) %>%
    drop_na() %>%
    mutate(Comparison = res_name)
}


# =========================================================================
# 3. Summarizing significant genes =========================================
# =========================================================================

# Get number of differential genes per comparison ---------------------------

  # Filters to get either up or down regulated genes, then turns that into 
  # a table of genes per comparison in that direction.

  table_sig_genes <- function(deseq_res, fc_direction) {
  
    if(fc_direction == 'up'){
      sig <- deseq_res %>%
        filter(Sig == 'padj <= 0.05' & log2FoldChange > 0) #& !grepl('Population', Comparison))
      table <- sig %>%
        group_by(Comparison) %>%
        summarise(sig_genes = n()) %>%
        mutate(Category = 'Up')
    } else if (fc_direction == 'down') {
      sig <- deseq_res %>%
        filter(Sig == 'padj <= 0.05' & log2FoldChange < 0) #& !grepl('Population', Comparison))
      table <- sig %>%
        group_by(Comparison) %>%
        summarise(sig_genes = n()) %>%
        mutate(Category = 'Down')
    }else {
      print('fc_direction must be either up or down')
    }
  
    return(table)
  }


# =========================================================================
# 4. Visualizations =======================================================
# =========================================================================

# Creates a bar graph of # genes up/down by comparison ----------------------

  # Bar graph of `table_sig_genes()` output. 

  bargraph_deres <- function(table_de_res, user_subtitle){
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


# Creates an MA plot and highlights several IFN-responsive genes ------------

  # MA plot (log2FoldChange vs. log10(baseMean)) for each individual,
  # colored by significance. I am specifying a small set of IFN-responsive
  # genes to plot. Maybe someday I should make the function take in a list 
  # of genes to use. These genes are circled and labeled. 

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
      theme_bw() +
      theme(
        plot.title = element_text(face="bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text( size = 14, hjust = 0.5),
        axis.title = element_text(face="bold", size = 14),
        legend.title = element_text(face="bold", size = 13),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      ) +
      ggtitle(paste0('DE of Genes in ', assay, '-seq'),
              subtitle = paste0(comparisons, ": IFN-beta vs Control")) +
      labs(x = expression(log[10](Average~Counts)),
           y = expression(log[2](Fold~Change)),
           color = 'Signficance') +
      geom_point(data = df_label, aes(y=log2FoldChange, x=log10(baseMean)),
                 shape = 1, color = 'black', size = 3) +
      geom_text(data = df_label, aes(y=log2FoldChange, x=log10(baseMean), label = gene_id),
                color = 'black', check_overlap = TRUE, vjust = 1.5,
                size = 4, fontface = 'bold') +
      annotate("text",
               x = c(-Inf, -Inf),
               y = c(Inf, -Inf),
               hjust = c(-0.02, -0.02),
               vjust = c(1.25, -1.25),
               label = c(paste("Up:", num_up),
                         paste("Down:", num_down)),
               fontface = "bold",
               size = 5)
    return(plot_res)
  }


# =========================================================================
# 5. Upset plots, and functions that use ComplexHeatmap ===================
# =========================================================================

# Create a list ready for the upset() function ------------------------------

  # Creates a matrix (gene x person) with T/F for if a gene is significant 
  # (padj <= 0.05) and moving in specified direction. 
  # This gets passed into `ComplexHeatmap::make_comb_mat()`/`UpSet()`. 
  # `population = 'no'` drops the pooled 'Population' comparison from the 
  # columns, because sometimes having Population is really boring. 

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
    } else {
      print('direction must be up or down')
    }
  
    upset_matrix <- list_to_matrix(upset_list)
  
    return(upset_matrix)
  }


# Make the UpSet plot, restricted to intersections above the minimum --------

  # Draws an UpSet plot (via `ComplexHeatmap`) from the matrix built by
  # `create_upset_list()`, keeping only intersections with more than
  # `min_num_genes` genes and ordering them largest-first.

  make_upset_plot <- function (upset_matrix, min_num_genes, user_title) {
  
    upset_combo <- make_comb_mat(upset_matrix)
  
    upset_p <- UpSet(
      upset_combo [comb_size(upset_combo ) > min_num_genes],
      comb_order = order(-comb_size(upset_combo [comb_size(upset_combo ) > min_num_genes])),
      column_title = user_title
    )
  
    draw(upset_p)
  
  }

# Count genes significant in all 9 individuals -------------------------------

  # Pulls out the size of the "all individuals" intersection (comb code
  # '111111111', i.e. all 9 columns) from an UpSet combination matrix and
  # returns it as a one-row tibble shaped to match `table_sig_genes()`
  # output (`Comparison = 'All Individuals'`), for appending to a
  # `bargraph_deres()` input table.

  num_genes_in_9 <- function (upset_matrix, direction) {
  
    upset_combo <- make_comb_mat(upset_matrix)
    num_genes <- length(extract_comb(upset_combo, '111111111'))
  
    temp <- tibble(Comparison = 'All Individuals',
                           sig_genes = num_genes,
                           Category = paste0(direction, ' in All'))
    return(temp)
  }


# Get the gene features behind a set of UpSet intersection codes ------------

  # For each `code` in `code_list` ( the 1/0 code, ie. '100000000'),
  # get the genes in that intersection from `upset_combo` and record
  # which individual(s) the code refers to. 
  # `num_de` specifies which codes to get, they are codes that add up 
  # to the given number. 
  # Could be used for any arbitrary number of genes, but because I'm focused
  # on 1 and 8, I specify that, and basically learn how to use match.arg

get_features_regressions <- function(upset_mat, upset_combo, code_list, num_de = c('1', '8')){

  num_de <- match.arg(num_de)

  # Create a tibble to hold data as it gets extracted
  features_df <- tibble()

  # Loop through all of the codes
  for (code in code_list){
    print(code)
    # split code into vector of 1's and 0's
    cols_code <- as.numeric(strsplit(code, "")[[1]])
    # Get the identity of the individual not significant by taking the logical
    # inverse of the code
    if (num_de == 1){
      name <- colnames(upset_mat)[as.logical(cols_code)]

    } else if (num_de == 8) {
      name <- colnames(upset_mat)[!as.logical(cols_code)]

    }


    features_temp <- data.frame(Feature = extract_comb(upset_combo, code),
                                Individual = name)
    features_df <- rbind(features_df, features_temp)
  }

  return(features_df)
}


# Summarize gene counts by intersection set size -----------------------------

  # Collapses an UpSet combination matrix's per-intersection up/down gene
  # counts (`combo_up`/`combo_down`, both from `make_comb_mat()`) down to
  # per-set-size totals (i.e. "how many genes are shared across exactly N
  # individuals"), with both raw counts and each size's percent of all
  # up/down genes.

  genes_per_setsize <- function( combo_up, combo_down, assay){
  
    temp_up <- comb_size(combo_up)
    temp_down <- comb_size(combo_down)
    temp_df_up <- tibble(combs = names(temp_up), num_up = temp_up)
    temp_df_down <- tibble(combs = names(temp_down), num_down = temp_down)
  
    temp_df <- full_join(temp_df_up, temp_df_down)
    temp_df[is.na(temp_df)] <- 0
  
    temp_df$`Number of Individuals` <- rowSums(do.call(rbind, strsplit(temp_df$combs, "")) == "1")
  
    temp_final <- temp_df %>%
      group_by(`Number of Individuals`) %>%
      summarise(`Up-regulated` = sum(num_up),
                `Down-regulated` = sum(num_down)) %>%
      mutate(`Percent Up` = `Up-regulated` / sum(`Up-regulated`),
             `Percent Down` = `Down-regulated` / sum(`Down-regulated`))
  
    return(temp_final)
  
  }


# Plot gene counts by intersection set size ----------------------------------

  # Bar graph of `genes_per_setsize()` output: number (or percentage, via
  # `type = 'count'`/`'percent'`) of up-/down-regulated genes. 
  # I pre-set colors depending on whether ths assay is PRO or RNA, to help 
  # visually show different assays when presenting.
  
  plot_genes_per_setsize <- function(summarized_df, assay, type){
  
    if (assay == 'PRO') {
      set_subtitle <- 'PRO-seq, 1hr of IFN-b Treatment'
      color_values <- c('#E1E5F2', '#7E90C8')
    } else if (assay == 'RNA') {
      set_subtitle <- 'RNA-seq, 3hr of IFN-b Treatment'
      color_values <- c('#BCCFF2', '#4A5EC8')
    } else {
      stop('Assay must be "PRO" or "RNA"')
    }
  
    if(type == 'count'){
      temp_plot <- summarized_df %>%
        select(`Number of Individuals`, `Up-regulated`, `Down-regulated`) %>%
        pivot_longer(cols = -`Number of Individuals`,
                     names_to = 'Regulation',
                     values_to = 'Number of Genes')
      y_val <- 'Number of Genes'
    } else if(type == 'percent'){
      temp_plot <- summarized_df %>%
        select(`Number of Individuals`, `Percent Up`, `Percent Down`) %>%
        pivot_longer(cols = -`Number of Individuals`,
                     names_to = 'Regulation',
                     values_to = 'Percentage of Genes')
      y_val <- 'Percentage of Genes'
    } else {
      print('Enter a type')
    }
  
  
    p <- ggplot(temp_plot,
                aes(x = `Number of Individuals`, y = .data[[y_val]], fill = Regulation)) +
      geom_bar(stat = 'identity', position = 'dodge') +
      scale_x_continuous(breaks = seq(0, 9, by = 1)) +
      ggtitle('Significant Genes By Number of Individuals',
              subtitle = set_subtitle) +
      theme_bw() +
      theme(
        plot.title = element_text(face="bold", size = 24, hjust = 0.5),
        plot.subtitle = element_text( size = 20, hjust = 0.5),
        axis.title = element_text(face="bold", size = 20),
        legend.title = element_text(face="bold", size = 18),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18)
      ) +
      scale_fill_manual(values = color_values)
  
    return(p)
  
  }

# Annotated UpSet plot, colored by IFN-responsive gene class -----------------

  # `ComplexUpset::upset()` plot of significant genes moving 'up' or 'down'
  # across individuals (always excludes 'Population'),
  # for intersections with at least `min_num_genes` genes. 
  # The bar plot part of the upset is stacked/colored by a `Class` column,
  # which here is IFN-responsive or not. 

  annotated_upset_deseq_res <- function(dataframe = de_newsig,
                              direction = c('up', 'down'),
                              min_num_genes = 10,
                              title){
  
    direction <- match.arg(direction)
    intersect_cols <- setdiff(unique(dataframe$Comparison), 'Population')
  
  
    if (direction == 'up') {
      sig_genes <- dataframe %>%
        filter(Comparison %in% intersect_cols,
               Sig == 'padj <= 0.05',
               log2FoldChange > 0)
    } else {
      sig_genes <- dataframe %>%
        filter(Comparison %in% intersect_cols,
               Sig == 'padj <= 0.05',
               log2FoldChange < 0)
    }
  
    plot_df <- sig_genes %>%
      mutate(present = 1) %>%
      pivot_wider(id_cols = c('gene_id', 'Class'),
                  names_from = Comparison,
                  values_from = present,
                  values_fill = 0)
  
  
    p <- upset( plot_df,
                   intersect = intersect_cols,
                   height_ratio = 1.2,
                   set_sizes = (
                     upset_set_size() +
                       ylab('Set Size') +
                       theme(axis.title = element_text(size = 20))
                   ),
                   base_annotations = list(
                     'Genes' = (
                       ggplot(mapping = aes(fill = Class)) +
                         geom_bar(stat = 'count', position = 'stack') +
                         scale_fill_manual(
                           values = c('Other' = 'gray60', 'IFN-responsive' = '#D6604D'),
                           name = 'Differentially\ntranscribed'
                         ) +
                         ylab('Number of genes') +
                         labs(title = title) +
                         theme(
                           plot.title = element_text(hjust = 0.5, size = 28),
                           axis.text = element_text(size = 24),
                           axis.title = element_text(size = 24),
                           legend.text = element_text(size = 20),
                           legend.title = element_text(size = 24)
                         )
                     )
                   ),
                   min_size = min_num_genes,
                   sort_sets = FALSE,
                   themes = upset_modify_themes(
                     list('default' = theme(text = element_text(size = 32)),
                          'intersections_matrix' = theme(text = element_text(size = 28),
                                                         axis.title.x = element_blank()),
                          'overall_sizes' = theme(axis.text = element_text(size = 18)))
                   )
    )
    return(p)
  }


# =========================================================================
# 6. Double check whether some gene counts correlate with library depth ===
# =========================================================================

# Theme settings  -----------------------------------------------------------

  theme_mcs <- theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )


# Get the depth of samples, by assay, condition, and person -----------------

  # Get TPM values, select condition and assay. 
  
  get_library_size <- function (df,
                                condition = c('BSA', 'IFN', 'BOTH'),
                                assay = c('PRO', 'RNA')) {

    condition <- match.arg(condition)
    assay <- match.arg(assay)

    if (condition  == 'BOTH'){
      df_temp <- df[df$assay == assay, ]

    } else {
      df_temp <- df[grepl(condition, df$samples) &
                      df$assay == assay, ]

    }

    if (nrow(df_temp) == 0) {
      stop(sprintf("Afer filtering for condition = '%s'
                   and assay = '%s', there are now rows left. ",
                   condition, assay))
    }

    split_cols <- do.call(rbind, strsplit(df_temp$samples, '_'))
    df_temp$treatment <- split_cols[, 1]
    df_temp$individual <- split_cols[, 2]
    df_temp$rep <- split_cols[, 3]

    df_avg <- aggregate(counts ~ individual, data = df_temp, FUN = mean)
    df_avg$mil <- df_avg$counts / 1000000
    df_avg$individual <- gsub(' ', '.', df_avg$individual)
    df_avg$individual <- gsub('-', '.', df_avg$individual)

    return(df_avg)
  }


# Join library depth and the df that has # of genes DE in 9 -----------------

  # Full-joins per-individual library depth (TPMs) with a per-individual 
  # gene-count df, and fill in 0 for any NA.

  join_depth_de <- function (depth, de) {
    joined <- full_join(depth, de)
    joined[is.na(joined)] <- 0

    return(joined)
  }


# Plot # of genes DE vs library depth ----------------------------------------

  # Scatter plot of library depth (TPM) vs. number of genes DE in all 
  # individuals (`df_in9`)

  plot_depth_de9 <- function (library_depth,
                              df_in9,
                              condition = c('BSA', 'IFN', 'BOTH'),
                              assay = c('PRO', 'RNA')){

    condition <- match.arg(condition)
    assay <- match.arg(assay)

    lib_depth_mean <- get_library_size(library_depth, condition, assay)

    depth_de <- join_depth_de(lib_depth_mean, df_in9)

    p <- ggplot(depth_de, aes(x = mil, y = num_genes,
                        color = individual)) +
      geom_point(size = 3) +
      scale_color_viridis_d() +
      theme_bw() +
      labs(title = 'Genes become DE in all individuals',
           subtitle = paste0(assay, ': ', condition),
           x = 'Million Mapped Reads',
           y = 'Genes DE in all individuals',
           color = 'Individual') +
      theme_mcs

    return(p)

  }
