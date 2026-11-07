

# -------------------------------------------------------------------------
# PURPOSE -------------------------------------------------------------
# -------------------------------------------------------------------------

  # Small, generic helpers with no ties to any one analysis domain
  # - loading/cleaning functions 
  # Is sourced by `deseq2_helper_functions.R` and 
  # `timeseries_helper_functions.R`
  

# -------------------------------------------------------------------------
# Set environment -------------------------------------------------------
# -------------------------------------------------------------------------

  library(tidyverse)


# Load longest-isoform annotation ------------------------------------------

  # Reads the bed file of longest-isoform annotations and builds a 
  # Gene/Transcript/Chr key. Used after `get_longest_isoform()` to collapse 
  # Rsubread featureCounts to a single (longest) isoform per gene.

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

  # Joins the key built by `load_longest_isoform_annotation()` with
  # count data so that only the longest-isoform transcript for
  # each gene is retained going forward.

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


# Clean column names --------------------------------------------------------

  # Column names have all this sequencing info attached. Get rid of all that
  # from sample column names. 
  # The time series ('ts') and population data need different cleaning, but 
  # ultimately end up with names that conform to a standard format. 
  # IMPORTANTLY: a previous lab member gave each cell line in the population 
  # dataset a pseudonym. Even though it's not their real name, I want to get 
  # rid of that to avoid any confusion, and because country of origin is more
  # informative anyways. 
  # Drop the `DR`, `Ursula` columns if they're still around. 

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
