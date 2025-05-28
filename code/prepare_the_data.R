library(tidyverse)
library(ggrepel)
library(plotly)

################################# LOAD DATA ####################################
rlang::inform("Data Loading...")

data_to_plot <- read_csv('./data/data_to_plot.csv', show_col_types = FALSE) # data to plotting

only_FC <- read.csv('./data/dataset1_only_FC.csv') # data to boost plotting

mitocarta <- read.csv('./data/Human.MitoCarta3.0.csv') |> # MitoCarta3.0 dataset
  select(Symbol, MitoCarta3.0_List, MitoCarta3.0_SubMitoLocalization, MitoCarta3.0_MitoPathways, UniProt)

dataset2 <- read_csv('./data/MB_triplicate_WT_FLAG-TOMM22 vs EV.csv', show_col_types = FALSE) |>
  select(Gene_names, Detected_yeast, Log2_enrichment_FLAG_EV, Log10_pvalue_FLAG_EV,
         `Protein IDs`, `Data imputed in EV (if <2vv)`, Significant_FLAG_EV) |>
  mutate(`Data imputed in EV (if <2vv)` = ifelse(`Data imputed in EV (if <2vv)` == 'Yes', TRUE, FALSE),
         Significant_FLAG_EV = ifelse(Significant_FLAG_EV == 'Yes', TRUE, FALSE))

############################### MITOCOP DATASET ################################

mitocop_A <- readxl::read_xlsx('./data/mitocop-dataset_copy.xlsx', sheet = '(A) All protein groups')
mitocop_B <- readxl::read_xlsx('./data/mitocop-dataset_copy.xlsx', sheet = '(B) MitoCoP (1,134 genes)')

# Mito-copies per cell Dataset

mitocopies_df <- mitocop_A |>
  select(`Protein IDs`, `Gene names`, `Simplified protein IDs`, `Log10 mean mito-copies per cell (≥2/3 Reps)`)

# Functional Dataset
functional_df <- mitocop_B |>
  select(`Simplified protein IDs`, `Gene name`,
         which(grepl('Morphology, dynamics & organization', colnames(mitocop_B))):which(grepl('Unknown', colnames(mitocop_B))))

functional_df$functional_class <- apply(functional_df, 1, function(row) {
  rel_names <- names(row)[which(row == " 1")]
  ifelse(length(rel_names > 0), paste0(rel_names, collapse = '; '), NA)
})
functional_df$functional_class <- replace_na(functional_df$functional_class, 'Unknown')

functional_df <- functional_df |>
  select(`Simplified protein IDs`, `Gene name`, functional_class)

################################ DATA CLEANING #################################
rlang::inform("Data Cleaning...")

cleaned_data <- data_to_plot |>
  filter(!is.na(Protein.IDs)) |>
  distinct() |>
  mutate(annotate = ifelse(
    ((FC_22_ev_XL >= 2.5 & FC_22_ev >= 2.5)) &
      (TOMM22_XL + TOMM22_NA) > 45 | 
      Gene == 'MUL1',
    TRUE, FALSE
  ),
  sig_22 = ifelse(FC_22_ev > 2 & p_22 < 0.05, TRUE, FALSE),
  sig_22_XL = ifelse(FC_22_ev_XL > 2 & p_22_XL < 0.05, TRUE, FALSE)
  ) |>
  separate(Protein.IDs, into=c('ID1', 'Simple_ID', 'ID3'), sep = '\\|') |>
  select(-ID1, -ID3) |>
  mutate(dataset2_id = sapply(Simple_ID, function(id) {
    matching_row <- sapply(dataset2$`Protein IDs`, function(id2) id %in% unlist(strsplit(id2, ";")))

    if (any(matching_row)) {
      return(dataset2$`Protein IDs`[which(matching_row)])
    } else {
      return(id)
    }
  })) |>
  mutate(mitocopies_id = sapply(Simple_ID, function(id) {
    matching_row <- sapply(mitocopies_df$`Protein IDs`, function(id2) id %in% unlist(strsplit(id2, ";")))
    
    if (any(matching_row)) {
      return(mitocopies_df$`Protein IDs`[which(matching_row)])
    } else {
      return(id)
    }
  }))
