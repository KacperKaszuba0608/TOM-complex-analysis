message("Loading Libraries...")
library(org.Hs.eg.db) |> suppressMessages()
library(tidyverse) |> suppressMessages()
library(ggrepel) |> suppressMessages()

################################# LOAD DATA ####################################
rlang::inform(strrep("#", 30))
message("Data Loading...")

data_to_plot <- read_csv('./data/data_to_plot.csv', show_col_types = FALSE) |> # data to plotting
  mutate(protein_names = sapply(Protein.IDs, function(ID) {
    ID <- unlist(str_split(ID, "_HUMAN"))[1] |> str_split("\\|") |> unlist()
    ID[3]
  }), .after = Protein.IDs)

only_FC <- read.csv('./data/dataset1_only_FC.csv') # data to boost plotting

mitocarta <- read.csv('./data/Human.MitoCarta3.0.csv') |> # MitoCarta3.0 dataset
  dplyr::select(Symbol, MitoCarta3.0_List, MitoCarta3.0_SubMitoLocalization, MitoCarta3.0_MitoPathways, UniProt) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, "Other"))

dataset2 <- read_csv('./data/MB_triplicate_WT_FLAG-TOMM22 vs EV.csv', show_col_types = FALSE) |>
  dplyr::select(Gene_names, Detected_yeast, Log2_enrichment_FLAG_EV, Log10_pvalue_FLAG_EV,
         `Protein IDs`, `Data imputed in EV (if <2vv)`, Significant_FLAG_EV) |>
  mutate(`Data imputed in EV (if <2vv)` = ifelse(`Data imputed in EV (if <2vv)` == 'Yes', TRUE, FALSE),
         Significant_FLAG_EV = ifelse(Significant_FLAG_EV == 'Yes', TRUE, FALSE),
         Gene_names = HGNChelper::checkGeneSymbols(Gene_names)$Suggested.Symbol |> suppressMessages())

dataset2_supp <- read.table("./data/raw_dataset2.txt", sep = "\t", header = TRUE) |>
  dplyr::select(Protein.IDs, Fasta.headers) 

dataset2 <- merge(dataset2, dataset2_supp, by.x = "Protein IDs", by.y = "Protein.IDs", all.x = TRUE) |>
  mutate(protein_name = sapply(Fasta.headers, function(header) {
    header <- unlist(str_split(header, "_HUMAN"))[1] |> str_split("\\|") |> unlist()
    header[3]
  })) |> 
  dplyr::select(-Fasta.headers)

tom20 <- readxl::read_xlsx("./data/Auswertung PR80 proteingroups.xlsx", sheet = "proteinGroups") |> 
  dplyr::select(`Gene names`, `T-test Difference`, `"-Log T-test p-value"`) |>
  dplyr::mutate(`Gene names` = gsub(";.*", "", `Gene names`),
    Gene_corrected = HGNChelper::checkGeneSymbols(`Gene names`)$Suggested.Symbol |> suppressMessages(),
    Gene_corrected = gsub(" .*", "", Gene_corrected)) |>
  group_by(Gene_corrected) |>
  mutate(p.value = 10^(-`"-Log T-test p-value"`)) |>
  summarise(`T-test Difference` = mean(`T-test Difference`, na.rm = TRUE),
            p.value = ifelse(length(p.value) > 1,
                             pchisq(-2 * sum(log(p.value)), df=2*length(p.value), lower.tail = FALSE), 
                             p.value))

############################### MITOCOP DATASET ################################

mitocop_A <- readxl::read_xlsx('./data/mitocop-dataset_copy.xlsx', sheet = '(A) All protein groups') |> suppressMessages()
mitocop_B <- readxl::read_xlsx('./data/mitocop-dataset_copy.xlsx', sheet = '(B) MitoCoP (1,134 genes)') |> suppressMessages()

# Half-lives Dataset
half_lifes <- mitocop_A |>
  dplyr::select(`Protein IDs`, `Gene names`, `Protein half-lives [h]`)
colnames(half_lifes)[3] <- 'halflife'

half_lifes <- half_lifes[-which(duplicated(half_lifes$`Gene names`)),]

# Mito-copies per cell Dataset
mitocopies_df <- mitocop_A |>
  dplyr::select(`Protein IDs`, `Gene names`, `Simplified protein IDs`, `Log10 mean mito-copies per cell (≥2/3 Reps)`, 
         `Mitochondrial based on all evidence sources`) |>
  mutate(`Mitochondrial based on all evidence sources` = ifelse(is.na(`Mitochondrial based on all evidence sources`), 'other', 'mito'),
         `Log10 mean mito-copies per cell (≥2/3 Reps)` = gsub("NaN", NA, `Log10 mean mito-copies per cell (≥2/3 Reps)`),
         `Log10 mean mito-copies per cell (≥2/3 Reps)` = as.numeric(`Log10 mean mito-copies per cell (≥2/3 Reps)`))

# Functional Dataset
functional_df <- mitocop_B |>
  dplyr::select(`Simplified protein IDs`, `Gene name`,
         which(grepl('Morphology, dynamics & organization', colnames(mitocop_B))):which(grepl('Unknown', colnames(mitocop_B))))

functional_df$functional_class <- apply(functional_df, 1, function(row) {
  rel_names <- names(row)[which(row == " 1")]
  ifelse(length(rel_names > 0), paste0(rel_names, collapse = '; '), NA)
})
functional_df$functional_class <- replace_na(functional_df$functional_class, 'Unknown')

functional_df <- functional_df |>
  dplyr::select(`Simplified protein IDs`, `Gene name`, functional_class)

################################ DATA CLEANING #################################
rlang::inform(strrep("#", 30))
message("Data Cleaning...")

cleaned_data <- data_to_plot |>
  filter(!is.na(Protein.IDs)) |>
  distinct() |>
  mutate(p_22.adj = p.adjust(p_22, method = "BH"),
  p_22_XL.adj = p.adjust(p_22_XL, method = "BH"),
  p_22_22_XL.adj = p.adjust(p_22_22_XL, method = "BH"),
  sig_22 = ifelse(FC_22_ev > 2 & p_22 < 0.05, TRUE, FALSE),
  sig_22_XL = ifelse(FC_22_ev_XL > 2 & p_22_XL < 0.05, TRUE, FALSE)
  ) |>
  separate(Protein.IDs, into=c('ID1', 'Simple_ID', 'ID3'), sep = '\\|') |>
  dplyr::select(-ID1, -ID3) |>
  filter(protein_names != "FLAG")

cleaned_data <- cleaned_data[-which(duplicated(cleaned_data$protein_names)),]

############################### DATA ANNOTATING ################################
rlang::inform(strrep("#", 30))
message("Annotating the Data...")

# Annotating dataset 1
UniProt <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = cleaned_data$Gene,
                                 keytype = "SYMBOL",
                                 columns = "UNIPROT") |>
  group_by(SYMBOL) |>
  summarise(UNIPROT = paste(UNIPROT, collapse = ";")) |>
  ungroup() |> suppressMessages()

cleaned_data <- merge(cleaned_data, UniProt, by.x = "Gene", by.y = "SYMBOL", all.x = TRUE) |>
  mutate(UNIPROT = ifelse(UNIPROT == "NA", Simple_ID, UNIPROT))

# Annotating dataset 2
UniProt <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = dataset2$Gene_names,
                                 keytype = "SYMBOL",
                                 columns = "UNIPROT") |>
  group_by(SYMBOL) |>
  summarise(UNIPROT = paste(UNIPROT, collapse = ";")) |>
  ungroup() |> suppressMessages()

dataset2 <- merge(dataset2, UniProt, by.x = "Gene_names", by.y = "SYMBOL", all.x = TRUE) |>
  mutate(UNIPROT = ifelse(UNIPROT == "NA", gsub(";.*", "", `Protein IDs`), UNIPROT))

# Annotating tom20 dataset
UniProt <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = tom20$Gene_corrected,
                                 keytype = "SYMBOL",
                                 columns = "UNIPROT", fuzzy = FALSE) |>
  group_by(SYMBOL) |>
  summarise(UNIPROT = paste(UNIPROT, collapse = ";")) |>
  ungroup() |> suppressMessages()

tom20 <- merge(tom20, UniProt, by.x = "Gene_corrected", by.y = "SYMBOL", all.x = TRUE)
tom20 <- tom20[-which(is.na(tom20$Gene_corrected)),]

# Annotating mitocopies dataset
UniProt <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = mitocopies_df$`Gene names` |> unique(),
                                 keytype = "SYMBOL",
                                 columns = "UNIPROT") |>
  group_by(SYMBOL) |>
  summarise(UNIPROT = paste(UNIPROT, collapse = ";")) |>
  ungroup() |> suppressMessages()

mitocopies_df <- merge(mitocopies_df, UniProt, by.x = "Gene names", by.y = "SYMBOL", all.x = TRUE) |>
  mutate(UNIPROT = ifelse(UNIPROT == "NA", `Protein IDs`, UNIPROT)) |>
  group_by(`Gene names`, UNIPROT, .add = TRUE) |>
  mutate(`Log10 mean mito-copies per cell (≥2/3 Reps)` = mean(`Log10 mean mito-copies per cell (≥2/3 Reps)`, na.rm = TRUE)) |>
  ungroup()

mitocopies_df <- mitocopies_df[-which(duplicated(mitocopies_df$UNIPROT)),]

# Annotating functional dataset
UniProt <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = functional_df$`Gene name`,
                                 keytype = "SYMBOL",
                                 columns = "UNIPROT") |>
  group_by(SYMBOL) |>
  summarise(UNIPROT = paste(UNIPROT, collapse = ";")) |>
  ungroup() |> suppressMessages()

functional_df <- merge(functional_df, UniProt, by.x = "Gene name", by.y = "SYMBOL", all.x = TRUE) |>
  mutate(UNIPROT = ifelse(UNIPROT == "NA", `Simplified protein IDs` , UNIPROT))

rm(UniProt)

rlang::inform(strrep("#", 30))
message("The data is ready for plotting!")
gc()
