suppressWarnings(source("prepare_the_data.R"))
# suppressWarnings(source(choose.files()))

# RAW DATA FRAMES
raw_dataset1 <- read_tsv("./data/raw_dataset1.tsv", show_col_types = FALSE)
raw_dataset2 <- read.table("./data/raw_dataset2.txt", sep = "\t", header = TRUE, check.names = FALSE)

# COMBINED DATA FRAME
combined_df <- merge(cleaned_data, dataset2, by.x="dataset2_id", by.y="Protein IDs", all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names)
  )

additional_yeast_homologs <- c("HSPA1A", "MTX3", "CYB5R1", "YME1L1", "QIL1", "HSPA8", 
                               "SAMM50", "RAB13", "MTX2", "MTX1", "MTCH1", "APOOL", 
                               "IMMT", "CHCHD3", "CHCHD6", "TMEM33", "ATAD1", "FAF2",
                               "USP30", "PTRH2", "RHOT1", "RHOT2", "MTCH2")

genes_to_translate <- unique(c(combined_df$Gene, combined_df$Gene_names))
translated_genes <- babelgene::orthologs(genes_to_translate, "Saccharomyces cerevisiae")
translated_genes <- translated_genes[-which(duplicated(translated_genes$ensembl)), ]

combined_df <- combined_df |>
  mutate(avg_fc = unlist(lapply(seq_along(FC_22_ev), function(i) {
    mean(c(FC_22_ev[i], Log2_enrichment_FLAG_EV[i]), na.rm = T)
  })),
  is_yeast_homolog = ifelse(Gene %in% translated_genes$human_symbol 
                            | Gene_names %in% translated_genes$human_symbol
                            | Gene %in% dataset2$Detected_yeast
                            | Gene_names %in% dataset2$Detected_yeast
                            | Gene %in% additional_yeast_homologs 
                            | Gene_names %in% additional_yeast_homologs,
                            TRUE, FALSE),
  colors_class = case_when(
    str_detect(Gene, "TOM") ~ "TOM subunits",
    is_yeast_homolog ~ "has yeast homolog",
    TRUE ~ "other"
  ),
  pvalue_dataset2 = 10^-Log10_pvalue_FLAG_EV)

combined_df$p_allmv <- unlist(sapply(1:nrow(combined_df), function(i) {
  row <- combined_df[i, c("p_22", "pvalue_dataset2")]
  ifelse(any(is.na(row)), row[which(!is.na(row))], metap::sumlog(as.numeric(row))$p)
}))

combined_df <- merge(combined_df, functional_df, by.x = "Gene", by.y = "Gene name", all.x = TRUE)

# STOCHIOMETRIC DATA FRAME
stochiometric_df <- merge(cleaned_data, mitocopies_df, by.x="mitocopies_id", by.y="Protein IDs") |>
  mutate(`Log10 mean mito-copies per cell (≥2/3 Reps)` = gsub("NaN", NA, `Log10 mean mito-copies per cell (≥2/3 Reps)`),
         `Log10 mean mito-copies per cell (≥2/3 Reps)` = as.numeric(`Log10 mean mito-copies per cell (≥2/3 Reps)`),
         UniProt = gsub("-\\d", "", Simple_ID, perl = TRUE))

stochiometric_df <- mitocarta |> filter(UniProt != "0") |>
  fuzzyjoin::fuzzy_right_join(stochiometric_df, by = c("UniProt" = "UniProt"),
                              match_fun = str_detect) |>
  mutate(MitoCarta3.0_List = ifelse(is.na(MitoCarta3.0_List), "other", "mito"),
         mito_copies_class = case_when(
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.8 & MitoCarta3.0_List == "mito" ~ "High abundant",
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.3 & MitoCarta3.0_List == "other" ~ "High abundant",
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.1 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.8 
           & MitoCarta3.0_List == "mito" ~ "Moderate abundant",
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 4.5 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.3 
           & MitoCarta3.0_List == "other" ~ "Moderate abundant",
           `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.1 & MitoCarta3.0_List == "mito" ~ "Low abundant",
           `Log10 mean mito-copies per cell (≥2/3 Reps)` < 4.5 & MitoCarta3.0_List == "other" ~ "Low abundant",
           TRUE ~ "Null value"
         ),
         is_significant = ifelse((FC_22_ev > 1 | FC_22_ev_XL > 1) & (p_22 < 0.05 & p_22_XL < 0.05),
                                 TRUE, FALSE),
         annotate = ifelse(FC_22_ev_XL >= 2.4 & is_significant & mito_copies_class != "Null value", TRUE, FALSE))

# BOOST DATA FRAME
only_FC <- only_FC |>
  mutate(FC_boost = FC_22_ev_XL - FC_22_ev,
         FC_boost_percent = round((FC_boost / FC_22_ev) * 100, 2),
         annotate = ifelse(FC_boost > 2, TRUE, FALSE)) |>
  separate(Protein.IDs, c("ID1", "Simple_ID", "Gene_name"), sep="\\|") |>
  mutate(Simple_ID = gsub("-\\d+", "", Simple_ID, perl = TRUE),
         Gene_name = gsub("_HUMAN", "", Gene_name)) |>
  select(-ID1, -Gene, -annotate)

# README
legend <- data.frame("Table_S1_Datasets" = c(
  "(A) Raw Dataset 1. \n
  Table containing the raw MS data from Msfragger.",
  "(B) Raw Dataset 2. \n
  Table containing the raw MS data from MaxQuant.",
  "(C) Related to Figures S1C, 2A, 4A, 5B. \n
  Data used to plot the figures with information from both dataset and additional databases (MitoCarta3.0, MitoCop, HCOP).",
  "(D) Related to Figure 5C. \n
  Data based on the Msfragger file containing log2 fold changes and boost of the crosslinked fold change values.",
  "(E) Related to Figure S3A. \n
  Data based on the Msfragger file containing and MitoCop mito-copies information."))

# SUPPLEMENTARY TABLE
list_of_sheets <- list("Legend" = legend,
                      "(A) Raw Dataset 1" = raw_dataset1,
                      "(B) Raw Dataset 2" = raw_dataset2,
                      "(C) Combined Dataset" = combined_df,
                      "(D) Figure 5C dataset" = only_FC,
                      "(E) Figure S3A dataset" = stochiometric_df)

openxlsx::write.xlsx(list_of_sheets, file = "supplementary_table_2.xlsx", keepNA=TRUE, na.string='NA')
