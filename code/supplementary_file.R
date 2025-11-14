suppressWarnings(source("prepare_the_data.R"))
# suppressWarnings(source(choose.files()))

# RAW DATA FRAMES
raw_dataset1 <- read_tsv("./data/raw_dataset1.tsv", show_col_types = FALSE)
raw_dataset2 <- read.table("./data/raw_dataset2.txt", sep = "\t", header = TRUE, check.names = FALSE)
raw_dataset3 <- data.frame()
fkbp8_dataset <- vroom::vroom("./data/proteinGroups_fkbp8.txt", show_col_types = FALSE)

# COMBINED DATA FRAME
combined_df <- merge(cleaned_data, dataset2, by = "UNIPROT", all = TRUE) |>
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

combined_df <- merge(combined_df, functional_df, by = "UNIPROT", all.x = TRUE)

# STOCHIOMETRIC DATA FRAME
stochiometric_df <- merge(cleaned_data, mitocopies_df, by = "UNIPROT", all.x = TRUE) |>
  mutate(
    mito_copies_class = case_when(
      `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.8 & `Mitochondrial based on all evidence sources` == "mito" ~ "High abundant",
      `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.2 & `Mitochondrial based on all evidence sources` == "other" ~ "High abundant",
      `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.1 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.8 
      & `Mitochondrial based on all evidence sources` == "mito" ~ "Moderate abundant",
      `Log10 mean mito-copies per cell (≥2/3 Reps)` > 4 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.2 
      & `Mitochondrial based on all evidence sources` == "other" ~ "Moderate abundant",
      `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.1 & `Mitochondrial based on all evidence sources` == "mito" ~ "Low abundant",
      `Log10 mean mito-copies per cell (≥2/3 Reps)` < 4 & `Mitochondrial based on all evidence sources` == "other" ~ "Low abundant",
      TRUE ~ "Null value"
    ),
    is_significant = ifelse(p_22_XL.adj < 0.05, TRUE, FALSE),
    annotate = ifelse(FC_22_ev_XL >= 3 & is_significant
                      | Gene %in% c("TOMM40L")
                      , TRUE, FALSE))

small_tims <- c("TIM8B", "TIM13", "T10B", "TIM8A")
stochiometric_df$UniProt <- gsub("-\\d", "", stochiometric_df$Simple_ID, perl = TRUE)

stochiometric_df <- mitocarta |> filter(UniProt != "0") |>
  fuzzyjoin::fuzzy_right_join(stochiometric_df, by = c("UniProt" = "UniProt"),
                              match_fun = str_detect) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, "Other"),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ "IMS",
           Gene %in% "TRABD" ~ "MOM",
           MitoCarta3.0_SubMitoLocalization == "unknown" ~ "Other",
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         ))

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
  "(A) Related to Figures S1c, 2a, 3b, 4a, S5a.
  Data used to plot the figures with information from both dataset and additional databases (MitoCarta3.0, MitoCop, HCOP).",
  "(B) Related to Figure 5. 
  Data based on the Msfragger file containing log2 fold changes and boost of the crosslinked fold change values.",
  "(C) Related to Figure S5b. 
  Data based on the Msfragger file containing and MitoCop mito-copies information.",
  "(D) Dataset 1,4. 
  This table contains the mitochondrial isolation data from HEK293T cells obtained from FragPipe.",
  "(E) Dataset 2. 
  This table contains the mitochondrial isolation data from HEK293T cells obtained from MaxQuant.",
  "(F) Dataset 3 
  This table contains the mitochondrial isolation data from HEK293T cells obtained from pLink.",
  "(G) Dataset from \u00d6zdemir et al. 
  This table contains the columns used from the dataset by \u00d6zdemir et al.",
  "(H) Dataset 5+6 
  Table containing the raw MS data from MaxQuant for fkbp8 KD.")
  )

# SUPPLEMENTARY TABLE
list_of_sheets <- list("Legend" = legend,
                      "(A) Dataset 1+2" = combined_df,
                      "(B) Figure 5" = only_FC,
                      "(C) Figure S5b" = stochiometric_df,
                      "(D) Dataset 1,4" = raw_dataset1,
                      "(E) Dataset 2" = raw_dataset2,
                      "(F) Dataset 3" = raw_dataset3,
                      "(G) \u00d6zdemir et al. Dataset" = tom20,
                      "(H) Dataset 5,6" = fkbp8_dataset
                      )

openxlsx::write.xlsx(list_of_sheets, file = "supplementary_table.xlsx", keepNA=TRUE, na.string='NA')
