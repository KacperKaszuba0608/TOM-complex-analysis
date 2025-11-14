library(tidyverse)
library(protti)
source("fxn.R")

# Seed for reproducibility
set.seed(123)

# Load the data
mitocarta <- vroom::vroom("./data/Human.MitoCarta3.0.csv", show_col_types = FALSE) |>
  select(UniProt, MitoCarta3.0_List, MitoCarta3.0_SubMitoLocalization) |>
  filter(UniProt != 0)

fkbp8_raw <- vroom::vroom("./data/proteinGroups_fkbp8.txt", show_col_types = FALSE) |>
  filter(is.na(`Only identified by site`),
         is.na(Reverse),
         is.na(`Potential contaminant`))

lfq <- fkbp8_raw |>
  select(`Protein IDs`, `Gene names`, starts_with("LFQ"))

lfq[lfq == 0] <- NA

metadata <- fkbp8_raw[-which(colnames(fkbp8_raw) %in% colnames(lfq)[3:14])]

# Separate for MITO and TOTAL
lfq_M <- lfq |>
  select(`Protein IDs`, `Gene names`, contains("M-"))

lfq_T <- lfq |>
  select(-contains("M-"))

# Filter out the empty rows
lfq_M <- lfq_M |>
  mutate(nas = apply(lfq_M, 1, function(row) {
    sum(is.na(row[3:8]))
  })) |>
  filter(nas < 5) |>
  select(-nas)

lfq_T <- lfq_T |>
  mutate(nas = apply(lfq_T, 1, function(row) {
    sum(is.na(row[3:8]))
    })) |>
  filter(nas < 5) |>
  select(-nas)

# Cleaning
colnames(lfq_M) <- gsub("LFQ intensity M-", "M_", colnames(lfq_M))
prots_to_rm_M <- lfq_M |>
  pivot_longer(3:8, names_to="sample", values_to="lfq") |>
  separate(sample, sep="-", c("sample", "rep")) |>
  group_by(`Protein IDs`, sample) |>
  summarise(nas = sum(is.na(lfq)), .groups = "drop") |>
  pivot_wider(id_cols = everything(), names_from = "sample", values_from = "nas") |>
  filter((M_NT == 2 & M_siFKBP8 == 1) | (M_NT == 1 & M_siFKBP8 == 2) |
           (M_NT == 2 & M_siFKBP8 == 2) | (M_NT == 3 & M_siFKBP8 > 0) |
           (M_NT > 0 & M_siFKBP8 == 3)) |>
  pull(`Protein IDs`)

lfq_M.clean <- lfq_M[-which(lfq_M$`Protein IDs` %in% prots_to_rm_M),]

colnames(lfq_T) <- gsub("LFQ intensity T-", "T_", colnames(lfq_T))
prots_to_rm_T <- lfq_T |>
  pivot_longer(3:8, names_to="sample", values_to="lfq") |>
  separate(sample, sep="-", c("sample", "rep")) |>
  group_by(`Protein IDs`, sample) |>
  summarise(nas = sum(is.na(lfq)), .groups = "drop") |>
  pivot_wider(id_cols = everything(), names_from = "sample", values_from = "nas") |>
  filter((T_NT == 2 & T_siFKBP8 == 1) | (T_NT == 1 & T_siFKBP8 == 2) |
           (T_NT == 2 & T_siFKBP8 == 2) | (T_NT == 3 & T_siFKBP8 > 0) |
           (T_NT > 0 & T_siFKBP8 == 3)) |>
  pull(`Protein IDs`)

lfq_T.clean <- lfq_T[-which(lfq_T$`Protein IDs` %in% prots_to_rm_T),]

# Imputation using MITO cells
lfq_M.miss <- lfq_M.clean |>
  mutate(prot.id = paste("prot",1:nrow(lfq_M.clean),sep="_")) |>
  pivot_longer(3:8, names_to = "Sample", values_to = "LFQvalue") |>
  separate(col=Sample, into=c("sampletype","rep"), sep = "-", remove = FALSE) |>
  mutate(sampletype_fac = as.factor(sampletype), rep = as.factor(rep),
         LFQvalue = log2(LFQvalue))

missingness <- assign_missing(
  protein.ids = lfq_M.miss$`Protein IDs`, 
  condition = lfq_M.miss$sampletype_fac,
  lfq_intensity = lfq_M.miss$LFQvalue
  )

temp <- missingness$df_wide |>
  select(prot.IDs, missingness)

lfq_M.protti <- assign_missingness(
    data = lfq_M.miss,
    sample = Sample,
    grouping = `Protein IDs`,
    condition = sampletype,
    intensity = LFQvalue
  )

lfq_M.protti <- merge(lfq_M.protti, temp, by.x = "Protein IDs", by.y = "prot.IDs",
                      suffixes = c("_protti", "_ours"))

# impute data with protti fxn using ludovic method
imputed_MITO <- impute(
  lfq_M.protti,
  sample = Sample,
  grouping = `Protein IDs`,
  intensity_log2 = LFQvalue,
  condition = sampletype,
  comparison = comparison,
  missingness = missingness_ours,
  method = "ludovic",
  skip_log2_transform_error = TRUE,
  retain_columns = missingness_protti
)

# Imputation using TOTAL cells
lfq_T.miss <- lfq_T.clean |>
  mutate(prot.id = paste("prot",1:nrow(lfq_T.clean),sep="_")) |>
  pivot_longer(3:8, names_to = "Sample", values_to = "LFQvalue") |>
  separate(col=Sample, into=c("sampletype","rep"), sep = "-", remove = FALSE) |>
  mutate(sampletype_fac = as.factor(sampletype), rep = as.factor(rep),
         LFQvalue = log2(LFQvalue))

missingness <- assign_missing(
  protein.ids = lfq_T.miss$`Protein IDs`, 
  condition = lfq_T.miss$sampletype_fac,
  lfq_intensity = lfq_T.miss$LFQvalue
)

temp <- missingness$df_wide |>
  select(prot.IDs, missingness)

lfq_T.protti <- assign_missingness(
  data = lfq_T.miss,
  sample = Sample,
  grouping = `Protein IDs`,
  condition = sampletype,
  intensity = LFQvalue
)

lfq_T.protti <- merge(lfq_T.protti, temp, by.x = "Protein IDs", by.y = "prot.IDs",
                      suffixes = c("_protti", "_ours"))

# impute data with protti fxn using ludovic method
imputed_TOTAL <- impute(
  lfq_T.protti,
  sample = Sample,
  grouping = `Protein IDs`,
  intensity_log2 = LFQvalue,
  condition = sampletype,
  comparison = comparison,
  missingness = missingness_ours,
  method = "ludovic",
  skip_log2_transform_error = TRUE,
  retain_columns = missingness_protti
)

rm(missingness, temp)

# Check density
imputed_MITO |>
  ggplot(aes(x = imputed_intensity, fill = imputed)) +
  geom_histogram()

imputed_TOTAL |>
  ggplot(aes(x = imputed_intensity, fill = imputed)) +
  geom_histogram()

# Check normality
qqnorm(imputed_TOTAL$imputed_intensity)
qqline(imputed_TOTAL$imputed_intensity)

# Fold Change MITO
colnames(imputed_MITO)[8] <- "LFQ"

means <- imputed_MITO |>
  group_by(`Protein IDs`, sampletype) |>
  summarise(FC_mean = mean(LFQ, na.rm=T), .groups = "drop") |>
  pivot_wider(id_cols = everything(), names_from = sampletype, values_from = FC_mean,
              names_prefix = "mean_") |>
  mutate(log2FC = mean_M_siFKBP8 - mean_M_NT)

imputed_MITO <- merge(imputed_MITO, means, by="Protein IDs") 

lfq_M.to_plot <- imputed_MITO |>
  select(-LFQvalue, -sampletype) |>
  pivot_wider(id_cols=everything(), names_from = Sample, values_from = c(imputed, LFQ))

lfq_M.to_plot <- lfq_M.to_plot |>
  mutate(imputed = apply(lfq_M.to_plot, 1, function(row) {
    any(row[8:13])
  }),
  p.value = apply(lfq_M.to_plot, 1, ttest,
                  grp1=grep("LFQ_M_NT", colnames(lfq_M.to_plot)),
                  grp2=grep("LFQ_M_si", colnames(lfq_M.to_plot))),
  p.adj = p.adjust(p.value, method="BH")) |>
  select(-missingness_protti)

lfq_M.to_plot <- merge(lfq_M.to_plot, metadata[c("Protein IDs", "Gene names")], by = "Protein IDs", all.x = TRUE)

# Fold change TOTAL
colnames(imputed_TOTAL)[8] <- "LFQ"

means <- imputed_TOTAL |>
  group_by(`Protein IDs`, sampletype) |>
  summarise(FC_mean = mean(LFQ, na.rm=T), .groups = "drop") |>
  pivot_wider(id_cols = everything(), names_from = sampletype, values_from = FC_mean,
              names_prefix = "mean_") |>
  mutate(log2FC = mean_T_siFKBP8 - mean_T_NT)

imputed_TOTAL <- merge(imputed_TOTAL, means, by="Protein IDs") 

lfq_T.to_plot <- imputed_TOTAL |>
  select(-LFQvalue, -sampletype) |>
  pivot_wider(id_cols=everything(), names_from = Sample, values_from = c(imputed, LFQ))

lfq_T.to_plot <- lfq_T.to_plot |>
  mutate(imputed = apply(lfq_T.to_plot, 1, function(row) {
    any(row[8:13])
  }),
  p.value = apply(lfq_T.to_plot, 1, ttest,
                  grp1=grep("LFQ_T_NT", colnames(lfq_T.to_plot)),
                  grp2=grep("LFQ_T_si", colnames(lfq_T.to_plot))),
  p.adj = p.adjust(p.value, method="BH")) |>
  select(-missingness_protti)

lfq_T.to_plot <- merge(lfq_T.to_plot, metadata[c("Protein IDs", "Gene names")], by = "Protein IDs", all.x = TRUE)

# Volcano MITO
lfq_M.to_plot <- lfq_M.to_plot |>
  mutate(UniProt = gsub(";", "|", `Protein IDs`))

to_volcano_M <- mitocarta |>
  fuzzyjoin::fuzzy_right_join(lfq_M.to_plot, by=c("UniProt" = "UniProt"),
                              match_fun = str_detect) |>
  mutate(annotate = ifelse(p.adj < 0.05 & abs(log2FC) > 1, TRUE, FALSE))

volcano_M <- to_volcano_M |>
  ggplot() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", alpha=0.4) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="black", alpha=0.4) +
  geom_point(data = subset(to_volcano_M, is.na(MitoCarta3.0_List)), 
             aes(
               x = log2FC, 
               y = -log10(p.value), 
               color = MitoCarta3.0_List,
               shape = imputed)) +
  geom_point(data = subset(to_volcano_M, !is.na(MitoCarta3.0_List)), 
             aes(
               x = log2FC, 
               y = -log10(p.value), 
               color = MitoCarta3.0_List,
               shape = imputed)) +
  ggrepel::geom_text_repel(data = subset(to_volcano_M, annotate), 
                           aes(
                             x = log2FC,
                             y = -log10(p.value),
                             color = MitoCarta3.0_List,
                             label = `Gene names`
                           ),
                           inherit.aes = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual("", values=c("MitoCarta3.0" = "darkgreen", "NA" = "darkgrey")) +
  scale_shape_manual("Imputed", values = c("TRUE"=1, "FALSE"=16))

volcano_M

# ggsave("./plots_fkbp/volcano_M.pdf", plot=volcano_M)

# Volcano TOTAL
lfq_T.to_plot <- lfq_T.to_plot |>
  mutate(UniProt = gsub(";", "|", `Protein IDs`))

to_volcano_T <- mitocarta |>
  fuzzyjoin::fuzzy_right_join(lfq_T.to_plot, by=c("UniProt" = "UniProt"),
                              match_fun = str_detect) |>
  mutate(annotate = ifelse(p.adj < 0.05 & abs(log2FC) > 1, TRUE, FALSE))

volcano_T <- to_volcano_T |>
  ggplot() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", alpha=0.4) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="black", alpha=0.4) +
  geom_point(data = subset(to_volcano_T, is.na(MitoCarta3.0_List)), 
             aes(
               x = log2FC, 
               y = -log10(p.value), 
               color = MitoCarta3.0_List,
               shape = imputed)) +
  geom_point(data = subset(to_volcano_T, !is.na(MitoCarta3.0_List)), 
             aes(
               x = log2FC, 
               y = -log10(p.value), 
               color = MitoCarta3.0_List,
               shape = imputed)) +
  ggrepel::geom_text_repel(data = subset(to_volcano_T, annotate), 
                           aes(
                             x = log2FC,
                             y = -log10(p.value),
                             color = MitoCarta3.0_List,
                             label = `Gene names`
                           ),
                           inherit.aes = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual("", values=c("MitoCarta3.0" = "darkgreen", "NA" = "darkgrey")) +
  scale_shape_manual("Imputed", values = c("TRUE"=1, "FALSE"=16))

volcano_T

# ggsave("./plots_fkbp/volcano_T.pdf", plot=volcano_T)







to_volcano_T |>
  filter(`Gene names` %in% c("FKBP8", "USP30", "BCL2L11")) |>
  select(`Gene names`, contains("LFQ"), imputed, missingness_ours) |> View()


















