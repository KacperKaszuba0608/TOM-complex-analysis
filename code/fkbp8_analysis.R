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

# Imputation of MNAR proteins
mnar_data <- imputed_MITO |>
  filter(missingness_ours == "MNAR")

mnar_data <- impute_data(mnar_data, col = "LFQvalue") |>
  select(`Protein IDs`, Sample, LFQvalue)

imputed_MITO <- merge(imputed_MITO, mnar_data, by=c("Protein IDs", "Sample"), all.x=TRUE)

imputed_MITO <- imputed_MITO |>
  mutate(imputed_intensity = ifelse(is.na(LFQvalue.x) & missingness_ours == "MNAR", LFQvalue.y, imputed_intensity)) |>
  select(-LFQvalue.y)

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

# Imputation of MNAR proteins
mnar_data <- imputed_TOTAL |>
  filter(missingness_ours == "MNAR")

mnar_data <- impute_data(mnar_data, col = "LFQvalue") |>
  select(`Protein IDs`, Sample, LFQvalue)

imputed_TOTAL <- merge(imputed_TOTAL, mnar_data, by=c("Protein IDs", "Sample"), all.x=TRUE) |>
  mutate(imputed_intensity = ifelse(is.na(LFQvalue.x) & missingness_ours == "MNAR", LFQvalue.y, imputed_intensity)) |>
  select(-LFQvalue.y)

rm(missingness, temp, mnar_data)

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
  select(-LFQvalue.x, -sampletype) |>
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
  select(-LFQvalue.x, -sampletype) |>
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

################################### PLOTTING ###################################
# Volcano MITO
lfq_M.to_plot <- lfq_M.to_plot |>
  mutate(UniProt = gsub(";", "|", `Protein IDs`))

to_volcano_M <- mitocarta |>
  fuzzyjoin::fuzzy_right_join(lfq_M.to_plot, by=c("UniProt" = "UniProt"),
                              match_fun = str_detect) 

to_volcano_M <- to_volcano_M |>
  mutate(annotate = ifelse(p.adj < 0.05 & abs(log2FC) > 1 | # p.adj p.value
                             `Gene names` == "FKBP8" |
                             grepl("TOMM", `Gene names`) & `Gene names` != "TOMM34", TRUE, FALSE),
         transparency = ifelse(p.adj < 0.05, 1, 0.2),
         MitoCarta3.0_List = ifelse(is.na(MitoCarta3.0_List), "Other", "Mitochondrial proteins"))

volcano_M <- to_volcano_M |>
  ggplot() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", alpha=0.4) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="black", alpha=0.4) +
  geom_point(data = subset(to_volcano_M, MitoCarta3.0_List != "Mitochondrial proteins"), 
             aes(
               x = log2FC, 
               y = -log10(p.value),
               alpha = transparency,
               color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours
               )) +
  geom_point(data = subset(to_volcano_M, MitoCarta3.0_List == "Mitochondrial proteins"), 
             aes(
               x = log2FC, 
               y = -log10(p.value),
               alpha = transparency,
               color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours,
               )) +
  ggrepel::geom_text_repel(data = subset(to_volcano_M, annotate), 
                           aes(
                             x = log2FC,
                             y = -log10(p.value),
                             color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours,
                             label = `Gene names`
                           ),
                           inherit.aes = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual("", values=c("Mitochondrial proteins" = "darkgreen", "Other" = "darkgrey"),
                     labels = c("Mitochondrial\nproteins", "Other")) +
  scale_alpha_identity() +
  labs(x = "Log2FC siFKBP8 (n=3) Dataset 6",
       y = "-Log10 p-value") +
  theme(
    legend.position = c(0.35, 1),
    legend.justification = c("right", "top"),
    legend.text.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_blank()
  ) +
  guides(alpha="none")

# volcano_M

ggsave("./plots_fkbp/volcano_M.pdf", plot=volcano_M, width = 3.6, height = 3.4)

# Volcano TOTAL
lfq_T.to_plot <- lfq_T.to_plot |>
  mutate(UniProt = gsub(";", "|", `Protein IDs`))

to_volcano_T <- mitocarta |>
  fuzzyjoin::fuzzy_right_join(lfq_T.to_plot, by=c("UniProt" = "UniProt"),
                              match_fun = str_detect)

to_volcano_T <- to_volcano_T |>
  mutate(annotate = ifelse(p.adj < 0.05 & abs(log2FC) > 1 | # p.adj p.value
                             `Gene names` == "FKBP8" |
                             grepl("TOMM", `Gene names`) & `Gene names` != "TOMM34", TRUE, FALSE),
         transparency = ifelse(p.adj < 0.05, 1, 0.2),
         MitoCarta3.0_List = ifelse(is.na(MitoCarta3.0_List), "Other", "Mitochondrial proteins"))

volcano_T <- to_volcano_T |>
  ggplot() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", alpha=0.4) +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="black", alpha=0.4) +
  geom_point(data = subset(to_volcano_T, MitoCarta3.0_List != "Mitochondrial proteins"), 
             aes(
               x = log2FC, 
               y = -log10(p.value),
               color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours,
               alpha = transparency)) +
  geom_point(data = subset(to_volcano_T, MitoCarta3.0_List == "Mitochondrial proteins"), 
             aes(
               x = log2FC, 
               y = -log10(p.value),
               color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours,
               alpha = transparency)) +
  ggrepel::geom_text_repel(data = subset(to_volcano_T, annotate), 
                           aes(
                             x = log2FC,
                             y = -log10(p.value),
                             color = MitoCarta3.0_List, #MitoCarta3.0_List, #missingness_ours,
                             label = `Gene names`
                           ),
                           inherit.aes = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual("", values=c("Mitochondrial proteins" = "darkgreen", "Other" = "darkgrey"),
                     labels = c("Mitochondrial\nproteins", "Other")) +
  scale_alpha_identity() +
  labs(x = "Log2FC siFKBP8 (n=3) Dataset t",
       y = "-Log10 p-value") +
  theme(
    legend.position = c(0.35, 1),
    legend.justification = c("right", "top"),
    legend.text.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_blank()
  ) +
  guides(alpha="none")

# volcano_T

ggsave("./plots_fkbp/volcano_T.pdf", plot=volcano_T, width = 3.6, height = 3.4)

################################ SUBMITO VIOLIN ################################
data_to_violin <- merge(lfq_M.to_plot |> select(`Protein IDs`, `Gene names`, log2FC), 
                        lfq_T.to_plot |> select(`Protein IDs`, `Gene names`, log2FC), 
                        by = c("Protein IDs", "Gene names"), all = TRUE, suffixes = c("_M", "_T")) |>
  mutate(UniProt = gsub(";", "|", `Protein IDs`)) |>
  fuzzyjoin::fuzzy_left_join(mitocarta, by = c("UniProt"="UniProt"), match_fun = str_detect) |>
  select(`Gene names`, log2FC_M, log2FC_T, MitoCarta3.0_SubMitoLocalization) |>
  pivot_longer(cols = 2:3, names_to = "Batch", values_to = "FC") 

data_to_violin <- data_to_violin |>
  mutate(MitoCarta3.0_SubMitoLocalization = case_when(
    MitoCarta3.0_SubMitoLocalization |> is.na() |
      MitoCarta3.0_SubMitoLocalization == "unknown" ~ "Other",
    TRUE ~ MitoCarta3.0_SubMitoLocalization),
         colors = case_when(
           MitoCarta3.0_SubMitoLocalization == "IMS" ~ "firebrick",
           MitoCarta3.0_SubMitoLocalization == "Matrix" ~ "orange",
           MitoCarta3.0_SubMitoLocalization == "Membrane" ~ "#8151A0",
           MitoCarta3.0_SubMitoLocalization == "MIM" ~ "red",
           MitoCarta3.0_SubMitoLocalization == "MOM" ~ "#3954A4",
           MitoCarta3.0_SubMitoLocalization == "Other" ~ "#484a4d",
           TRUE ~ "grey"
         ))

mito_prot <- c("MT-CO1", "MT-ATP8", "MT-CO2", "MT-ND1")

# Calculate the outlier thresholds for each group
outlier_thresholds <- data_to_violin |>
  group_by(MitoCarta3.0_SubMitoLocalization) |>
  summarise(
    mean = mean(FC, na.rm = TRUE),
    sd = sd(FC, na.rm = TRUE),
    .groups = 'drop'
  )

# Join these thresholds back and filter main data
outlier_data <- data_to_violin |>
  left_join(outlier_thresholds, by = "MitoCarta3.0_SubMitoLocalization") |>
  filter((abs(FC - mean) > 2 * sd & MitoCarta3.0_SubMitoLocalization != 'Other') |
           (`Gene names` %in% mito_prot & !is.na(FC)) |
           (grepl("TOMM", `Gene names`) & `Gene names` != "TOMM34" & !is.na(FC))) |>
  mutate(order = case_when(
    MitoCarta3.0_SubMitoLocalization == "Matrix" ~ 1,
    MitoCarta3.0_SubMitoLocalization == "IMS" ~ 2,
    MitoCarta3.0_SubMitoLocalization == "MIM" ~ 3,
    MitoCarta3.0_SubMitoLocalization == "Membrane" ~ 4,
    MitoCarta3.0_SubMitoLocalization == "MOM" ~ 5,
    TRUE ~ 6
  ))

submito_violin <- ggplot() +
  theme_classic() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "darkgrey") +
  annotate("rect", ymin = -1, ymax = 1, xmin = -Inf, xmax = Inf, fill = "grey", alpha = 0.5) +
  geom_text(aes(x = 6.1, y = 1.1, label = "Whole Cell Lysates", color="#484a4d"), hjust=0) +
  geom_text(aes(x = 5.9, y = 1.1, label = "Mitochondrial fractions", color="#484a4d"), alpha = 0.5, hjust=0) +
  gghalves::geom_half_violin(data = data_to_violin,
    aes(
      x = MitoCarta3.0_SubMitoLocalization,
      y = FC,
      split = Batch,
      fill = colors,
      colour = ifelse(Batch == "log2FC_M", colors, 'black'),
      alpha = ifelse(Batch == "log2FC_M", 0.5, 1)
    ),
    position = "identity", draw_quantiles = 0.5, show.legend = FALSE) +
  geom_jitter(data = filter(outlier_data, Batch != "log2FC_M"),
              aes(
                x = MitoCarta3.0_SubMitoLocalization,
                y = FC,
                colour = "black",
                alpha = 1), 
              position = position_nudge(x = 0.1), show.legend = FALSE) +
  geom_jitter(data = filter(outlier_data, Batch == "log2FC_M"),
              aes(
                x = MitoCarta3.0_SubMitoLocalization,
                y = FC,
                colour = colors,
                alpha = 0.5),
              position = position_nudge(x = -0.1), show.legend = FALSE) +
  ggrepel::geom_text_repel(data = filter(outlier_data, Batch != "log2FC_M"),
                           aes(
                             x = order + 0.1,
                             y = FC,
                             label = `Gene names`,
                             color = "black"), nudge_x = 0.17,
                           show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0) +
  ggrepel::geom_text_repel(data = filter(outlier_data, Batch == "log2FC_M"),
                           aes(
                             x = order - 0.1,
                             y = FC,
                             label = `Gene names`,
                             color = colors), nudge_x = -0.17,
                           show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0) +
  scale_alpha_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  labs(y = "log2FC siFKBP8/Ctrl (n=3), Dataset 5 (tops) and Dataset 6 (bottoms)") +
  scale_x_discrete(limits = c("Matrix", "IMS", "MIM", "Membrane", "MOM", "Other"))

# submito_violin

ggsave("plots_fkbp//split_violin_fkbp8_submito.pdf", plot = submito_violin, width = 10, height = 10)

################################### BAR PLOT ###################################
# ONLY FOR MITO WITH HITS from the figure 3d

hits <- c("FKBP8", "TIM17A", "TOM22", "COA7", "TOM40", "TIM23", "SDHA", "TOM70",
          "TOM20", "MUL1", "ATP5B", "SAMM50", "VDAC1", 
          lfq_M.to_plot[grepl("TOMM", lfq_M.to_plot$`Gene names`), "Gene names"])
hits <- HGNChelper::checkGeneSymbols(hits)$Suggested.Symbol |> unique()

bar_plot <- lfq_M.to_plot |>
  filter(`Gene names` %in% hits & `Gene names` != "TOMM34") |>
  select(`Gene names`, log2FC, contains("LFQ")) |>
  mutate(sd = lapply(seq_along(`Gene names`), function(i) {
    v1 <- .data[["LFQ_M_siFKBP8-1"]][i]
    v2 <- .data[["LFQ_M_siFKBP8-2"]][i]
    v3 <- .data[["LFQ_M_siFKBP8-3"]][i]
    sd(c(v1,v2,v3))
    
  }) |> unlist()) |>
  ggplot(aes(y = reorder(`Gene names`, -log2FC),
             x = log2FC, fill="grey", color="#484a4d")) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(xmin = log2FC - sd, xmax = log2FC + sd),
                width = 0.5, ) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  labs(x = "Log2FC siFKBP8/Ctrl (n=3, Mito), Dataset 6")

# bar_plot

ggsave("plots_fkbp/bar_plot_siFKBP8.pdf", plot = bar_plot, width = 5, height = 7)

