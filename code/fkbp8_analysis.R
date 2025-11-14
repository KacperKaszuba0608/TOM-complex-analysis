# loading libraries
library(org.Hs.eg.db) |> suppressMessages()
library(tidyverse) |> suppressMessages()
library(biomaRt) |> suppressMessages()
library(enrichplot) |> suppressMessages()
library(clusterProfiler) |> suppressMessages()

################################ FKBP8 ANALYSIS ################################

# Reading in the data
fkbp8_kd <- vroom::vroom("./data/ProteinGroups_Mayra-sep2025-revision.txt", show_col_types = FALSE) |>
  mutate(p.value_M = 10^(-`-Log10 p-value M-siFKBP8 vs M-NT`),
         p.value_T = 10^(-`-Log10 p-value T-siFKBP8 vs T-NT`))

# Extracting the necessary columns

data_mito <- fkbp8_kd |>
  dplyr::select(`Gene names`, `Protein IDs`, p.value_M, contains("M-"), -starts_with("LFQ"), -starts_with("Razor")) |>
  filter(!is.na(`Log2FC M-siFKBP8 vs M-NT`) & !is.na(`-Log10 p-value M-siFKBP8 vs M-NT`)) |>
  mutate(p.adj_M = p.adjust(p.value_M, method="fdr"), 
         sig_bool = p.adj_M < 0.05)

data_total <- fkbp8_kd |>
  dplyr::select(`Gene names`, `Protein IDs`, p.value_T, contains("T-"), -starts_with("LFQ"), -starts_with("Razor")) |>
  filter(!is.na(`Log2FC T-siFKBP8 vs T-NT`) & !is.na(`-Log10 p-value T-siFKBP8 vs T-NT`)) |>
  mutate(p.adj_T = p.adjust(p.value_T, method="fdr"), 
         sig_bool = p.adj_T < 0.05)

  # downloading all annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", version = "Ensembl Genes 115")
uniprots <- Rkeys(org.Hs.egUNIPROT)
ENSEMBL_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENSEMBL", "UNIPROT")

################################## EXPLORING ###################################

mitocarta <- vroom::vroom("./data/Human.MitoCarta3.0.csv", show_col_types = FALSE) |>
  select(UniProt, MitoCarta3.0_List) |>
  filter(UniProt != 0)

df_explor <- fkbp8_kd |>
  select(`Gene names`, `Protein IDs`, `Log2FC M-siFKBP8 vs M-NT`, `Log2FC T-siFKBP8 vs T-NT`, `MitoCarta3.0 SubMito Localization`)

df_explor$UniProt <- gsub(";", "|", df_explor$`Protein IDs`)

df_explor <- mitocarta |>
  fuzzyjoin::fuzzy_right_join(df_explor, by = c("UniProt" = "UniProt"),
                              match_fun = str_detect)

df_explor <- df_explor |>
  filter(!is.na(MitoCarta3.0_List)) |>
  mutate(fc_diff = `Log2FC M-siFKBP8 vs M-NT` - `Log2FC T-siFKBP8 vs T-NT`) |>
  dplyr::arrange(desc(abs(fc_diff)))

df_explor |>
  select(`Gene names`, `Log2FC M-siFKBP8 vs M-NT`, `Log2FC T-siFKBP8 vs T-NT`, fc_diff, `MitoCarta3.0 SubMito Localization`) |>
  write.csv("mito_diff.csv", row.names = FALSE)

################################ SUBMITO VIOLIN ################################
data_to_violin <- fkbp8_kd |>
  select(`Log2FC M-siFKBP8 vs M-NT`, `Log2FC T-siFKBP8 vs T-NT`, 
         `MitoCarta3.0 SubMito Localization`, `Gene names`) |>
  pivot_longer(cols = 1:2, names_to = "Batch", values_to = "FC") |>
  mutate(`MitoCarta3.0 SubMito Localization` = replace_na(`MitoCarta3.0 SubMito Localization`, 'unknown'),
         colors = case_when(
    `MitoCarta3.0 SubMito Localization` == "IMS" ~ "firebrick",
    `MitoCarta3.0 SubMito Localization` == "Matrix" ~ "orange",
    `MitoCarta3.0 SubMito Localization` == "Membrane" ~ "darkgreen",
    `MitoCarta3.0 SubMito Localization` == "MIM" ~ "#30D5C8",
    `MitoCarta3.0 SubMito Localization` == "MOM" ~ "skyblue",
    `MitoCarta3.0 SubMito Localization` == "unknown" ~ "#484a4d",
    TRUE ~ "grey"
  ))

mito_prot <- c("MT-CO1", "MT-ATP8", "MT-CO2", "MT-ND1")

# Calculate the outlier thresholds for each group
outlier_thresholds <- data_to_violin |>
  group_by(`MitoCarta3.0 SubMito Localization`) |>
  summarise(
    mean = mean(FC, na.rm = TRUE),
    sd = sd(FC, na.rm = TRUE),
    .groups = 'drop'
  )

# Join these thresholds back and filter main data
outlier_data <- data_to_violin |>
  left_join(outlier_thresholds, by = "MitoCarta3.0 SubMito Localization") |>
  filter((abs(FC - mean) > 2 * sd & `MitoCarta3.0 SubMito Localization` != 'unknown') |
           (`Gene names` %in% mito_prot & !is.na(FC)) |
           (grepl("TOMM", `Gene names`) & `Gene names` != "TOMM34" & !is.na(FC)))

submito_violin <- ggplot(data = data_to_violin) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "darkgrey") +
  gghalves::geom_half_violin(
    aes(
      x = reorder(`MitoCarta3.0 SubMito Localization`, FC, FUN = function(x) -median(x)),
      y = FC,
      split = Batch,
      fill = colors,
      colour = ifelse(Batch == "Log2FC M-siFKBP8 vs M-NT", colors, 'black'),
      alpha = ifelse(Batch == "Log2FC M-siFKBP8 vs M-NT", 0.5, 1)
    ),
    position = "identity", draw_quantiles = 0.5, show.legend = FALSE) +
  geom_jitter(data = filter(outlier_data, Batch != "Log2FC M-siFKBP8 vs M-NT"),
              aes(
                x = `MitoCarta3.0 SubMito Localization`,
                y = FC,
                colour = "black",
                alpha = 1), 
             position = position_nudge(x = 0.1), show.legend = FALSE) +
  geom_jitter(data = filter(outlier_data, Batch == "Log2FC M-siFKBP8 vs M-NT"),
              aes(
                x = `MitoCarta3.0 SubMito Localization`,
                y = FC,
                colour = colors,
                alpha = 0.5),
              position = position_nudge(x = -0.1), show.legend = FALSE) +
  ggrepel::geom_text_repel(data = filter(outlier_data, Batch != "Log2FC M-siFKBP8 vs M-NT"),
            aes(
              x = `MitoCarta3.0 SubMito Localization` |> as.factor() |> as.numeric() + 0.1,
              y = FC,
              label = `Gene names`,
              color = "black"), nudge_x = 0.2,
            show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0) +
  ggrepel::geom_text_repel(data = filter(outlier_data, Batch == "Log2FC M-siFKBP8 vs M-NT"),
                           aes(
                             x = `MitoCarta3.0 SubMito Localization` |> as.factor() |> as.numeric() - 0.1,
                             y = FC,
                             label = `Gene names`,
                             color = colors), nudge_x = -0.2,
                           show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0) +
  scale_alpha_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(limits = c("IMS", "Matrix", "Membrane", "MIM", "MOM", "unknown"))

submito_violin

ggsave("plots2/split_violin_fkbp8_submito.pdf", plot = submito_violin, width = 10, height = 10)

################################### BAR PLOT ###################################
# ONLY FOR MITO WITH HITS from the figure 3d

hits <- c("FKBP8", "TIM17A", "TOM22", "COA7", "TOM40", "TIM23", "SDHA", "TOM70",
          "TOM20", "MUL1", "ATP5B", "SAMM50", "VDAC1", 
          fkbp8_kd[grep("TOMM", fkbp8_kd$`Gene names`), "Gene names"] |> pull())
hits <- HGNChelper::checkGeneSymbols(hits)$Suggested.Symbol |> unique()

bar_plot <- fkbp8_kd |>
  filter(`Gene names` %in% hits & `Gene names` != "TOMM34") |>
  select(`Gene names`, `Log2FC M-siFKBP8 vs M-NT`, contains("LFQ intensity M-si")) |>
  mutate(sd = lapply(seq_along(`Gene names`), function(i) {
    v1 <- .data[["LFQ intensity M-siFKBP8-1"]][i]
    v2 <- .data[["LFQ intensity M-siFKBP8-2"]][i]
    v3 <- .data[["LFQ intensity M-siFKBP8-3"]][i]
    sd(c(v1,v2,v3))
    
  }) |> unlist()) |>
  ggplot(aes(y = reorder(`Gene names`, -`Log2FC M-siFKBP8 vs M-NT`),
             x = `Log2FC M-siFKBP8 vs M-NT`)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(xmin = `Log2FC M-siFKBP8 vs M-NT` - sd, xmax = `Log2FC M-siFKBP8 vs M-NT` + sd),
                width = 0.5) +
  theme_bw() +
  theme(axis.title.y = element_blank())

bar_plot

ggsave("plots2/bar_plot_siFKBP8.pdf", plot = bar_plot, width = 5, height = 7)

############################## GO TERM ENRICHMENT ##############################
Entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")

data_mito <- data_mito |> 
  mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) |>
  mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
data_mito <- data_mito |> left_join(Entrez_ids, by= "UNIPROT")
# 840 still un-annotated
data_mito |> filter(is.na(ENTREZID)) |> distinct(`Protein IDs`) |> dim()

data_total <- data_total |> 
  mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) |>
  mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
data_total <- data_total |> left_join(Entrez_ids, by= "UNIPROT")
# 840 still un-annotated
data_total |> filter(is.na(ENTREZID)) |> distinct(`Protein IDs`) |> dim()

################################### GO MITO ####################################
# GO enrichment up

this.group.up <- data_mito |> 
  filter(sig_bool & `Log2FC M-siFKBP8 vs M-NT` > (1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.up <- enrichGO(gene = this.group.up, OrgDb="org.Hs.eg.db", 
                         pvalueCutoff = 0.05, pAdjustMethod="none", 
                         ont=this.ont, universe = data_mito$ENTREZID)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_M <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- data_mito |> 
  filter(sig_bool & `Log2FC M-siFKBP8 vs M-NT` < (-1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           pvalueCutoff = 0.05, pAdjustMethod="none", 
                           ont=this.ont, universe = data_mito$ENTREZID)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.down_M <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_mitos_up <- go_enrichx2.up_M@result
df_go_mitos_down <- go_enrichx2.down_M@result

################################### GO TOTAL ###################################
# GO enrichment up
this.group.up <- data_total |> 
  filter(sig_bool & `Log2FC T-siFKBP8 vs T-NT` > (1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.up <- enrichGO(gene = this.group.up, OrgDb="org.Hs.eg.db", 
                         pvalueCutoff = 0.05, pAdjustMethod="none",
                         ont=this.ont, universe = data_total$ENTREZID)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_T <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- data_total |> 
  filter(sig_bool & `Log2FC T-siFKBP8 vs T-NT` < (-1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           pvalueCutoff = 0.05, pAdjustMethod="none", 
                           ont=this.ont, universe = data_total$ENTREZID)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
# go_enrichx2.down_T <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_totals_up <- go_enrichx2.up_T@result
# df_go_totals_down <- go_enrichx2.down_T@result

################################### PLOTTING ###################################

# MITO
geneList <- setNames(data_mito$`Log2FC M-siFKBP8 vs M-NT`, data_mito$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", 
                                    midpoint = 0, limits = c(-5, 5))

# Enriched Up
p_go_heat_up_M <- heatplot(go_enrichx2.up_M, foldChange = geneList) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0, vjust = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low="grey",high="green") +
  scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
  coord_flip() +
  labs(title = "p-value cutoff 0.05 & none adjustment method")

p_go_heat_up_M

ggsave("plots_GO/heatmap_up_MITO.pdf", plot = p_go_heat_up_M, width = 20, height = 7)

# Enriched Down
p_go_heat_down_M <- heatplot(go_enrichx2.down_M, foldChange = geneList) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0, vjust = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low="red",high="grey") +
  scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
  coord_flip() +
  labs(title = "p-value cutoff 0.05 & none adjustment method")

p_go_heat_down_M

ggsave("plots_GO/heatmap_down_MITO.pdf", plot = p_go_heat_down_M, width = 20, height = 7)

# TOTAL
geneList <- setNames(data_total$`Log2FC T-siFKBP8 vs T-NT`, data_total$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", 
                                    midpoint = 0, limits = c(-5, 5))
# Enriched Up
p_go_heat_up_T <- heatplot(go_enrichx2.up_T, foldChange = geneList) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0, vjust = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low="grey",high="green") +
  scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
  coord_flip() +
  labs(title = "p-value cutoff 0.05 & none adjustment method")

p_go_heat_up_T

ggsave("plots_GO/heatmap_up_TOTAL.pdf", plot = p_go_heat_up_T, width = 15, height = 7)

############################ ONE GO TERM ENRICHMENT ############################
# ONLY MITOCHONDRIUM AND OTHER ADD ALSO THE HALF VIOLIN UNDER THE VOLCANO
# COLOR THE MITOCHDONRIUM WITH DARK GREEN
# LABEL SIG PROTEIN WHICH P < 0.05 AND |FC| > 1

# Extracting necessary columns
data_one_term <- fkbp8_kd |>
  select(`Protein IDs`, `Gene names`, contains("Log"), starts_with("p."), contains("q-value", ignore.case = FALSE)) |>
  mutate(
    sig_bool_M = p.adj_M < 0.05 & abs(`Log2FC M-siFKBP8 vs M-NT`) > 1,
    sig_bool_T = p.adj_T < 0.05 & abs(`Log2FC T-siFKBP8 vs T-NT`) > 1,
    UNIPROT = gsub(";.*", "", `Protein IDs`),
    UNIPROT = gsub("-.", "", UNIPROT)
  )

# Creating a vector of needed GO terms and map vector
this.GO.mito <- "GO:0005739" # Mitochondrium
organelle <- "Mitochondrium"
names(organelle) <- this.GO.mito

# Retrieving all possible UNIPROT IDs for our GO terms
retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", 
                                   keys = this.GO.mito, 
                                   columns = c("ENSEMBL", "UNIPROT")) |>
  select(UNIPROT, GOALL) |>
  mutate(Organelles = organelle[GOALL]) |> # mapping the GO IDs with the names
  distinct(UNIPROT, Organelles) |> # Removing noise
  group_by(UNIPROT) |>
  mutate(Organelles = paste(Organelles, collapse = " | ")) |> # Collapsing more than 1 organelle
  ungroup()

# Merging with our original dataset
data_to_plot <- merge(data_one_term, retrieved, by = "UNIPROT", all.x = TRUE) |>
  mutate(Organelles = replace_na(Organelles, "Other"),
         transparency_M = ifelse(sig_bool_M, 1, 0.2),
         transparency_T = ifelse(sig_bool_T, 1, 0.2)) |>
  distinct()

# Check if NA still exist
is.na(data_to_plot$Organelles) |> sum() == 0

# MITO volcano plot
volcano_M <- ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(data = subset(data_to_plot, Organelles == "Other"), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = -log10(p.adj_M),
    color = Organelles,
    alpha = transparency_M)) +
  geom_point(data = subset(data_to_plot, Organelles != "Other"), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = -log10(p.adj_M),
    color = Organelles,
    alpha = transparency_M)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot, sig_bool_M), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = -log10(p.adj_M),
    color = Organelles,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_alpha_identity() +
  scale_color_manual("Organelles", values = c(
    'Mitochondrium' = 'darkgreen',
    'Other' = 'darkgrey'
  ), breaks = c(organelle |> unname(), 'Other'),
  labels = gsub(" ", "\n", c(organelle |> unname(), 'Other'))) +
  xlim(-3,3) +
  theme(legend.position = c(0.4, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_blank()) +
  guides(alpha = 'none')

volcano_M

density_M <- data_to_plot |>
  select(`Log2FC M-siFKBP8 vs M-NT`, Organelles) |>
  ggplot() +
  gghalves::geom_half_violin(aes(
    y = `Log2FC M-siFKBP8 vs M-NT`, 
    x = '1', 
    fill = Organelles),
    position = "identity", draw_quantiles = 0.5, side = c('r', 'l')) + 
  scale_fill_manual("", values=c("Mitochondrium" = "darkgreen",
                             "Other" = "darkgrey")) +
  guides(fill='none') +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.line = element_blank()) +
  ylim(-3,3) +
  coord_flip()

# density_M

merged_plot <- ggpubr::ggarrange(volcano_M, density_M, nrow = 2, heights = c(4,1),
                                 align = 'v')

ggsave("plots2/volcano_MITO.pdf", plot = merged_plot, width = 4, height = 6)

# TOTAL volcano plot
volcano_T <- ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(data = subset(data_to_plot, Organelles == "Other"), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = -log10(p.adj_T),
    color = Organelles,
    alpha = transparency_T)) +
  geom_point(data = subset(data_to_plot, Organelles != "Other"), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = -log10(p.adj_T),
    color = Organelles,
    alpha = transparency_T)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot, sig_bool_T), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = -log10(p.adj_T),
    color = Organelles,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_alpha_identity() +
  scale_color_manual("Organelles", values = c(
    'Mitochondrium' = 'darkgreen',
    'Other' = 'darkgrey'
  ), breaks = c(organelle |> unname(), 'Other'),
  labels = gsub(" ", "\n", c(organelle |> unname(), 'Other'))) +
  xlim(-5.5, 4) +
  theme(legend.position = c(0.4, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_blank()) +
  guides(alpha = 'none')

volcano_T
# ggsave("plots_GO/volcano_one_term_TOTAL.pdf", plot = volcano_T, height = 7)

density_T <- data_to_plot |>
  select(`Log2FC T-siFKBP8 vs T-NT`, Organelles) |>
  ggplot() +
  gghalves::geom_half_violin(aes(
    y = `Log2FC T-siFKBP8 vs T-NT`, 
    x = '', 
    fill = Organelles),
    position = "identity", draw_quantiles = 0.5, side = c('r', 'l')) + 
  scale_fill_manual("", values=c("Mitochondrium" = "darkgreen",
                                 "Other" = "darkgrey")) +
  guides(fill='none') +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.line = element_blank()) +
  coord_flip() +
  ylim(-5.5,4)

# density_T

merged_plot_T <- ggpubr::ggarrange(volcano_T, density_T, nrow = 2, heights = c(4,1),
                                 align = 'v')

ggsave("plots2/volcano_TOTAL.pdf", plot = merged_plot_T, width = 4, height = 6)

# SPLIT VIOLIN PLOT
data_to_plot_copy <- data_to_plot |>
  mutate(color = case_when(
    Organelles2 == 'Mitochondrium' ~ 'orange',
    Organelles2 == 'Endoplasmatic Reticulum' ~ 'purple',
    Organelles2 == 'Golgi aparatus' ~ 'blue',
    Organelles2 == 'Lysosome' ~ 'magenta',
    Organelles2 == 'Peroxisome' ~ 'darkgreen',
    TRUE ~ 'darkgrey'
  ))

# Pivoting the df so that the Batch and LFQs are their own columns
data_to_plot_copy <- data_to_plot_copy |> 
  tidyr::pivot_longer(c(`Log2FC M-siFKBP8 vs M-NT`, `Log2FC T-siFKBP8 vs T-NT`), 
                      names_to = "Batch", values_to = "LFQ") |>
  mutate(transparency = ifelse(Batch == 'Log2FC M-siFKBP8 vs M-NT', 'MITO', 'TOTAL'))

# create custom legend
legend_grob <- grid::grobTree(
  # Title
  grid::textGrob("Intensity\ndistribution in", x = 0.1, y = 0.9, just = "left",
           gp = grid::gpar(col = "black", fontsize = 9, fontface = "bold")),
  
  # Whole Cells color square
  grid::rectGrob(x = 0.12, y = 0.77, width = 0.05, height = 0.1,
           gp = grid::gpar(fill = "darkgrey", col = 'black')),
  grid::rectGrob(x = 0.17, y = 0.77, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "darkgreen", col = 'black')),
  grid::rectGrob(x = 0.22, y = 0.77, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "magenta", col = 'black')),
  grid::rectGrob(x = 0.27, y = 0.77, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "blue", col = 'black')),
  grid::rectGrob(x = 0.32, y = 0.77, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "purple", col = 'black')),
  grid::rectGrob(x = 0.37, y = 0.77, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "orange", col = 'black')),
  
  # Whole Cells text
  grid::textGrob("Whole Cells", x = 0.42, y = 0.77, just = "left",
           gp = grid::gpar(col = "black", fontsize = 9)),
  
  # Mitochondrial Fraction color square
  grid::rectGrob(x = 0.12, y = 0.61, width = 0.05, height = 0.1,
           gp = grid::gpar(fill = "darkgrey", col = "darkgrey", alpha=0.5)),
  grid::rectGrob(x = 0.17, y = 0.61, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "darkgreen", col = "darkgreen", alpha=0.5)),
  grid::rectGrob(x = 0.22, y = 0.61, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "magenta", col = "magenta", alpha=0.5)),
  grid::rectGrob(x = 0.27, y = 0.61, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "blue", col = "blue", alpha=0.5)),
  grid::rectGrob(x = 0.32, y = 0.61, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "purple", col = "purple", alpha=0.5)),
  grid::rectGrob(x = 0.37, y = 0.61, width = 0.05, height = 0.1, 
           gp = grid::gpar(fill = "orange", col = "orange", alpha=0.5)),
  grid::textGrob("Mitochondrial\nFraction", x = 0.42, y = 0.61, just = "left",
           gp = grid::gpar(col = "black", fontsize = 9))
)

# half violin plot
organelle_violin_split <- ggplot(data_to_plot_copy, 
       aes(x = reorder(`Organelles2`, LFQ, FUN = function(x) -median(x)), 
           y = LFQ, 
           fill = Organelles2)) +
  gghalves::geom_half_violin(
    aes(alpha = transparency, #split = Batch, 
        color = ifelse(Batch == 'Log2FC M-siFKBP8 vs M-NT', color, 'black')),
    draw_quantiles = c(0.5),
    side = c(rep('l', 6), rep('r', 6)),
    position = 'identity') + # stacking the plots into one
  ylab('Log2 FC siFKBP8 vs NT') +
  scale_color_identity() +
  scale_fill_manual(values = c(
    'Mitochondrium' = 'orange',
    'Endoplasmatic Reticulum' = 'purple',
    'Golgi aparatus' = 'blue',
    'Lysosome' = 'magenta',
    'Peroxisome' = 'darkgreen',
    'Other' = 'grey'
  )) +
  scale_alpha_manual(values = c("MITO"=0.5, "TOTAL"=1)) +
  coord_flip() +
  annotation_custom(legend_grob, xmin = 2, xmax = 6, ymin = -5, ymax = -3) +
  guides(fill = "none", alpha = 'none') +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 10)) +
  scale_x_discrete(limits = c(organelles |> unname(), "Other"), 
                   labels = gsub(" ", "\n", c(organelles |> unname(), "Other")))

organelle_violin_split

ggsave("plots_GO/half_violin_organelles_v2.pdf")

############################ Mitophagy / Apoptosis #############################

# Creating a vector of needed GO terms and map vector
this.GO.mitophagy <- "GO:0000423" # Mitophagy
this.GO.apoptosis <- "GO:0006915" # Apoptosis

GO_terms <- c(this.GO.mitophagy, this.GO.apoptosis)
bio_processes <- c('Mitophagy', 'Apoptosis')
names(bio_processes) <- GO_terms

# Retrieving all possible UNIPROT IDs for our GO terms
retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", 
                                   keys = GO_terms, 
                                   columns = c("ENSEMBL", "UNIPROT")) |>
  select(UNIPROT, GOALL) |>
  mutate(bio_processes = bio_processes[GOALL]) |> # mapping the GO IDs with the names
  distinct(UNIPROT, bio_processes) |> # Removing noise
  group_by(UNIPROT) |>
  mutate(bio_processes = paste(bio_processes, collapse = " | ")) |> # Collapsing more than 1 organelle
  ungroup()

# Merging with our original dataset
data_to_plot2 <- merge(data_one_term, retrieved, by = "UNIPROT", all.x = TRUE) |>
  mutate(bio_processes = replace_na(bio_processes, "Other"),
         annotate = ifelse(`Gene names` == 'FKBP8' | grepl('TOMM', `Gene names`), TRUE, FALSE)) |>
  distinct()

# Check if NA still exist
is.na(data_to_plot2$bio_processes) |> sum() == 0

# Density
data_to_plot2 |>
  pivot_longer(c(`Log2FC M-siFKBP8 vs M-NT`, `Log2FC T-siFKBP8 vs T-NT`), 
               names_to = "Batch", values_to = "LFQ") |>
  ggplot(aes(x = LFQ, color = Batch)) +
  geom_density()

ggsave("plots_GO/fc_density_plot.pdf")

# MITO volcano plot
volcano_M2 <- data_to_plot2 |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(data = subset(data_to_plot2, bio_processes == "Other"), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = bio_processes,
    label = `Gene names`)) +
  geom_point(data = subset(data_to_plot2, bio_processes != "Other"), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = bio_processes,
    label = `Gene names`)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot2, annotate), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = ifelse(bio_processes == "Other", 'black', bio_processes),
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_color_manual("Organelles", values = c(
    'Apoptosis' = 'orange',
    'Mitophagy' = 'purple',
    'Mitophagy | Apoptosis' = 'darkgreen',
    'Other' = 'grey',
    'black' = 'black'
  ), breaks = c("Apoptosis", "Mitophagy", "Mitophagy | Apoptosis", "Other"))

volcano_M2

ggsave("plots_GO/volcano_bio_processes_MITO_v2.pdf", plot = volcano_M2, height = 7)

# TOTAL volcano plot
volcano_T2 <- data_to_plot2 |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(data = subset(data_to_plot2, bio_processes == "Other"), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = `-Log10 p-value T-siFKBP8 vs T-NT`,
    color = bio_processes,
    label = `Gene names`)) +
  geom_point(data = subset(data_to_plot2, bio_processes != "Other"), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = `-Log10 p-value T-siFKBP8 vs T-NT`,
    color = bio_processes,
    label = `Gene names`)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot2, sig_bool_T), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = `-Log10 p-value T-siFKBP8 vs T-NT`,
    color = bio_processes,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_color_manual("Organelles", values = c(
    'Apoptosis' = 'orange',
    'Mitophagy' = 'purple',
    'Mitophagy | Apoptosis' = 'darkgreen',
    'Other' = 'grey'
  ))

volcano_T2

ggsave("plots_GO/volcano_bio_processes_TOTAL.pdf", plot = volcano_T2, height = 7)
