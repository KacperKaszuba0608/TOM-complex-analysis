# loading libraries
library(org.Hs.eg.db) |> suppressMessages()
library(tidyverse) |> suppressMessages()
library(biomaRt) |> suppressMessages()
library(enrichplot) |> suppressMessages()
library(clusterProfiler) |> suppressMessages()

# Reading in the data
fkbp8_kd <- vroom::vroom("./data/ProteinGroups_Mayra-sep2025-revision.txt", show_col_types = FALSE)

# Extracting the necessary columns

data_mito <- fkbp8_kd |>
  dplyr::select(`Gene names`, `Protein IDs`, contains("M-"), -starts_with("LFQ"), -starts_with("Razor")) |>
  mutate(sig_bool = ifelse(is.na(`Significant M-siFKBP8 vs M-NT`), FALSE, TRUE))

data_total <- fkbp8_kd |>
  dplyr::select(`Gene names`, `Protein IDs`, contains("T-"), -starts_with("LFQ"), -starts_with("Razor")) |>
  mutate(sig_bool = ifelse(is.na(`Significant T-siFKBP8 vs T-NT`), FALSE, TRUE))

# downloading all annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", version = "Ensembl Genes 115")
uniprots <- Rkeys(org.Hs.egUNIPROT)
ENSEMBL_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENSEMBL", "UNIPROT")

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
                         pvalueCutoff = 0.05, pAdjustMethod="fdr",
                         ont=this.ont)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_M <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- data_mito |> 
  filter(sig_bool & `Log2FC M-siFKBP8 vs M-NT` < (-1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           pvalueCutoff = 0.05, pAdjustMethod="fdr",
                           ont=this.ont)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
# Don't have enriched terms
# go_enrichx2.down_M <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_mitos_up <- go_enrichx2.up_M@result
# df_go_mitos_down <- go_enrichx2.down_M@result

################################### GO TOTAL ###################################
# GO enrichment up
this.group.up <- data_total |> 
  filter(sig_bool & `Log2FC T-siFKBP8 vs T-NT` > (1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.up <- enrichGO(gene = this.group.up, OrgDb="org.Hs.eg.db", 
                         pvalueCutoff = 0.05, pAdjustMethod="fdr",
                         ont=this.ont)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_T <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- data_total |> 
  filter(sig_bool & `Log2FC T-siFKBP8 vs T-NT` < (-1)) |> 
  dplyr::select(ENTREZID) |> pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           pvalueCutoff = 0.05, pAdjustMethod="fdr",
                           ont=this.ont)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
# Don't have enriched terms
# go_enrichx2.down_T <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_totals_up <- go_enrichx2.up_T@result
# df_go_totals_down <- go_enrichx2.down_T@result

################################### PLOTTING ###################################

# MITO
geneList <- setNames(data_mito$`Log2FC M-siFKBP8 vs M-NT`, data_mito$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", 
                                    midpoint = 0, limits = c(-5, 5))

p_go_heat_up_M <- heatplot(go_enrichx2.up_M, foldChange = geneList) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0, vjust = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low="grey",high="green") +
  scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
  coord_flip()

p_go_heat_up_M

ggsave("plots_GO/heatmap_up_MITO.pdf", plot = p_go_heat_up_M, width = 20, height = 7)

# TOTAL
geneList <- setNames(data_total$`Log2FC T-siFKBP8 vs T-NT`, data_total$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", 
                                    midpoint = 0, limits = c(-5, 5))

p_go_heat_up_T <- heatplot(go_enrichx2.up_T, foldChange = geneList) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low="grey",high="green") +
  scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
  coord_flip()

p_go_heat_up_T

ggsave("plots_GO/heatmap_up_TOTAL.pdf", plot = p_go_heat_up_T, width = 10, height = 5)

############################ ONE GO TERM ENRICHMENT ############################

# Extracting necessary columns
data_one_term <- fkbp8_kd |>
  select(`Protein IDs`, `Gene names`, contains("Log"), starts_with("Sig"), starts_with("q-value")) |>
  mutate(
    sig_bool_M = ifelse(is.na(`Significant M-siFKBP8 vs M-NT`), FALSE, TRUE),
    sig_bool_T = `q-value T-siFKBP8 vs T-NT` < 0.05 & abs(`Log2FC T-siFKBP8 vs T-NT`) > 1,
    UNIPROT = gsub(";.*", "", `Protein IDs`),
    UNIPROT = gsub("-.", "", UNIPROT)
  )

# Creating a vector of needed GO terms and map vector
this.GO.mito <- "GO:0005739" # Mitochondrium
this.GO.er <- "GO:0005783" # ER
this.GO.golgi <- "GO:0005794" # Golgi aparatus
this.GO.lyso <- "GO:0005764" # Lysosome
this.GO.perox <- "GO:0005777" # Peroxisome (cellular component)

GO_terms <- c(this.GO.mito, this.GO.er, this.GO.golgi, this.GO.lyso, this.GO.perox)
organelles <- c('Mitochondrium', 'Endoplasmatic Reticulum', 'Golgi aparatus', 'Lysosome', 'Peroxisome')
names(organelles) <- GO_terms

# Retrieving all possible UNIPROT IDs for our GO terms
retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", 
                                   keys = GO_terms, 
                                   columns = c("ENSEMBL", "UNIPROT")) |>
  select(UNIPROT, GOALL) |>
  mutate(Organelles = organelles[GOALL]) |> # mapping the GO IDs with the names
  distinct(UNIPROT, Organelles) |> # Removing noise
  group_by(UNIPROT) |>
  mutate(Organelles = paste(Organelles, collapse = " | ")) |> # Collapsing more than 1 organelle
  ungroup()

# Merging with our original dataset
data_to_plot <- merge(data_one_term, retrieved, by = "UNIPROT", all.x = TRUE) |>
  mutate(Organelles = replace_na(Organelles, "Other"),
         Organelles2 = gsub(" \\|.*", "", Organelles)) |>
  distinct()

# Check if NA still exist
is.na(data_to_plot$Organelles) |> sum() == 0

# MITO volcano plot
volcano_M <- data_to_plot |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = Organelles2,
    label = `Gene names`)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot, sig_bool_M), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = Organelles2,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_color_manual("Organelles", values = c(
    'Mitochondrium' = 'orange',
    'Endoplasmatic Reticulum' = 'purple',
    'Golgi aparatus' = 'blue',
    'Lysosome' = 'magenta',
    'Peroxisome' = 'darkgreen',
    'Other' = 'grey'
  ), breaks = c(organelles |> unname(), 'Other'),
  labels = gsub(" ", "\n", c(organelles |> unname(), 'Other')))

volcano_M

ggsave("plots_GO/volcano_one_term_MITO.pdf", plot = volcano_M, height = 7)

# TOTAL volcano plot
volcano_T <- data_to_plot |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = `-Log10 p-value T-siFKBP8 vs T-NT`,
    color = Organelles2,
    label = `Gene names`)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot, sig_bool_T), aes(
    x = `Log2FC T-siFKBP8 vs T-NT`,
    y = `-Log10 p-value T-siFKBP8 vs T-NT`,
    color = Organelles2,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_color_manual("Organelles", values = c(
    'Mitochondrium' = 'orange',
    'Endoplasmatic Reticulum' = 'purple',
    'Golgi aparatus' = 'blue',
    'Lysosome' = 'magenta',
    'Peroxisome' = 'darkgreen',
    'Other' = 'grey'
  ), breaks = c(organelles |> unname(), 'Other'),
  labels = gsub(" ", "\n", c(organelles |> unname(), 'Other')))

volcano_T

ggsave("plots_GO/volcano_one_term_TOTAL.pdf", plot = volcano_T, height = 7)

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
  mutate(bio_processes = replace_na(bio_processes, "Other")) |>
  distinct()

# Check if NA still exist
is.na(data_to_plot2$bio_processes) |> sum() == 0

# MITO volcano plot
volcano_M2 <- data_to_plot2 |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = bio_processes,
    label = `Gene names`)) +
  ggrepel::geom_text_repel(data = subset(data_to_plot2, sig_bool_M), aes(
    x = `Log2FC M-siFKBP8 vs M-NT`,
    y = `-Log10 p-value M-siFKBP8 vs M-NT`,
    color = bio_processes,
    label = `Gene names`
  ), show.legend = F, max.overlaps = Inf, verbose = TRUE, min.segment.length = 0) +
  scale_color_manual("Organelles", values = c(
    'Apoptosis' = 'orange',
    'Mitophagy' = 'purple',
    'Mitophagy | Apoptosis' = 'darkgreen',
    'Other' = 'grey'
  ))

volcano_M2

ggsave("plots_GO/volcano_bio_processes_MITO.pdf", plot = volcano_M2, height = 7)

# TOTAL volcano plot
volcano_T2 <- data_to_plot2 |>
  ggplot() +
  theme_classic() +
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  geom_point(aes(
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
