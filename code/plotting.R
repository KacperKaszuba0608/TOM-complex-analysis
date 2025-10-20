suppressWarnings(source("prepare_the_data.R"))

################################## BASE PLOTS ##################################
FC.cutoff = 1.5

main_plot <- ggplot()+
  geom_hline(yintercept=FC.cutoff, linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  geom_vline(xintercept=FC.cutoff, linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_vline(xintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  theme_classic()

######################### FIGURE S1 C ##########################

combined_df <- merge(cleaned_data, dataset2, by.x="dataset2_id", by.y="Protein IDs", all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         annotate = ifelse(abs(FC_22_ev - Log2_enrichment_FLAG_EV) > 1.5, TRUE, annotate),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names),
         protein_names = ifelse(!is.na(protein_names), protein_names, protein_name),
         pvalue_dataset2 = 10^-Log10_pvalue_FLAG_EV
  )

# Calculating confidence interval
df_to_regression <- subset(combined_df, FC_22_ev != -5 & Log2_enrichment_FLAG_EV != -5)
model <- lm(Log2_enrichment_FLAG_EV~FC_22_ev, data = df_to_regression)

combined_df <- combined_df |>
  mutate(is_significant = case_when(
    p_22.adj < 0.05 & pvalue_dataset2 < 0.05 & FC_22_ev > FC.cutoff & Log2_enrichment_FLAG_EV > FC.cutoff ~ "Both Datasets",
    (p_22.adj < 0.05 & FC_22_ev > FC.cutoff) | (pvalue_dataset2 < 0.05 & Log2_enrichment_FLAG_EV > FC.cutoff)  ~ "One Dataset",
    TRUE ~ "Not Significant"
  ),
  annotate = ifelse(FC_22_ev > FC.cutoff & Log2_enrichment_FLAG_EV > FC.cutoff & is_significant == "Both Datasets", TRUE, FALSE))

combined_plot <- main_plot +
  geom_abline(slope = 1, intercept = 0, color = "green", linewidth = 1) +
  guides(linetype = guide_legend(override.aes = list(color = "green"))) +
  geom_point(data = subset(combined_df, FC_22_ev >= 0 & Log2_enrichment_FLAG_EV >= 0), aes(
    x = FC_22_ev,
    y = Log2_enrichment_FLAG_EV,
    color = is_significant,
    label = Gene
  )) +
  scale_color_manual("Significant Class", values=c(
      "Not Significant" = "grey",
      "Both Datasets" = "darkgreen",
      "One Dataset" = "lightgreen"),
      breaks = c("Both Datasets", "One Dataset", "Not Significant")) +
  geom_text_repel(data = subset(combined_df, annotate), mapping = aes(
    x = FC_22_ev, y = Log2_enrichment_FLAG_EV, color = is_significant,
    label = protein_names),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = "Figure S1. High-confidence human TOM interactome.",
       subtitle = paste0("R\u00b2 = ", round(summary(model)$r.squared, 4)),
       x = "Dataset 1 Log2 FC TOMM22-FLAG / EV (n=3)",
       y = "Dataset 2 Log2 FC TOMM22-FLAG / EV (n=3)") +
  theme(
    legend.position = c(.97, .6),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.spacing.y = unit(0, "cm")
    ) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 1.5), 
                     labels = c("0", "2", "4", "6", "8", "10", "1.5"), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 1.5),
                     labels = c("0", "2", "4", "6", "8", "10", "1.5"), expand = c(0,0)) + 
  guides(color = guide_legend(override.aes = list(lty = NA)))

# combined_plot
# ggplotly(combined_plot)

ggsave(plot = combined_plot, filename = "plots/figure_S1C.png", width = 6, height = 6)
ggsave(plot = combined_plot, filename = "plots/figure_S1C.pdf", width = 6, height = 6)

############################# FIGURE 2 A ##############################

additional_yeast_homologs <- c("HSPA1A", "MTX3", "CYB5R1", "YME1L1", "QIL1", "HSPA8", 
                               "SAMM50", "RAB13", "MTX2", "MTX1", "MTCH1", "APOOL", 
                               "IMMT", "CHCHD3", "CHCHD6", "TMEM33", "ATAD1", "FAF2",
                               "USP30", "PTRH2", "RHOT1", "RHOT2", "MTCH2")

combined_df <- merge(cleaned_data, dataset2, by.x="dataset2_id", by.y="Protein IDs", all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names),
         protein_names = ifelse(!is.na(protein_names), protein_names, protein_name)
  )

genes_to_translate <- unique(c(combined_df$Gene, combined_df$Gene_names))
translated_genes <- babelgene::orthologs(genes_to_translate, "Saccharomyces cerevisiae")
translated_genes <- translated_genes[-which(duplicated(translated_genes$ensembl)), ]

average_df <- combined_df |>
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

average_df$p_allmv <- unlist(sapply(1:nrow(average_df), function(i) {
  row <- average_df[i, c("p_22", "pvalue_dataset2")]
  ifelse(any(is.na(row)), row[which(!is.na(row))], metap::sumlog(as.numeric(row))$p)
}))

average_df$p_allmv.adj <- p.adjust(average_df$p_allmv, method = 'BH')

average_df <- average_df |>
  mutate(annotate = ifelse(p_allmv < 0.05 & avg_fc > 2 & is_yeast_homolog, TRUE, FALSE),
         transparency = case_when(
           p_allmv.adj < 0.05 & avg_fc > 2 ~ "no_transparent",
           p_allmv < 0.05 & avg_fc > 2 ~ "less_transparent",
           TRUE ~ "high_transparent"
         ),
         annotate = ifelse(protein_names %in% c("FKBP8", "MARH5"), TRUE, annotate))

average_plot <- average_df |>
  ggplot() +
  geom_vline(xintercept=c(2), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  theme_classic() +
  geom_point(data = average_df, aes(
    x = avg_fc,
    y = -log10(p_allmv.adj),
    color = colors_class,
    alpha = ifelse(protein_names == "MARH5", "less_transparent", transparency),
    label = Gene
  )) +
  scale_color_manual("", values=c(
    "TOM subunits" = "purple",
    "has yeast homolog" = "orange",
    "other" = "darkgrey"), 
    breaks = c("TOM subunits", "has yeast homolog", "other")) +
  scale_alpha_manual("", values = c("no_transparent" = 1, 
                                    "less_transparent" = 0.6,
                                    "high_transparent" = 0.1)) +
  geom_text_repel(data = subset(average_df, annotate), mapping = aes(
    x = avg_fc, y = -log10(p_allmv.adj), color = colors_class, label = protein_names,
    alpha = ifelse(protein_names == "MARH5", "less_transparent", transparency)),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE, point.size = 0.05
  ) + 
  labs(title = "Figure 2A.",
       x = "Log2 FC TOMM22-FLAG vs EV (n=6)",
       y = "-log10 adj. p-value") +
  theme(
    legend.position = c(.2, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-4,10), breaks = c(seq(-2,10,2)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4.1), breaks = c(seq(0, 4, 1), round(-log10(0.05), 2)), 
                     labels = c(seq(0, 4, 1), round(-log10(0.05), 2)) |> as.character(), expand = c(0,0)) + 
  guides(alpha = "none")

# average_plot

ggsave(plot = average_plot, filename = "plots/figure_2A.png", width = 10, height = 5)
ggsave(plot = average_plot, filename = "plots/figure_2A.pdf", width = 10, height = 5)

############################# FIGURE S2 A ##############################

average_df <- average_df |>
  mutate(annotate = ifelse(p_allmv < 0.05 & avg_fc > 2, TRUE, FALSE))

average_plot_v2 <- average_df |>
  ggplot() +
  geom_vline(xintercept=c(2), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  theme_classic() +
  geom_point(data = average_df, aes(
    x = avg_fc,
    y = -log10(p_allmv.adj),
    color = colors_class,
    alpha = transparency,
    label = Gene
  )) +
  scale_color_manual("", values=c(
    "TOM subunits" = "purple",
    "has yeast homolog" = "orange",
    "other" = "darkgrey"), 
    breaks = c("TOM subunits", "has yeast homolog", "other")) +
  scale_alpha_manual("", values = c("no_transparent" = 1, 
                                    "less_transparent" = 0.6,
                                    "high_transparent" = 0.1)) +
  geom_text_repel(data = subset(average_df, annotate), mapping = aes(
    x = avg_fc, y = -log10(p_allmv.adj), color = colors_class, label = protein_names,
    alpha = transparency),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = "FIGURE S2A.",
       x = "Log2 FC TOMM22-FLAG vs EV (n=6)",
       y = "-log10 adj. p-value") +
  theme(
    legend.position = c(.2, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-4,10), breaks = c(seq(-2,10,2)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4.1), breaks = c(seq(0, 4, 1), round(-log10(0.05), 2)), 
                     labels = c(seq(0, 4, 1), round(-log10(0.05), 2)) |> as.character(), expand = c(0,0)) + 
  guides(alpha = "none")

# average_plot_v2
# ggplotly(average_plot_v2)

ggsave(plot = average_plot_v2, filename = "plots/figure_S2A.png", width = 10, height = 8)
ggsave(plot = average_plot_v2, filename = "plots/figure_S2A.pdf", width = 10, height = 8)

########################### FIGURE 4 B ############################

average_func_df <- merge(average_df, functional_df, by.x = "Gene", by.y = "Gene name", all.x = TRUE) |>
  mutate(annotate = ifelse(!is_yeast_homolog & p_allmv < 0.05 & avg_fc > FC.cutoff | Gene == "PTRH2", TRUE, FALSE),
         functional_class = replace_na(functional_class, "Other"),
         functional_class = case_when(
           is_yeast_homolog  ~ "Has yeast homolog",
           p_allmv < 0.05 & avg_fc > FC.cutoff ~ functional_class,
           TRUE ~ functional_class
         ))

average_func_df <- merge(average_func_df, half_lifes, by.x = "Gene", by.y = "Gene names", all.x = TRUE) |>
  mutate(halflife_cat = case_when(
    halflife < 51 ~ 'short',
    halflife > 176 ~ 'long',
    TRUE ~ 'medium'
  ),
  transparency = case_when(
    p_allmv.adj < 0.05 & avg_fc > FC.cutoff & !is_yeast_homolog ~ "no_transparent",
    p_allmv < 0.05 & avg_fc > FC.cutoff ~ "less_transparent",
    TRUE ~ "high_transparent"
  ))

average_plot_v3 <- ggplot() +
  geom_vline(xintercept=c(FC.cutoff), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_point(data = average_func_df, aes(
    x = avg_fc,
    y = -log10(p_allmv.adj),
    color = functional_class,
    size = halflife_cat,
    alpha = transparency,
    label = Gene
  )) +
  scale_alpha_manual("", values = c("no_transparent" = 1, 
                                    "less_transparent" = 0.6,
                                    "high_transparent" = 0.1)) +
  scale_size_manual("Half-life (MitoCoP)", values = c(6, 3.5, 2, 3.5),
                    breaks = c("long", "medium", "short")) +
  scale_color_manual("Functional class (MitoCoP)",
                     values = c(
                       "Amino acid metabolism" = "darkgray",
                       "Carbohydrate & energy metabolism" = "darkgray",
                       "Cardiolipin biosynthesis & remodeling" = "darkgray",
                       "Carriers, channels & transporters" = "black",
                       "Cell death/ apoptosis" = "red",
                       "Cofactor metabolism" = "darkgray",
                       "DNA-related processes" = "darkgray",
                       "Fatty acid β-oxidation" = "darkgray",
                       "Fe-S protein biogenesis" = "darkgray",
                       "Lipid metabolism" = "darkgray",
                       "Mitochondrial ribosomes" = "darkgray",
                       "Morphology, dynamics & organization" = "darkgreen",
                       "Nucleotide/ nucleoside metabolism" = "darkgray",
                       "Other" = "darkgray",
                       "Other known processes" = "darkgray",
                       "Oxidative stress response" = "purple",
                       "OXPHOS assembly & stability" = "magenta",
                       "OXPHOS complex subunits" = "darkgray",
                       "PDH complex & PDH regulation" = "darkgray",
                       "Protein import & biogenesis" = "brown",
                       "Protein maturation & folding" = "blue",
                       "Regulation & signaling" = "green",
                       "Ribosome assembly & translation" = "darkgray",
                       "RNA-related processes" = "darkgray",
                       "TCA cycle" = "darkgray",
                       "tRNA-related processes" = "darkgray",
                       "Unknown" = "darkgray",
                       "Has yeast homolog" = "orange"
                     ),
    breaks = c("Has yeast homolog", "Cell death/ apoptosis",
                "Carriers, channels & transporters", "Morphology, dynamics & organization",
                "Oxidative stress response", "OXPHOS assembly & stability",
                "Protein import & biogenesis", "Protein maturation & folding",
                "Regulation & signaling", "Other", "Unknown")) +
  theme_classic() +
  geom_text_repel(data = subset(average_func_df, annotate), mapping = aes(
    x = avg_fc, y = -log10(p_allmv.adj), color = functional_class, label = protein_names,
    alpha = transparency),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE, 
    point.size = 0.2
  ) + 
  labs(title = "Figure 4B.",
       x = "Log2 FC TOMM22-FLAG vs EV (n=6)",
       y = "-log10 adj. p-value") +
  theme(
    legend.position = c(.3,0.98),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.margin = margin(b=0)
  ) +
  scale_x_continuous(limits = c(-4,10), breaks = c(seq(-2,10,2),FC.cutoff), 
                     labels = c(seq(-2,10,2),FC.cutoff) |> as.character(), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4.1), breaks = c(seq(0, 4, 1), round(-log10(0.05), 2)), 
                     labels = c(seq(0, 4, 1), round(-log10(0.05), 2)) |> as.character(), 
                     expand = c(0,0)) + 
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 4)))

# average_plot_v3

ggsave(plot = average_plot_v3, filename = "plots/figure_4B.png", 
       width = 4*800, height = 4*467, units = "px")
ggsave(plot = average_plot_v3, filename = "plots/figure_4B.pdf", 
       width = 4*800, height = 4*467, units = "px")

##################### FIGURE 5 A v2 ######################
gene_names <- c(
  "TOMM40", "TOMM70", "TOMM5", "TOMM20", "MARCHF5", "VDAC1", "VDAC2", "FKBP8", 
  "NDKA", "NME2", "AKAP1", "MCL1", "RAB10", "KIF11", "RIF1", "BAX", "ATAD1", 
  "DARS2", "SYNJ2BP"
)

cleaned_data$UniProt <- gsub("-\\d", "", cleaned_data$Simple_ID, perl = TRUE)

xlvsnxl_DF <- mitocarta |> filter(UniProt != "0") |>
  fuzzyjoin::fuzzy_right_join(cleaned_data, by = c("UniProt" = "UniProt"),
                              match_fun = str_detect) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, "Other"),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ "IMS",
           Gene %in% "TRABD" ~ "MOM",
           MitoCarta3.0_SubMitoLocalization == "unknown" ~ "Other",
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         ),
         is_significant = ifelse((FC_22_ev > 1 | FC_22_ev_XL > 1) & (p_22 < 0.05 & p_22_XL < 0.05),
                                 TRUE, FALSE),
         annotate = ifelse(Gene %in% c("PTRH2", "TOMM70") | (FC_22_ev_XL >= 2.4 & is_significant), TRUE, annotate)#,
         )

# 2nd Version
tom20 <- readxl::read_xlsx("./data/Auswertung PR80 proteingroups.xlsx", sheet = "proteinGroups") |> 
  select(`Gene names`, `T-test Difference`, `"-Log T-test p-value"`) |>
  mutate(Gene_corrected = HGNChelper::checkGeneSymbols(`Gene names`)$Suggested.Symbol)

to_plot <- merge(xlvsnxl_DF, tom20, by.x = "Gene", by.y = "Gene_corrected", all.x = TRUE)

xlvsnxl_plot_v2 <- ggplot() +
  geom_point(data = to_plot, aes(
    x = FC_22_ev_XL,
    y = FC_22_ev,
    color = `T-test Difference`,
    shape = ifelse(is.na(`T-test Difference`), "No", "Yes"),
    label = Gene
  )) +
  geom_text_repel(data = subset(to_plot, annotate), mapping = aes(
    y = FC_22_ev, x = FC_22_ev_XL, color = `T-test Difference`,
    label = protein_names),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = "Figure",
       y = "Log2 FC (TOMM22-FLAG / empty vector) (n=3)",
       x = "Log2 FC (XL-TOMM22-FLAG / XL-empty vector) (n=3)") +
  scale_shape_manual("Was detected by Özdemir et al.?", values = c("No" = 1, "Yes" = 16),
                     breaks = c("Yes", "No")) +
  scale_color_gradientn("T-test Diffrence\nÖzdemir et al.", 
                        colors = c("darkgrey", "darkgrey", "orange", "darkgreen"),
                        values = scales::rescale(c(-2, 0, 2, 4)),
                        na.value = "darkgrey") +
  theme_classic() +
  theme(
    legend.position = c(.4, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0,8,2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(2,8,2), expand = c(0,0))

# xlvsnxl_plot_v2
# ggplotly(xlvsnxl_plot_v2)

ggsave(plot = xlvsnxl_plot_v2, filename = "plots/figure_5A_V2.png", width = 10, height = 10)
ggsave(plot = xlvsnxl_plot_v2, filename = "plots/figure_5A_V2.pdf", width = 10, height = 10)

############################### FIGURE 5 B ###############################
boost.cutoff <- 1
FC_XL.cutoff <- 2

only_FC <- only_FC |>
  mutate(FC_boost = FC_22_ev_XL - FC_22_ev,
         FC_boost_percent = round((FC_boost / FC_22_ev) * 100, 2),
         annotate = ifelse(FC_boost > 2, TRUE, FALSE)) |>
  separate(Protein.IDs, c("ID1", "Simple_ID", "Gene_name"), sep="\\|") |>
  mutate(Simple_ID = gsub("-\\d+", "", Simple_ID, perl = TRUE),
         Gene_name = gsub("_HUMAN", "", Gene_name)) |>
  select(-ID1)

only_FC_filtered <- only_FC |>
  filter(FC_boost > boost.cutoff) |>
  filter(FC_22_ev_XL > FC_XL.cutoff)

only_FC_with_mito <- merge(only_FC_filtered, mitocarta, by.x="Simple_ID", by.y="UniProt", all.x=TRUE) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, "Other"),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ "IMS",
           Gene == "TRABD" ~ "MOM",
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         ))

boosted_bar_plot <- only_FC_with_mito |>
  ggplot() +
  geom_bar(aes(
    x = FC_22_ev_XL,
    y = reorder(Gene_name, FC_boost_percent),
    fill = MitoCarta3.0_SubMitoLocalization,
  ), stat="identity", alpha=0.5) +
  geom_bar(aes(
    x = FC_22_ev,
    y = reorder(Gene_name, FC_boost_percent),
    fill = MitoCarta3.0_SubMitoLocalization
  ), stat="identity") +
  scale_fill_manual(
    "Mitochondrial Location",
    values = c(
      "Other" = "grey",
      "IMS" = "brown",
      "Matrix" = "orange",
      "Membrane" = "purple",
      "MIM" = "red",
      "MOM" = "blue"
    ),
    breaks = c("MOM", "Membrane", "IMS", "MIM", "Matrix", "Other")
  ) +
  geom_text(aes(
    x = FC_22_ev_XL,
    y = Gene_name,
    label = paste0(Gene_name, ": ", FC_boost_percent, "% boosted")
  ), hjust = -0.25) +
  xlim(0,8) +
  labs(
    title = paste0("BAR PLOT OF THE LOG2 FOLD CHANGE BOOST (FC-XL - FC) > ",boost.cutoff," and FC-XL > ", FC_XL.cutoff),
    x = "Log2 FC XL-TOMM22-FLAG / XL-EV",
    y = "Proteins boosted by cross-linking"
  ) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) +
  guides(fill = "none")

# boosted_bar_plot
# ggplotly(boosted_bar_plot, tooltip = c("y", "label"))

ggsave(plot = boosted_bar_plot, filename = "plots/figure_5B.png", width = 10, height = 15)
ggsave(plot = boosted_bar_plot, filename = "plots/figure_5B.pdf", width = 10, height = 15)

############################## FIGURE S5 B #############################

stochiometric_df <- merge(cleaned_data, mitocopies_df, by.x="mitocopies_id", by.y="Protein IDs", all.x = TRUE) |>
  mutate(`Log10 mean mito-copies per cell (≥2/3 Reps)` = gsub("NaN", NA, `Log10 mean mito-copies per cell (≥2/3 Reps)`),
         `Log10 mean mito-copies per cell (≥2/3 Reps)` = as.numeric(`Log10 mean mito-copies per cell (≥2/3 Reps)`),
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
         is_significant = ifelse(FC_22_ev_XL > 1 & p_22_XL.adj < 0.05,
                                 TRUE, FALSE),
         annotate = ifelse(FC_22_ev_XL >= 2.4 & is_significant & mito_copies_class != "Null value" 
                           | Gene %in% c("TOMM40L")
                           , TRUE, FALSE))

stochiometric_df |>
  ggplot(aes(x=`Log10 mean mito-copies per cell (≥2/3 Reps)`, colour = `Mitochondrial based on all evidence sources`)) +
  geom_density()

stochiometric_df |>
  filter(`Mitochondrial based on all evidence sources` == "mito") |>
  select(`Log10 mean mito-copies per cell (≥2/3 Reps)`) |> summary()

stochiometric_df |>
  filter(Gene == "GDAP1") |>
  select(Gene, protein_names, annotate, contains("FC"), contains("p_"))

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

stochiometric_plot <- ggplot() +
  theme_classic() +
  geom_hline(yintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  geom_vline(xintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  geom_point(data = stochiometric_df, aes(
    x = FC_22_ev_XL,
    y = FC_22_ev,
    size = mito_copies_class,
    color = MitoCarta3.0_SubMitoLocalization,
    shape = is_significant,
    label = Gene
  )) + 
  geom_text_repel(data = subset(stochiometric_df, annotate), mapping = aes(
    x = FC_22_ev_XL, y = FC_22_ev, color = MitoCarta3.0_SubMitoLocalization, label = protein_names),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  scale_shape_manual("adj. p-value (XL) < 0.05 \n& FC (XL) > 1", 
                     values = c("TRUE" = 16, "FALSE" = 1), 
                     breaks = c("TRUE", "FALSE")) +
  scale_color_manual(
    "Mitochondrial Location",
    values = c(
      "Other" = "darkgrey",
      "IMS" = "brown",
      "Matrix" = "orange",
      "Membrane" = "purple",
      "MIM" = "red",
      "MOM" = "blue"
    ),
    breaks = c("MOM", "Membrane", "IMS", "MIM", "Matrix", "Other")
  ) +
  scale_size_manual("Mito-copies class", values = c(
    "High abundant" = 4,
    "Moderate abundant" = 3,
    "Low abundant" = 2,
    "Null value" = 1
  ), breaks = c("High abundant", "Moderate abundant", "Low abundant", "Null value")) +
  labs(y = "Log2 FC TOMM22-FLAG / EV (n=3)",
       x = "Log2 FC XL-TOMM22-FLAG / XL-EV (n=3)") +
  theme(
    legend.position = c(.3, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-6, 8), breaks = seq(-6,8,2)) +
  scale_y_continuous(limits = c(-3, 8), breaks = seq(-2,8,2))

# stochiometric_plot
# ggplotly(stochiometric_plot)

ggsave(plot = stochiometric_plot, filename = "plots/figure_S5B.png", width = 10, height = 10)
ggsave(plot = stochiometric_plot, filename = "plots/figure_S5B.pdf", width = 10, height = 10)

############################## FIGURE S5 C #############################

additional_yeast_homologs <- c("HSPA1A", "MTX3", "CYB5R1", "YME1L1", "QIL1", "HSPA8", 
                               "SAMM50", "RAB13", "MTX2", "MTX1", "MTCH1", "APOOL", 
                               "IMMT", "CHCHD3", "CHCHD6", "TMEM33", "ATAD1", "FAF2",
                               "USP30", "PTRH2", "RHOT1", "RHOT2", "MTCH2")

genes_to_translate <- unique(cleaned_data$Gene)
translated_genes <- babelgene::orthologs(genes_to_translate, "Saccharomyces cerevisiae")
translated_genes <- translated_genes[-which(duplicated(translated_genes$ensembl)), ]

crosslinked_df <- cleaned_data |>
  mutate(
    is_yeast_homolog = ifelse(Gene %in% translated_genes$human_symbol 
                              | Gene %in% dataset2$Detected_yeast
                              | Gene %in% additional_yeast_homologs,
                              TRUE, FALSE),
    colors_class = case_when(
      str_detect(Gene, "TOM") ~ "TOM subunits",
      is_yeast_homolog ~ "has yeast homolog",
      TRUE ~ "other"
    ),
    annotate = ifelse(str_detect(Gene, "TOMM") 
                      | str_detect(Gene, "VDAC") 
                      | Gene %in% c("MARCHF5", "FKBP8", "PTRH2")
                      | (is_yeast_homolog & p_22_XL < 0.05 & FC_22_ev_XL > 2) 
                      , TRUE, FALSE),
    annotate = ifelse(Gene == "HM13", FALSE, annotate),
    transparency = case_when(
      p_22_XL.adj < 0.05 & FC_22_ev_XL > 2 ~ "no_transparent",
      p_22_XL < 0.05 & FC_22_ev_XL > 2 ~ "less_transparent",
      TRUE ~ "high_transparent"
    ))

crosslinked_volcano <- ggplot() +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_vline(xintercept=2, linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  theme_classic() +
  geom_point(data = crosslinked_df, aes(
    x = FC_22_ev_XL,
    y = -log10(p_22_XL.adj),
    color = colors_class,
    alpha = transparency
  )) +
  scale_color_manual("", values=c(
    "TOM subunits" = "purple",
    "has yeast homolog" = "orange",
    "other" = "darkgrey"), 
    breaks = c("TOM subunits", "has yeast homolog", "other")) +
  scale_alpha_manual("", values = c("no_transparent" = 1, 
                                    "less_transparent" = 0.6,
                                    "high_transparent" = 0.1)) +
  geom_text_repel(data = subset(crosslinked_df, annotate), 
                  aes(x = FC_22_ev_XL, y = -log10(p_22_XL.adj), color = colors_class, 
                      label = protein_names, alpha = transparency),
                  verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  labs(title = "Crosslinked Volcano",
       y = "-log10 adj. p-value",
       x = "Log2 FC XL-TOMM22-FLAG / XL-EV (n=3)") +
  guides(alpha="none") +
  theme(
    legend.position = c(.2, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-6,7), breaks = c(seq(-6,7,2)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 2), breaks = c(1, 2, round(-log10(0.05),2)),
                     labels = c(1, 2, round(-log10(0.05),2)) |> as.character(), expand = c(0,0))

# crosslinked_volcano
# ggplotly(crosslinked_volcano)

ggsave(plot = crosslinked_volcano, filename = "plots/figure_S5C.png", width = 10, height = 5)
ggsave(plot = crosslinked_volcano, filename = "plots/figure_S5C.pdf", width = 10, height = 5)

############ NEW PLOT Dataset1 VS Dataset 2 colored by Özdemir ################

combined_df <- merge(cleaned_data, dataset2, by.x="dataset2_id", by.y="Protein IDs", all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         annotate = ifelse(abs(FC_22_ev - Log2_enrichment_FLAG_EV) > 1.5, TRUE, annotate),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names),
         protein_names = ifelse(!is.na(protein_names), protein_names, protein_name),
         pvalue_dataset2 = 10^-Log10_pvalue_FLAG_EV
  )

df_d1d2 <- combined_df |>
  select(protein_names, Gene, FC_22_ev, Log2_enrichment_FLAG_EV) |>
  mutate(Gene_corrected = HGNChelper::checkGeneSymbols(Gene)$Suggested.Symbol)

tom20 <- readxl::read_xlsx("./data/Auswertung PR80 proteingroups.xlsx", sheet = "proteinGroups") |> 
  select(`Gene names`, `T-test Difference`, `"-Log T-test p-value"`) |>
  mutate(Gene_corrected = HGNChelper::checkGeneSymbols(`Gene names`)$Suggested.Symbol)

to_plot <- merge(df_d1d2, tom20, by = "Gene_corrected", all.x = TRUE) |>
  mutate(annotate = FC_22_ev > FC.cutoff & Log2_enrichment_FLAG_EV > FC.cutoff)

# Plot
combined_plot <- main_plot +
  geom_point(data = subset(to_plot, FC_22_ev >= 0 & Log2_enrichment_FLAG_EV >= 0), aes(
    x = FC_22_ev,
    y = ,
    color = `T-test Difference`,
    shape = ifelse(is.na(`T-test Difference`), 'No', 'Yes'),
    label = protein_names
  )) +
  geom_text_repel(data = subset(to_plot, annotate), mapping = aes(
    x = FC_22_ev, y = Log2_enrichment_FLAG_EV, color = `T-test Difference`,
    label = protein_names), point.size = 0.2,
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  scale_shape_manual("Was detected by\nÖzdemir et al.?", values = c("No" = 1, "Yes" = 16),
                     breaks = c("Yes", "No")) +
  scale_color_gradientn("T-test Diffrence\nÖzdemir et al.", 
                        colors = c("darkgrey", "darkgrey", "orange", "darkgreen"),
                        values = scales::rescale(c(-2, 0, 2, 4)),
                        na.value = "darkgrey") +
  labs(title = "Figure ",
       x = "Dataset 1 Log2 FC TOMM22-FLAG / EV (n=3)",
       y = "Dataset 2 Log2 FC TOMM22-FLAG / EV (n=3)") +
  theme(
    legend.position = c(.97, .77),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12, 1.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12, 1.5), expand = c(0,0))

# combined_plot

ggsave(plot = combined_plot, filename = "plots/figure_for_point_major_4.png", width = 6, height = 6)
ggsave(plot = combined_plot, filename = "plots/figure_for_point_major_4.pdf", width = 6, height = 6)


# List of hits
# tom20_sig <- tom20 |>
#   filter(`"-Log T-test p-value"` > -log10(0.05)) |>
#   filter(!(Gene_corrected %in% df_d1d2$Gene_corrected)) |>
#   select(`Gene names`, Gene_corrected, `T-test Difference`, `log2 T-test Difference`, `"-Log T-test p-value"`) |>
#   write.csv("sig_hits_of_Özdemir.csv", row.names = FALSE)

######################## Pearson Plots ########################

# Plot 1
df_to_plot1 <- cleaned_data |>
  select(protein_names, Gene, FC_22_ev, p_22, p_22.adj) |>
  mutate(Gene_corrected = HGNChelper::checkGeneSymbols(Gene)$Suggested.Symbol)

df_to_plot1 <- merge(df_to_plot1, tom20, by = "Gene_corrected", all.x = TRUE)

test1 <- cor.test(~ FC_22_ev + `T-test Difference`, data = df_to_plot1)
cor1 <- round(test1$estimate, 2)

pear_plot1 <- df_to_plot1 |>
  ggplot() +
  theme_classic() +
  geom_point(aes(y = FC_22_ev, x = `T-test Difference`)) +
  labs(subtitle = paste0("Pearson corr coef: ", cor1),
       y = "Dataset1 log2 FC TOM22 / EV (n=3)",
       x = "T-test Difference Ozdemir et al.")

pear_plot1

# Plot 2
df_to_plot2 <- dataset2 |>
  select(protein_name, Gene_names, Log2_enrichment_FLAG_EV, Log10_pvalue_FLAG_EV) |>
  mutate(Gene_corrected = HGNChelper::checkGeneSymbols(Gene_names)$Suggested.Symbol)

df_to_plot2 <- merge(df_to_plot2, tom20, by = "Gene_corrected", all.x = TRUE)

test2 <- cor.test(~ Log2_enrichment_FLAG_EV + `T-test Difference`, data = df_to_plot2)
cor2 <- round(test2$estimate, 2)

pear_plot2 <- df_to_plot2 |>
  ggplot() +
  theme_classic() +
  geom_point(aes(y = Log2_enrichment_FLAG_EV, x = `T-test Difference`)) +
  labs(subtitle = paste0("Pearson corr coef: ", cor2),
       y = "Dataset2 log2 FC TOM22 / EV (n=3)",
       x = "T-test Difference Ozdemir et al.")

pear_plot2

# Plot 3
df_to_plot3 <- combined_df |>
  filter(FC_22_ev != -5 & Log2_enrichment_FLAG_EV != -5)

test3 <- cor.test(~ FC_22_ev + Log2_enrichment_FLAG_EV, df_to_plot3)
cor3 <- round(test3$estimate, 2)

pear_plot3 <- df_to_plot3 |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = Log2_enrichment_FLAG_EV, y = FC_22_ev)) +
  labs(subtitle = paste0("Pearson corr coef: ", cor3),
       y = "Dataset1 log2 FC TOM22 / EV (n=3)",
       x = "Dataset2 log2 FC TOM22 / EV (n=3)")

pear_plot3

pear_plots <- ggpubr::ggarrange(pear_plot1, pear_plot2, pear_plot3, ncol = 3, nrow = 1)
pear_plots

ggsave("plots/pearson_plots.png", pear_plots, width = 10, height = 5)
ggsave("plots/pearson_plots.pdf", pear_plots, width = 10, height = 5)

######################## DIGGING #########################
small_tims <- c("TIM8B", "TIM13", "T10B", "TIM8A")

df_to_digg <- cleaned_data |>
  merge(mitocarta, by.x = "Simple_ID", by.y = "UniProt", all.x = TRUE) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, "Other"),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ "IMS",
           Gene == "TRABD" ~ "MOM",
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         )) |>
  select(protein_names, MitoCarta3.0_SubMitoLocalization, starts_with("FC_"), 
         starts_with("p_"), -contains(".adj"))

submito_no <- table(df_to_digg$MitoCarta3.0_SubMitoLocalization) |> as.data.frame()
submito_no_mitocarta <- table(mitocarta$MitoCarta3.0_SubMitoLocalization) |> as.data.frame()
levels(submito_no_mitocarta$Var1)[levels(submito_no_mitocarta$Var1) == ""] <- "Other"
submito_no <- merge(submito_no, submito_no_mitocarta, by = "Var1", 
                    suffixes = c("_all", "_mitocarta"), all = TRUE)
submito_no[is.na(submito_no)] <- 0

# ratio of non-crosslinked proteins
ev_22 <- df_to_digg |>
  filter(p_22 < 0.05 & FC_22_ev > 1.5) # keep only significant proteins

ev_22 <- table(ev_22$MitoCarta3.0_SubMitoLocalization) |> as.data.frame()

ev_22 <- merge(submito_no, ev_22, by = "Var1", suffixes = c("", "_22"), all.x = TRUE) |>
  mutate(Freq = replace_na(Freq, 0), # replacing NA with 0
         ratio_22 = round(Freq / Freq_all, 2),
         ratio_22_mitocarta = round(Freq / Freq_mitocarta, 2)) |> # calculating ratios
  select(Var1, ratio_22, ratio_22_mitocarta)

# ratio of crosslinked proteins
ev_22_xl <- df_to_digg |>
  filter(p_22_XL < 0.05 & FC_22_ev_XL > 1.5) # keep only significant proteins

ev_22_xl <- table(ev_22_xl$MitoCarta3.0_SubMitoLocalization) |> as.data.frame()

ev_22_xl <- merge(submito_no, ev_22_xl, by = "Var1", suffixes = c("", "_22xl"), all.x = TRUE) |>
  mutate(Freq_22xl = replace_na(Freq, 0), # replacing NA with 0
         ratio_22xl = round(Freq / Freq_all, 2),
         ratio_22xl_mitocarta = round(Freq / Freq_mitocarta, 2)) |> # calculating ratios
  select(Var1, ratio_22xl, ratio_22xl_mitocarta)

# export results
ratios <- merge(ev_22, ev_22_xl, by="Var1", all = TRUE)
colnames(ratios) <- c("MitoCarta3.0_SubMitoLocalization", "ratio 22 [%]", "ratio 22 [%] (mitocarta)", 
                      "ratio 22-xl [%]", "ratio 22-xl [%] (mitocarta)")
  
# write.csv(ratios, "submitolocalization_RATIOS.csv", row.names = FALSE)

# siginificant protein after FDR correction
tom20 <- tom20 |>
  mutate(p.value = 10^-(`"-Log T-test p-value"`))

sig_avg_df <- average_df[(which(average_df$p_allmv.adj < 0.05)), "Gene"]
sig_tom20 <- tom20[(which(tom20$p.value < 0.05)), "Gene names"] |> unlist() |> unname()
sig_tom20 <- sapply(sig_tom20, function(gene) str_split(gene, ";")) |> unlist()
sig_tom20 <- HGNChelper::checkGeneSymbols(sig_tom20)$Suggested.Symbol

# Prepare a Venn Diagram
venn <- ggvenn::ggvenn(
  data = list("Borrero, Linke et al." = sig_avg_df,
              "Özdemir et al." = sig_tom20),
  fill_color = c("orange", "grey"),
  show_percentage = FALSE,
  text_size = 8,
  fill_alpha = 0.8,
  set_name_size = 5
)

ggsave(plot = venn, filename = "plots/venn_diagram.png")
ggsave(plot = venn, filename = "plots/venn_diagram.pdf")

intersect(sig_avg_df, sig_tom20)
sig_avg_df[-which(sig_avg_df %in% sig_tom20)]
