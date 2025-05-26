suppressWarnings(source("prepare_the_data.R"))

################################## BASE PLOTS ##################################
FC.cutoff = 1.5

main_plot <- ggplot(show_legend = FALSE)+
  geom_hline(yintercept=FC.cutoff, linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  geom_vline(xintercept=FC.cutoff, linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_vline(xintercept=0, linewidth = 0.2, color = "black", alpha=0.5) +
  theme_classic()

######################### DATASET 2 VS DATASET 1 PLOT ##########################

combined_df <- merge(cleaned_data, dataset2, by.x='dataset2_id', by.y='Protein IDs', all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         annotate = ifelse(abs(FC_22_ev - Log2_enrichment_FLAG_EV) > 1.5, TRUE, annotate),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names),
         pvalue_dataset2 = 10^-Log10_pvalue_FLAG_EV
  )

# Calculating confidence interval
df_to_regression <- subset(combined_df, FC_22_ev != -5 & Log2_enrichment_FLAG_EV != -5)
model <- lm(Log2_enrichment_FLAG_EV~FC_22_ev, data = df_to_regression)

predictions <- predict(model, newdata = df_to_regression, interval = "confidence")

# Add predictions to new data frame
df_to_regression$fit <- predictions[, "fit"]
df_to_regression$lwr <- predictions[, "lwr"]
df_to_regression$upr <- predictions[, "upr"]

combined_df <- combined_df |>
  mutate(is_significant = case_when(
    p_22 < 0.05 & pvalue_dataset2 < 0.05 & FC_22_ev > FC.cutoff & Log2_enrichment_FLAG_EV > FC.cutoff ~ "Both Datasets",
    (p_22 < 0.05 & FC_22_ev > FC.cutoff) | (pvalue_dataset2 < 0.05 & Log2_enrichment_FLAG_EV > FC.cutoff)  ~ "One Dataset",
    TRUE ~ "Not Significant"
  ),
  annotate = ifelse(FC_22_ev > FC.cutoff & Log2_enrichment_FLAG_EV > FC.cutoff & is_significant == "Both Datasets", TRUE, FALSE))

combined_plot <- main_plot +
  geom_line(data = df_to_regression, 
            aes(x = FC_22_ev, y = fit, linetype = paste0("R^2 = ", round(summary(model)$r.squared, 4))),
            color = "green", size = 1) +
  geom_ribbon(data = df_to_regression,
              aes(x = FC_22_ev, ymin = lwr, ymax = upr),
              alpha = 0.3, fill = 'green') +
  scale_linetype_manual(name = "Linear Regression", values = 1) +
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
      breaks = c('Both Datasets', 'One Dataset', 'Not Significant')) +
  geom_text_repel(data = subset(combined_df, annotate), mapping = aes(
    x = FC_22_ev, y = Log2_enrichment_FLAG_EV, color = is_significant,
    label = ifelse(annotate, Gene, '')),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = 'Figure S1. High-confidence human TOM interactome.',
       x = "Dataset 1 Log2 FC (TOMM22-FLAG / empty vector) (n=3)",
       y = "Dataset 2 Log2 FC (TOMM22-FLAG / empty vector) (n=3)") +
  theme(
    legend.position = c(.97, .6),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.spacing.y = unit(0, "cm")
    ) +
  scale_x_continuous(limits = c(0, 11), breaks = c(0, 2, 4, 6, 8, 10, 1.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 11), breaks = c(0, 2, 4, 6, 8, 10, 1.5), expand = c(0,0)) + 
  guides(color = guide_legend(override.aes = list(lty = NA)))

combined_plot
# ggplotly(combined_plot)

ggsave(plot = combined_plot, filename = 'updated_plots/dataset2_vs_dataset1_plot_colored_by_significant.png', width = 6, height = 6)
ggsave(plot = combined_plot, filename = 'updated_plots/dataset2_vs_dataset1_plot_colored_by_significant.pdf', width = 6, height = 6)

############################# YEAST HOMOLOGS PLOT ##############################

additional_yeast_homologs <- c("HSPA1A", "MTX3", "CYB5R1", "YME1L1", "QIL1", "HSPA8", 
                               "SAMM50", "RAB13", "MTX2", "MTX1", "MTCH1", "APOOL", 
                               "IMMT", "CHCHD3", "CHCHD6", "TMEM33", "ATAD1", "FAF2",
                               "USP30", "PTRH2", "RHOT1", "RHOT2", "MTCH2")

combined_df <- merge(cleaned_data, dataset2, by.x='dataset2_id', by.y='Protein IDs', all = TRUE) |>
  mutate(FC_22_ev = ifelse(Log2_enrichment_FLAG_EV > 2 & is.na(FC_22_ev), -5, FC_22_ev),
         Log2_enrichment_FLAG_EV = ifelse(FC_22_ev > 2 & is.na(Log2_enrichment_FLAG_EV), -5, Log2_enrichment_FLAG_EV),
         Gene = ifelse(!is.na(Gene), Gene, Gene_names)
  )

genes_to_translate <- unique(c(combined_df$Gene, combined_df$Gene_names))
translated_genes <- babelgene::orthologs(genes_to_translate, 'Saccharomyces cerevisiae')
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
    str_detect(Gene, 'TOM') ~ 'TOM subunits',
    is_yeast_homolog ~ 'has yeast homolog',
    TRUE ~ 'other'
  ),
  pvalue_dataset2 = 10^-Log10_pvalue_FLAG_EV)

average_df$p_allmv <- unlist(sapply(1:nrow(average_df), function(i) {
  row <- average_df[i, c("p_22", "pvalue_dataset2")]
  ifelse(any(is.na(row)), row[which(!is.na(row))], metap::sumlog(as.numeric(row))$p)
}))

average_df <- average_df |>
  mutate(annotate = ifelse(str_detect(Gene, 'TOMM') 
                           | Gene == 'MARCHF5' 
                           | str_detect(Gene, 'VDAC') 
                           | Gene == 'FKBP8'
                           | Gene == "PTRH2"
                           | (is_yeast_homolog & p_allmv < 0.05 & avg_fc > FC.cutoff) 
                           , TRUE, FALSE),
         annotate = ifelse(Gene == "HM13", FALSE, annotate))

average_plot <- average_df |>
  ggplot() +
  geom_vline(xintercept=c(FC.cutoff), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_point(data = average_df, aes(
    x = avg_fc,
    y = -log10(p_allmv),
    color = colors_class,
    alpha = ifelse((p_allmv < 0.05 & avg_fc > 1.5) | Gene == "PTRH2" & !is.na(Gene), 1, 0.95),
    label = Gene
  )) +
  scale_color_manual("", values=c(
    'TOM subunits' = 'purple',
    'has yeast homolog' = 'orange',
    'other' = 'darkgrey'), 
    breaks = c('TOM subunits', 'has yeast homolog', 'other')) +
  theme_classic() +
  geom_text_repel(data = subset(average_df, annotate), mapping = aes(
    x = avg_fc, y = -log10(p_allmv), color = colors_class, label = Gene),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = 'A.',
       x = "AVG Log2 AVG FC (Dataset 1, Dataset 2) TOMM22-FLAG vs EV  (n=6)",
       y = "-log10 p-value") +
  theme(
    legend.position = c(.2, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right"
  ) +
  scale_x_continuous(limits = c(-4,10), breaks = c(seq(-2,10,2),FC.cutoff), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(0, 8, 2), round(-log10(0.05), 2)), expand = c(0,0)) + 
  guides(alpha = 'none')

average_plot
# ggplotly(average_plot)

ggsave(plot = average_plot, filename = 'Updated_plots/average_volcano_yeast.png', width = 10, height = 5)
ggsave(plot = average_plot, filename = 'updated_plots/average_volcano_yeast.pdf', width = 10, height = 5)

########################### NON-YEAST HOMOLOGS PLOT ############################

average_func_df <- merge(average_df, functional_df, by.x = 'Gene', by.y = 'Gene name', all.x = TRUE) |>
  mutate(annotate = ifelse(!is_yeast_homolog & p_allmv < 0.05 & avg_fc > FC.cutoff | Gene == 'PTRH2', TRUE, FALSE),
         functional_class = replace_na(functional_class, 'Other'),
         functional_class = case_when(
           is_yeast_homolog  ~ 'Has yeast homolog',
           p_allmv < 0.05 & avg_fc > FC.cutoff ~ functional_class,
           TRUE ~ functional_class
         ))

average_plot <- ggplot() +
  geom_vline(xintercept=c(FC.cutoff), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linewidth = 0.2, linetype="dashed", color = "black", alpha=0.5) +
  geom_point(data = average_func_df, aes(
    x = avg_fc,
    y = -log10(p_allmv),
    color = functional_class,
    alpha = ifelse(functional_class == "Has yeast homolog", 0.5, 1),
    label = Gene
  )) +
  scale_color_manual("Functional class (MitoCoP)",
                     values = c(
                       'Amino acid metabolism' = '#E45F9C',
                       'Carbohydrate & energy metabolism' = '#FF7F00',
                       'Cardiolipin biosynthesis & remodeling' = '#FFFF33',
                       'Carriers, channels & transporters' = 'black',
                       'Cell death/ apoptosis' = 'red',
                       'Cofactor metabolism' = '#A65628',
                       'DNA-related processes' = '#F781BF',
                       'Fatty acid β-oxidation' = '#999999',
                       'Fe-S protein biogenesis' = '#4DAF4A',
                       'Lipid metabolism' = '#984EA3',
                       'Mitochondrial ribosomes' = '#377EB8',
                       'Morphology, dynamics & organization' = 'darkgreen',
                       'Nucleotide/ nucleoside metabolism' = '#A6761D',
                       'Other' = 'darkgray',
                       'Other known processes' = '#66C2A5',
                       'Oxidative stress response' = 'purple',
                       'OXPHOS assembly & stability' = 'magenta',
                       'OXPHOS complex subunits' = '#80B1D3',
                       'PDH complex & PDH regulation' = '#8DD3C7',
                       'Protein import & biogenesis' = 'brown',
                       'Protein maturation & folding' = 'blue',
                       'Regulation & signaling' = 'green',
                       'Ribosome assembly & translation' = '#FB8072',
                       'RNA-related processes' = '#BEBADA',
                       'TCA cycle' = '#FDB462',
                       'tRNA-related processes' = '#B3DE69',
                       'Unknown' = 'darkgray',
                       'Has yeast homolog' = 'orange'
                     ),
    breaks = c('Has yeast homolog', 'Cell death/ apoptosis',
                'Carriers, channels & transporters', 'Morphology, dynamics & organization',
                'Oxidative stress response', 'OXPHOS assembly & stability',
                'Protein import & biogenesis', 'Protein maturation & folding',
                'Regulation & signaling', 'Other', 'Unknown')) +
  theme_classic() +
  geom_text_repel(data = subset(average_func_df, annotate), mapping = aes(
    x = avg_fc, y = -log10(p_allmv), color = functional_class, label = Gene),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = 'Figure 4. The human TOM complex interacts with components that do not have yeast homologs.\nA.',
       x = "AVG Log2 FC (Log2 FC Dataset 1 + Log2 FC Dataset 2) TOMM22-FLAG vs EV (n=3)",
       y = "-log10 p-value") +
  theme(
    legend.position = c(.3, 1.03),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-4,10), breaks = c(seq(-2,10,2),FC.cutoff), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(0, 8, 2), round(-log10(0.05), 2)), expand = c(0,0)) + 
  guides(alpha = 'none')

average_plot

ggsave(plot = average_plot, filename = 'updated_plots/average_volcano_not_yeast.png', width = 10, height = 5)
ggsave(plot = average_plot, filename = 'updated_plots/average_volcano_not_yeast.pdf', width = 10, height = 5)

##################### CORSSLINKED VS NON-CROSSLINKED PLOT ######################
small_tims <- c('TIM8B', 'TIM13', 'T10B', 'TIM8A')
gene_names <- c(
  "TOMM40", "TOMM70", "TOMM5", "TOMM20", "MARCHF5", "VDAC1", "VDAC2", "FKBP8", 
  "NDKA", "NME2", "AKAP1", "MCL1", "RAB10", "KIF11", "RIF1", "BAX", "ATAD1", 
  "DARS2", "SYNJ2BP"
)

cleaned_data$UniProt <- gsub('-\\d', '', cleaned_data$Simple_ID, perl = TRUE)

xlvsnxl_DF <- mitocarta |> filter(UniProt != "0") |>
  fuzzyjoin::fuzzy_right_join(cleaned_data, by = c('UniProt' = 'UniProt'),
                              match_fun = str_detect) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, 'Other'),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ 'IMS',
           Gene %in% 'TRABD' ~ 'MOM',
           MitoCarta3.0_SubMitoLocalization == 'unknown' ~ 'Other',
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         ),
         is_significant = ifelse((FC_22_ev > 1 | FC_22_ev_XL > 1) & (p_22 < 0.05 & p_22_XL < 0.05),
                                 TRUE, FALSE),
         annotate = ifelse(Gene %in% c('PTRH2', 'TOMM70') | (FC_22_ev_XL >= 2.4 & is_significant), TRUE, annotate)#,
         )

xlvsnxl_plot <- main_plot +
  geom_abline(intercept = c(-1.5, 0, 1.5), slope = 1, color = "grey", linetype = "dashed") +
  geom_point(data = xlvsnxl_DF, aes(
    x = FC_22_ev_XL,
    y = FC_22_ev,
    color = MitoCarta3.0_SubMitoLocalization,
    shape = is_significant,
    label = Gene
  )) +
  geom_text_repel(data = subset(xlvsnxl_DF, annotate), mapping = aes(
    y = FC_22_ev, x = FC_22_ev_XL, color = MitoCarta3.0_SubMitoLocalization,
    label = ifelse(annotate, Gene, '')),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  labs(title = 'Figure 5. Crosslink-amplified immunopurification of human TOM complex.\nB. (Log2 FC_22_ev > 0 and Log2 FC_22_ev_XL > 0)',
       y = "Log2 FC (TOMM22-FLAG / empty vector) (n=3)",
       x = "Log2 FC (XL-TOMM22-FLAG / XL-empty vector) (n=3)") +
  scale_color_manual(
    'Mitochondrial Location',
    values = c(
      'Other' = 'darkgrey',
      'IMS' = 'brown',
      'Matrix' = 'orange',
      'Membrane' = 'purple',
      'MIM' = 'red',
      'MOM' = 'blue'
    ),
    breaks = c('MOM', 'Membrane', 'IMS', 'MIM', 'Matrix', 'Other')
  ) + 
  scale_shape_manual('Both p-value < 0.05 & at least one FC > 1', 
                     values = c('TRUE' = 16, 'FALSE' = 1), breaks = c(TRUE, FALSE)) +
  theme(
    legend.position = c(.5, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-6, 8), breaks = c(seq(-6,8,2), FC.cutoff)) +
  scale_y_continuous(limits = c(-3, 8), breaks = c(seq(-2,8,2), FC.cutoff))

xlvsnxl_plot
# ggplotly(xlvsnxl_plot)

ggsave(plot = xlvsnxl_plot, filename = 'updated_plots/xl_vs_notxl_plot_with_negative.png', width = 7, height = 7)
ggsave(plot = xlvsnxl_plot, filename = 'updated_plots/xl_vs_notxl_plot_with_negative.pdf', width = 7, height = 7)

############################### BOOSTED BAR PLOT ###############################
boost.cutoff <- 1
FC_XL.cutoff <- 2

only_FC <- only_FC |>
  mutate(FC_boost = FC_22_ev_XL - FC_22_ev,
         FC_boost_percent = round((FC_boost / FC_22_ev) * 100, 2),
         annotate = ifelse(FC_boost > 2, TRUE, FALSE)) |>
  separate(Protein.IDs, c('ID1', 'Simple_ID', 'Gene_name'), sep='\\|') |>
  mutate(Simple_ID = gsub('-\\d+', '', Simple_ID, perl = TRUE),
         Gene_name = gsub('_HUMAN', '', Gene_name)) |>
  select(-ID1)

only_FC_filtered <- only_FC |>
  filter(FC_boost > boost.cutoff) |>
  filter(FC_22_ev_XL > FC_XL.cutoff)

only_FC_with_mito <- merge(only_FC_filtered, mitocarta, by.x='Simple_ID', by.y='UniProt', all.x=TRUE) |>
  mutate(MitoCarta3.0_SubMitoLocalization = replace_na(MitoCarta3.0_SubMitoLocalization, 'Other'),
         MitoCarta3.0_SubMitoLocalization = case_when(
           Gene %in% small_tims ~ 'IMS',
           Gene == 'TRABD' ~ 'MOM',
           TRUE ~ MitoCarta3.0_SubMitoLocalization
         ))

boosted_bar_plot <- only_FC_with_mito |>
  ggplot() +
  geom_bar(aes(
    x = FC_22_ev_XL,
    y = reorder(Gene_name, FC_boost_percent),
    fill = MitoCarta3.0_SubMitoLocalization,
  ), stat='identity', alpha=0.5) +
  geom_bar(aes(
    x = FC_22_ev,
    y = reorder(Gene_name, FC_boost_percent),
    fill = MitoCarta3.0_SubMitoLocalization
  ), stat='identity') +
  scale_fill_manual(
    'Mitochondrial Location',
    values = c(
      'Other' = 'grey',
      'IMS' = 'brown',
      'Matrix' = 'orange',
      'Membrane' = 'purple',
      'MIM' = 'red',
      'MOM' = 'blue'
    ),
    breaks = c('MOM', 'Membrane', 'IMS', 'MIM', 'Matrix', 'Other')
  ) +
  geom_text(aes(
    x = FC_22_ev_XL,
    y = Gene_name,
    label = paste0(Gene, ': ', FC_boost_percent, '% boosted')
  ), hjust = -0.25) +
  xlim(0,8) +
  labs(
    title = paste0('BAR PLOT OF THE LOG2 FOLD CHANGE BOOST (FC-XL - FC) > ',boost.cutoff,' and FC-XL > ', FC_XL.cutoff),
    x = 'Log2 FC (XL-TOMM22-FLAG / XL-empty vector)',
    y = 'Proteins, sorted by boost from crosslinking'
  ) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

boosted_bar_plot
# ggplotly(boosted_bar_plot, tooltip = c('y', 'label'))

ggsave(plot = boosted_bar_plot, filename = 'plots/boosted_bar_plot.png', width = 15, height = 10)
ggsave(plot = boosted_bar_plot, filename = 'plots/boosted_bar_plot.pdf', width = 15, height = 10)

############################## STOCHIOMETRICS PLOT #############################

stochiometric_df <- merge(cleaned_data, mitocopies_df, by.x='mitocopies_id', by.y='Protein IDs') |>
  mutate(`Log10 mean mito-copies per cell (≥2/3 Reps)` = gsub('NaN', NA, `Log10 mean mito-copies per cell (≥2/3 Reps)`),
         `Log10 mean mito-copies per cell (≥2/3 Reps)` = as.numeric(`Log10 mean mito-copies per cell (≥2/3 Reps)`),
         UniProt = gsub('-\\d', '', Simple_ID, perl = TRUE))

stochiometric_df <- mitocarta |> filter(UniProt != "0") |>
  fuzzyjoin::fuzzy_right_join(stochiometric_df, by = c('UniProt' = 'UniProt'),
                              match_fun = str_detect) |>
  mutate(MitoCarta3.0_List = ifelse(is.na(MitoCarta3.0_List), 'other', 'mito'),
         mito_copies_class = case_when(
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.8 & MitoCarta3.0_List == 'mito' ~ 'High abundant',
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.3 & MitoCarta3.0_List == 'other' ~ 'High abundant',
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 5.1 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.8 
           & MitoCarta3.0_List == 'mito' ~ 'Moderate abundant',
           `Log10 mean mito-copies per cell (≥2/3 Reps)` > 4.5 & `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.3 
           & MitoCarta3.0_List == 'other' ~ 'Moderate abundant',
           `Log10 mean mito-copies per cell (≥2/3 Reps)` < 5.1 & MitoCarta3.0_List == 'mito' ~ 'Low abundant',
           `Log10 mean mito-copies per cell (≥2/3 Reps)` < 4.5 & MitoCarta3.0_List == 'other' ~ 'Low abundant',
           TRUE ~ 'Null value'
         ),
         is_significant = ifelse((FC_22_ev > 1 | FC_22_ev_XL > 1) & (p_22 < 0.05 & p_22_XL < 0.05),
                                 TRUE, FALSE),
         annotate = ifelse(FC_22_ev_XL >= 2.4 & is_significant & mito_copies_class != "Null value", TRUE, FALSE))

stochiometric_df |>
  ggplot(aes(x=`Log10 mean mito-copies per cell (≥2/3 Reps)`, colour = MitoCarta3.0_List)) +
  geom_density()

stochiometric_df |>
  filter(MitoCarta3.0_List == 'other') |>
  select(`Log10 mean mito-copies per cell (≥2/3 Reps)`) |> summary()

stochiometric_plot <- main_plot +
  geom_abline(intercept = c(-1.5, 0, 1.5), slope = 1, color = "grey", linetype = "dashed") +
  geom_point(data = stochiometric_df, aes(
    x = FC_22_ev_XL,
    y = FC_22_ev,
    color = mito_copies_class,
    shape = MitoCarta3.0_List,
    label = Gene
  )) + 
  geom_text_repel(data = subset(stochiometric_df, annotate), mapping = aes(
    x = FC_22_ev_XL, y = FC_22_ev, color = mito_copies_class,
    label = ifelse(annotate, Gene, '')),
    verbose = TRUE, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
  ) + 
  # scale_shape_manual('Both p-value < 0.05 & at least one FC > 1', 
  #                    values = c('TRUE' = 16, 'FALSE' = 1), breaks = c(TRUE, FALSE)) +
  scale_color_manual("Mito-copies class", values = c(
    "High abundant" = 'orange',
    "Moderate abundant" = "darkgreen",
    "Low abundant" = "blue",
    "Null value" = "grey"
  ), breaks = c("High abundant", "Moderate abundant", "Low abundant", "Null value")) +
  labs(title = 'FC vs FC colored by stochiometrics (Log10 mean mito-copies per cell)',
       y = "Log2 FC (TOMM22-FLAG / empty vector) (n=3)",
       x = "Log2 FC (XL-TOMM22-FLAG / XL-empty vector) (n=3)") +
  theme(
    legend.position = c(.5, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(-6, 8), breaks = c(seq(-6,8,2), FC.cutoff)) +
  scale_y_continuous(limits = c(-3, 8), breaks = c(seq(-2,8,2), FC.cutoff))

stochiometric_plot
# ggplotly(stochiometric_plot)

ggsave(plot = stochiometric_plot, filename = 'updated_plots/xl_vs_nxl_with_mitocopies.png', width = 7, height = 7)
ggsave(plot = stochiometric_plot, filename = 'updated_plots/xl_vs_nxl_with_mitocopies.pdf', width = 7, height = 7)
