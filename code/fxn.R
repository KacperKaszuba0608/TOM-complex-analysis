# Function to perform missingness type assigning
assign_missing <- function(protein.ids, condition, lfq_intensity) {
  # NOTE: data frame must be in long format
  # protein.ids = column with protein IDs
  # condition = column conditions of samples (e.g. KO, WT)
  # lfq_intensity = column with LFQ intensity
  
  # verification of the fxn assumptions
  if (!is.factor(condition)) {rlang::abort("The `condition` columns is not a factor type!")}
  
  # occurences of protein IDs per condition
  occur <- dplyr::tibble(protein.ids) |>
    dplyr::group_by(protein.ids) |>
    dplyr::summarise(len = length(protein.ids)) |>
    dplyr::distinct(len) |>
    as.integer()
  
  # occurences of condition
  occur_per_cond <- dplyr::tibble(protein.ids, condition) |>
    dplyr::group_by(protein.ids, condition) |>
    dplyr::summarise(len = length(protein.ids))
  occur_per_cond <- unique(occur_per_cond$len)
  
  if (length(protein.ids)/occur != length(unique(protein.ids))) {rlang::abort("The protein IDs are not unique!")}
  if (!is.numeric(lfq_intensity)) {rlang::abort("The lfq intensities are not numeric!")}
  
  df <- data.frame("prot.IDs"=protein.ids, "condition"=condition, "lfq"=lfq_intensity)
  
  # Calculating number of missing values for each protein for each condition
  number_missing1 <- df |>
    dplyr::group_by(prot.IDs, condition) |>
    dplyr::summarise(no_NAs = sum(is.na(lfq))) |>
    dplyr::ungroup() 
  # |>
  #     dplyr::arrange(prot.IDs, desc(condition))
  
  # Assigning the missing type per condition
  df_miss <- merge(df, number_missing1, by=c('prot.IDs', 'condition'))
  
  df_miss <- df_miss |>
    mutate(missingness_per_cond = sapply(no_NAs, function(no) {
      dplyr::case_when(
        no == occur_per_cond ~ "all_NA",
        no == occur_per_cond-1 ~ "MNAR",
        no < occur_per_cond-1 & no != 0 ~ "MAR",
        no == 0 ~ "complete",
        TRUE ~ NA
      )
    })) |>
    select(-no_NAs)
  
  # Calculating number of missing values for each protein
  number_missing2 <- df |>
    dplyr::group_by(prot.IDs) |>
    dplyr::summarise(no_NAs = sum(is.na(lfq))) |>
    dplyr::ungroup()
  
  df_miss <- merge(df_miss, number_missing2, by='prot.IDs')
  
  # Assigining the missing type per protein
  df_miss <- df_miss |>
    mutate(missingness_per_prot = sapply(no_NAs, function(no) {
      dplyr::case_when(
        no == occur ~ "all_NA",
        no == occur-1 ~ "MNAR",
        no < occur-1 & no != 0 ~ "MAR",
        no == 0 ~ "complete",
        TRUE ~ NA
      )
    })) |>
    select(-no_NAs)
  
  # Final missing assigning
  suppressWarnings(
    missingness_final <- df_miss |>
      dplyr::select(prot.IDs, condition, missingness_per_cond) |>
      tidyr::pivot_wider(names_from = "condition", 
                         values_from="missingness_per_cond")
  )
  
  missingness_final[, 2] <- apply(missingness_final[,2], 1, function(type) unique(type |> unlist()))
  missingness_final[, 3] <- apply(missingness_final[,3], 1, function(type) unique(type |> unlist()))
  missingness_final$missingness <- apply(missingness_final, 1, function(row) {
    dplyr::case_when(
      (row[2] == "complete" & row[3] == "all_NA") | 
        (row[3] == "complete" & row[2] == "all_NA") | 
        (row[2] == "complete" & row[3] =="MNAR") | 
        (row[3] == "complete" & row[2] =="MNAR") ~ "MNAR",
      all(row[2:3] == "MAR") ~ "MAR",
      all(row[2:3] == "complete") ~ "complete",
      TRUE ~ NA
    )
  })
  
  df <- merge(df, missingness_final[,c('prot.IDs', "missingness")], by = 'prot.IDs')
  
  # return
  ret_list <- list(df_long = df, df_wide = missingness_final, 
                   missingness_long = df$missingness,
                   missingness_wide = missingness_final$missingness)
}

ttest <- function(df, grp1, grp2){ 
  x = df[grp1]
  y = df[grp2]
  x = as.numeric((x))
  y = as.numeric((y))
  results = t.test(x,y, 
                   alternative = "two.sided", #one-sided: "greater" is x > y
                   paired = F,
                   na.action=na.omit)
  results$p.value
}
