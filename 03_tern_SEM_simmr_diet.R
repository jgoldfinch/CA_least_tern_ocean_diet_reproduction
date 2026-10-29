# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: run_simmr_diet_pipeline
# Description:  Runs SIMMR mixing models, performs a-posteriori source combining, 
#               generates bar plots, and outputs dietary proportions for SEM input.

################################## Function #####################################
run_simmr_diet_pipeline <- function(all_eggs_df, sources_df, tefs_df) {
  if (!requireNamespace("simmr", quietly = TRUE)) stop("Package 'simmr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")
  
  mcbc_eggs <- subset(all_eggs_df, site == "MCBCP")
  nbc_eggs  <- subset(all_eggs_df, site == "NBC")
  
  targets_all   <- all_eggs_df[c(1:2, 10, 9)]
  targets_MCBCP <- mcbc_eggs[c(1:2, 10, 9)]
  targets_NBC   <- nbc_eggs[c(1:2, 10, 9)]
  
  process_site_simmr <- function(targets, site_label) {
    simmr_in <- simmr::simmr_load(
      mixtures         = targets[, 3:4],
      source_names     = sources_df[, 1],
      source_means     = sources_df[, 2:3],
      source_sds       = sources_df[, 4:5],
      correction_means = tefs_df[, 2:3],
      correction_sds   = tefs_df[, 4:5],
      group            = as.factor(targets$year)
    )
    
    simmr_out <- simmr::simmr_mcmc(simmr_in)
    
    comb1 <- simmr::combine_sources(simmr_out, to_combine = c("anchovy", "topsmelt", "killifish", "perch"), new_source_name = "carbon-enriched")
    comb2 <- simmr::combine_sources(comb1, to_combine = c("juvenile rockfish", "saury", "krill"), new_source_name = "carbon-depleted")
    
    num_grps <- length(unique(targets$year))
    summ <- summary(comb2, type = "statistics", group = 1:num_grps)
    
    stats_df <- as.data.frame(summ$statistics) %>% tibble::rownames_to_column(var = "source")
    subset_df <- stats_df[2:3, ]
    
    mean_cols <- grep("group_\\d+\\.mean", names(subset_df), value = TRUE)
    cluster_df <- subset_df[, c("source", mean_cols)]
    
    long_df <- cluster_df %>% tidyr::pivot_longer(cols = -source, names_to = "group", values_to = "proportion")
    group_numbers <- gsub("group_|\\.mean", "", unique(long_df$group))
    years <- levels(as.factor(targets$year))
    
    map_df <- tibble::tibble(group = paste0("group_", group_numbers, ".mean"), year = years)
    plot_df <- dplyr::left_join(long_df, map_df, by = "group")
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = year, y = proportion, fill = source)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(x = "Year", y = "Proportion", fill = "Prey Cluster", title = paste("Source Proportions:", site_label)) +
      ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    ggplot2::ggsave(paste0("simmr_stacked_bar_", site_label, ".png"), plot = p, width = 8, height = 6, dpi = 300)
    return(plot_df)
  }
  
  df_mcbc  <- process_site_simmr(targets_MCBCP, "MCBCP")
  df_nbc   <- process_site_simmr(targets_NBC, "NBC")
  df_pooled<- process_site_simmr(targets_all, "Sites_Combined")
  
  format_site_props <- function(df, site_label) {
    df %>%
      dplyr::mutate(
        year = as.integer(as.character(year)),
        diet = dplyr::case_when(
          source == "carbon-enriched" ~ "C-enriched diet prop",
          source == "carbon-depleted" ~ "C-depleted diet prop",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(diet)) %>%
      dplyr::select(year, diet, proportion) %>%
      tidyr::pivot_wider(names_from = diet, values_from = proportion) %>%
      dplyr::rename_with(~ paste(.x, site_label), -year)
  }
  
  pooled   <- format_site_props(df_pooled, "Sites pooled")
  pendleton<- format_site_props(df_mcbc, "Pendleton")
  coronado <- format_site_props(df_nbc, "Coronado")
  
  carbon_proportions_table <- dplyr::full_join(pooled, pendleton, by = "year") %>%
    dplyr::full_join(coronado, by = "year") %>%
    dplyr::arrange(year)
  
  write.csv(carbon_proportions_table, "carbon_proportions_table.csv", row.names = FALSE)
  return(carbon_proportions_table)
}