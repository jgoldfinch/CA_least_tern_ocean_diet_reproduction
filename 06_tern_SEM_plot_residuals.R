# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: plot_sem_residuals
# Description:  Plots time-series and ACF residual diagnostic plots.

################################## Function #####################################
plot_sem_residuals <- function(resids, var_labels = NULL, max_lag = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  melt_to_long <- function(sub_obj, site_name) {
    mat <- as.matrix(sub_obj)
    v_names <- colnames(mat)
    do.call(rbind, lapply(v_names, function(v) {
      data.frame(obs_id = seq_len(nrow(mat)), site = site_name, variable = v, residual = as.numeric(mat[, v]), stringsAsFactors = FALSE)
    }))
  }
  
  df_long <- do.call(rbind, lapply(seq_along(resids), function(i) melt_to_long(resids[[i]], names(resids)[i])))
  df_long <- df_long %>% dplyr::filter(abs(residual) > 1e-5)
  
  p_raw <- ggplot2::ggplot(df_long, ggplot2::aes(x = obs_id, y = residual, color = site)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::facet_grid(variable ~ site, scales = "free_y") + ggplot2::theme_minimal()
  
  acf_df <- df_long %>%
    dplyr::group_by(site, variable) %>%
    dplyr::reframe(lag = 0:min(max_lag, dplyr::n() - 1), acf = as.numeric(stats::acf(residual, lag.max = max_lag, plot = FALSE)$acf), n = dplyr::n()) %>%
    dplyr::filter(lag > 0)
  
  ci_val <- 1.96 / sqrt(mean(acf_df$n))
  p_acf  <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf, fill = site)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = c(ci_val, -ci_val), linetype = "dashed", color = "blue") +
    ggplot2::geom_col(width = 0.3) + ggplot2::facet_grid(variable ~ site) + ggplot2::theme_minimal()
  
  return(list(raw_plot = p_raw, acf_plot = p_acf))
}