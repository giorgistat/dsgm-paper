rm(list = ls())
library(RiskMap)
anpit <- readRDS("anpit_hk_mult.rds")

p_anpit <- plot_AnPIT(anpit, mode = "average", which = "conditional")
p_anpit <- p_anpit + ggplot2::ggtitle("") +
                     ggplot2::ylab("AnPIT(u)") +
                     ggplot2::xlab("u")

pdf("hk_anpit.pdf")
p_anpit
dev.off()

plot_zero_calibration <- function(object,
                                  model_name = NULL,
                                  n_bins = 10,
                                  combine_panels = FALSE,
                                  by_location = FALSE,
                                  title1 = "Positivity-gate reliability: predicted vs observed P(Y > 0)",
                                  xlab1 = "Mean predicted P(Y > 0) in bin",
                                  ylab1 = "Observed fraction Y > 0 in bin",
                                  ylim1 = c(0, 1),
                                  title2 = "Fraction of zeros: observed vs model-implied, by test fold",
                                  xlab2 = "Test fold",
                                  ylab2 = "Fraction Y = 0",
                                  ylim2 = NULL,
                                  title3 = "Observed vs predicted prevalence by location",
                                  xlab3 = "Predicted prevalence (location mean)",
                                  ylab3 = "Observed prevalence (location mean)",
                                  xlim3 = c(0, 1),
                                  ylim3 = c(0, 1)) {
  missing_title1 <- missing(title1)
  missing_title2 <- missing(title2)

  if (!inherits(object, "RiskMap.spatial.cv"))
    stop("`object` must be a 'RiskMap.spatial.cv' produced by assess_pp().")

  all_models <- names(object$model)
  if (!is.null(model_name)) {
    if (!model_name %in% all_models)
      stop("Model name '", model_name, "' not found in `object$model`.")
    all_models <- model_name
  }

  ## pool obs_ind / pred_prob / loc_id across test-set folds, per model
  make_pool <- function(mname) {
    m <- object$model[[mname]]
    if (is.null(m$pos_cal))
      stop("Model '", mname, "' has no `pos_cal` — this diagnostic is only ",
           "produced for DSGM/intprev models by the updated assess_pp().")

    obs_ind   <- unlist(lapply(m$pos_cal, `[[`, "obs_ind"))
    pred_prob <- unlist(lapply(m$pos_cal, `[[`, "pred_prob"))
    fold_id   <- rep(seq_along(m$pos_cal),
                     vapply(m$pos_cal, function(x) length(x$obs_ind), integer(1)))

    if (by_location) {
      if (is.null(m$pos_cal[[1]]$loc_id))
        stop("Model '", mname, "' has no `loc_id` in `pos_cal` — re-run ",
             "assess_pp() with the location-ID edit to use by_location = TRUE.")
      loc_id <- unlist(lapply(m$pos_cal, `[[`, "loc_id"))
      ## loc_id is only unique *within* a fold, so combine with fold_id
      ## to get a globally unique location key across folds.
      loc_key <- paste(fold_id, loc_id, sep = "_")
    } else {
      loc_key <- NA_character_
    }

    data.frame(model = mname, obs = obs_ind, pred = pred_prob,
               fold = fold_id, loc_key = loc_key)
  }

  pooled <- do.call(rbind, lapply(all_models, make_pool))

  ## If by_location: collapse to one row per (model, fold, location) by
  ## averaging obs/pred within location first, so every location counts
  ## once regardless of how many individuals were sampled there.
  if (by_location) {
    pooled <- pooled %>%
      dplyr::group_by(model, fold, loc_key) %>%
      dplyr::summarize(
        obs  = mean(obs,  na.rm = TRUE),
        pred = mean(pred, na.rm = TRUE),
        .groups = "drop"
      )
  }

  ## per-fold summary (headline numbers: observed vs model-implied zero fraction).
  ## Computed from `pooled` so it respects by_location -- the obs_frac/pred_frac
  ## stored directly in pos_cal are always individual-weighted, so we recompute
  ## here rather than reuse them when by_location = TRUE.
  fold_summary <- pooled %>%
    dplyr::group_by(model, fold) %>%
    dplyr::summarize(
      obs_frac  = mean(obs,  na.rm = TRUE),
      pred_frac = mean(pred, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    as.data.frame()
  fold_summary$obs_zero_frac  <- 1 - fold_summary$obs_frac
  fold_summary$pred_zero_frac <- 1 - fold_summary$pred_frac

  ## reliability bins: bin by predicted P(Y>0), compare to observed frequency
  pooled$bin <- cut(pooled$pred, breaks = seq(0, 1, length.out = n_bins + 1),
                    include.lowest = TRUE)

  reliability <- pooled %>%
    dplyr::group_by(model, bin) %>%
    dplyr::summarize(
      pred_mean = mean(pred, na.rm = TRUE),
      obs_mean  = mean(obs,  na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    )

  if (by_location) {
    if (missing_title1) title1 <- paste(title1, "(location-level)")
    if (missing_title2) title2 <- paste(title2, "(location-level)")
  }

  id_line <- ggplot2::geom_abline(intercept = 0, slope = 1,
                                  linetype = "dashed", colour = "red")

  ## Location-level scatter: observed vs predicted prevalence per location,
  ## coloured by fold/test set. Only meaningful when by_location = TRUE,
  ## since `pooled` then has exactly one row per (model, fold, location).
  p_location <- NULL
  if (by_location) {
    p_location <- ggplot2::ggplot(pooled,
                                  ggplot2::aes(pred, obs, colour = factor(fold))) +
      ggplot2::geom_point(alpha = 0.8) +
      id_line +
      ggplot2::coord_cartesian(xlim = xlim3, ylim = ylim3) +
      ggplot2::labs(title = title3,
                    x = xlab3,
                    y = ylab3,
                    colour = "Test set") +
      ggplot2::theme_minimal() +
      { if (length(all_models) > 1) ggplot2::facet_wrap(~model) else NULL }
  }

  ## Only connect points with a line for models that have >=2 populated
  ## bins -- geom_line() warns (and draws nothing useful) for a group with
  ## a single observation, which can happen with small held-out samples.
  reliability_line <- reliability %>%
    dplyr::group_by(model) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::ungroup()

  p_reliability <- ggplot2::ggplot(reliability,
                                   ggplot2::aes(pred_mean, obs_mean,
                                                size = n,
                                                colour = if (length(all_models) > 1) model else NULL)) +
    ggplot2::geom_point(alpha = 0.8) +
    { if (nrow(reliability_line) > 0)
      ggplot2::geom_line(data = reliability_line,
                         ggplot2::aes(group = model), alpha = 0.5, linewidth = 0.6)
      else NULL } +
    id_line +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim1) +
    ggplot2::labs(title = title1,
                  x = xlab1,
                  y = ylab1,
                  size = "n points") +
    ggplot2::theme_minimal() +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Model"))

  ## per-fold check on the zero fraction specifically
  p_fold <- ggplot2::ggplot(fold_summary,
                            ggplot2::aes(x = factor(fold), colour = model)) +
    ggplot2::geom_point(ggplot2::aes(y = obs_zero_frac, shape = "Observed"), size = 3) +
    ggplot2::geom_point(ggplot2::aes(y = pred_zero_frac, shape = "Predicted (model)"), size = 3) +
    ggplot2::geom_segment(ggplot2::aes(xend = factor(fold), y = obs_zero_frac, yend = pred_zero_frac),
                          alpha = 0.4) +
    ggplot2::labs(title = title2,
                  x = xlab2, y = ylab2, shape = "") +
    ggplot2::theme_minimal() +
    { if (!is.null(ylim2)) ggplot2::coord_cartesian(ylim = ylim2) else NULL }

  if (combine_panels && requireNamespace("gridExtra", quietly = TRUE)) {
    if (!is.null(p_location)) {
      gridExtra::grid.arrange(p_reliability, p_fold, p_location, ncol = 2)
    } else {
      gridExtra::grid.arrange(p_reliability, p_fold, ncol = 2)
    }
  } else {
    print(p_reliability)
    print(p_fold)
    if (!is.null(p_location)) print(p_location)
  }

  invisible(list(reliability = reliability, fold_summary = fold_summary,
                 pooled = pooled, by_location = by_location,
                 p_reliability = p_reliability, p_fold = p_fold,
                 p_location = p_location))
}
library(dplyr)
res <- plot_zero_calibration(anpit, by_location = TRUE,
                             xlim3 = c(0, 0.25), ylim3 = c(0, 0.3),
                             combine_panels = FALSE,
                             xlab3 = "Predicted hookworm prevalence",
                             ylab3 = "Empirical hookworm prevalence",
                             title3 = "")

p_zero <- res$p_location   # the location-scatter plot object

pdf("hk_zero_plot.pdf")
print(p_zero)
dev.off()
pdf("hk_zero_plot.pdf")
p_zero
dev.off()

library(sf)
library(ggplot2)

test_sf <- do.call(rbind, lapply(seq_along(anpit$test_set), function(i) {
  d <- anpit$test_set[[i]]
  d$test_set <- factor(i)
  d
}))

## one row per unique (location, test_set) -- deliberately keeps a location
## in every fold it belongs to, since we're faceting rather than colouring
test_sf_unique <- test_sf[
  !duplicated(paste(sf::st_as_text(sf::st_geometry(test_sf)), test_sf$test_set)),
]

ggplot() +
  geom_sf(data = test_sf_unique, colour = "steelblue", size = 2, alpha = 0.85) +
  facet_wrap(~test_set) +
  theme_minimal() +
  labs(title = "Held-out locations by test set")
