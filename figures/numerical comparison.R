# =============================================================================
# Approximation comparison: True CDF vs theoretical approximations
# Families: G (Gamma), NB (NegBin), LG (Log-Gaussian), IG (Inverse-Gaussian)
# =============================================================================

beta0      <- -2.0
alpha      <- 1.0
sigma2     <- 1.5
omega_vals <- c(0.05, 0.1, 1.0, 2.0)
n_sim      <- 20000
x_max      <- 30

library(ggplot2)
library(dplyr)
library(patchwork)
library(statmod)
library(grid)

set.seed(42)

# =============================================================================
# Core functions
# =============================================================================

prev_fun <- function(mu, omega, alpha)
  1 - (omega / (omega + mu * (1 - exp(-alpha))))^omega

marginal_moments <- function(beta0, sigma2, omega, alpha, n_mc = 50000) {
  S  <- rnorm(n_mc, 0, sqrt(sigma2))
  mu <- exp(beta0 + S)
  p  <- prev_fun(mu, omega, alpha)
  w  <- p / mean(p)
  
  mu_C     <- weighted.mean(alpha * mu / p, w)
  sigma2_C <- weighted.mean(
    alpha * mu * (1 + alpha) / p +
      alpha^2 * mu^2 / p * (1/omega + 1 - 1/p), w
  )
  
  list(mu_C = mu_C, sigma2_C = sigma2_C)
}

# Always return all families; invalid ones are NULL
build_cdfs <- function(mu_C, sigma2_C) {
  mu_Z  <- mu_C - 1
  var_Z <- sigma2_C
  
  out <- list(G = NULL, NB = NULL, LG = NULL, IG = NULL)
  if (!is.finite(mu_Z) || !is.finite(var_Z) || mu_Z <= 0 || var_Z <= 0) return(out)
  
  # G: Gamma (continuous)
  nu     <- var_Z / mu_Z
  lambda <- mu_Z^2 / var_Z
  if (is.finite(nu) && is.finite(lambda) && nu > 0 && lambda > 0) {
    out$G <- function(y) pgamma(pmax(y - 1.5, 0), shape = lambda, scale = nu)
  }
  
  # NB: NegBin (discrete)
  denom <- var_Z - mu_Z
  if (is.finite(denom) && denom > 0) {
    r <- mu_Z^2 / denom
    if (is.finite(r) && r > 0) {
      out$NB <- function(y) pnbinom(pmax(y - 1, 0), size = r, mu = mu_Z)
    }
  }
  
  # LG: Lognormal (continuous)
  cv2 <- var_Z / mu_Z^2
  if (is.finite(cv2) && cv2 > -1) {
    s2 <- log(1 + cv2)
    m  <- log(mu_Z) - s2/2
    if (is.finite(s2) && s2 >= 0 && is.finite(m)) {
      out$LG <- function(y) plnorm(pmax(y - 1.5, 1e-10), meanlog = m, sdlog = sqrt(s2))
    }
  }
  
  # IG: Inverse-Gaussian (continuous)
  lambda_ig <- mu_Z^3 / var_Z
  if (is.finite(lambda_ig) && lambda_ig > 0) {
    out$IG <- function(y) statmod::pinvgauss(pmax(y - 1.5, 1e-10),
                                             mean = mu_Z, shape = lambda_ig)
  }
  
  out
}

sim_marginal <- function(n, beta0, sigma2, omega, alpha) {
  out <- integer()
  while (length(out) < n) {
    S <- rnorm(n * 4, 0, sqrt(sigma2))
    W <- rnbinom(length(S), mu = exp(beta0 + S), size = omega)
    Y <- rpois(length(W), alpha * W)
    out <- c(out, Y[Y > 0])
  }
  out[1:n]
}

# =============================================================================
# Aesthetics
# =============================================================================

fam_levels <- c("True", "G", "NB", "LG", "IG")

fam_labels <- c(
  True = "True CDF",
  G    = "G",
  NB   = "NB",
  LG   = "LG",
  IG   = "IG"
)

fam_cols <- c(
  True = "grey20",
  G    = "#2166AC",
  NB   = "#E66101",
  LG   = "#1B7837",
  IG   = "#762A83"
)

fam_lty <- c(
  True = "solid",
  G    = "dashed",
  NB   = "solid",
  LG   = "dotdash",
  IG   = "dotted"
)

theme_paper <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.key.width = unit(1.6, "cm")
  )

# =============================================================================
# Main loop
# =============================================================================

plots <- list()
stat_rows <- list()
y_grid <- 1:x_max

for (omega in omega_vals) {
  
  y_sim <- sim_marginal(n_sim, beta0, sigma2, omega, alpha)
  Fn    <- ecdf(y_sim)
  
  mm   <- marginal_moments(beta0, sigma2, omega, alpha)
  cdfs <- build_cdfs(mm$mu_C, mm$sigma2_C)
  
  # True CDF from simulated data
  ecdf_df <- data.frame(
    y = y_grid,
    cdf = Fn(y_grid),
    family = factor("True", levels = fam_levels)
  )
  
  # Theoretical curves (only for non-NULL families)
  theo_df <- bind_rows(lapply(names(cdfs), function(fam) {
    if (is.null(cdfs[[fam]])) return(NULL)
    data.frame(
      y = y_grid,
      cdf = cdfs[[fam]](y_grid),   # <-- IMPORTANT: correct evaluation
      family = factor(fam, levels = fam_levels)
    )
  }))
  
  # Plot data combined (this also guarantees legend training if needed)
  plot_df <- bind_rows(ecdf_df, theo_df)
  
  p <- ggplot(plot_df, aes(x = y, y = cdf, colour = family, linetype = family)) +
    geom_step(data = subset(plot_df, family %in% c("True","NB")), linewidth = 0.8) +
    geom_line(data = subset(plot_df, family %in% c("G","LG","IG")), linewidth = 0.7) +
    scale_colour_manual(values = fam_cols, breaks = fam_levels, labels = fam_labels, drop = FALSE) +
    scale_linetype_manual(values = fam_lty, breaks = fam_levels, labels = fam_labels, drop = FALSE) +
    scale_x_continuous(limits = c(1, x_max), breaks = pretty(1:x_max, n = 5)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = bquote(omega == .(omega)), x = "y", y = "CDF") +
    theme_paper
  
  # We'll collect legend globally, so remove it per-panel
  plots[[as.character(omega)]] <- p + theme(legend.position = "none")
  
  # Stats: ALWAYS include G / NB / LG / IG rows
  for (fam in names(cdfs)) {
    if (is.null(cdfs[[fam]])) {
      stat_rows[[paste(omega, fam)]] <- data.frame(omega = omega, family = fam,
                                                   KS = NA_real_, MAD = NA_real_, CvM = NA_real_)
    } else {
      cdf_vals  <- cdfs[[fam]](y_grid)
      ecdf_vals <- Fn(y_grid)
      diff <- abs(ecdf_vals - cdf_vals)
      stat_rows[[paste(omega, fam)]] <- data.frame(omega = omega, family = fam,
                                                   KS = max(diff), MAD = median(diff), CvM = mean(diff^2))
    }
  }
}

# =============================================================================
# Combine with a guaranteed legend (patchwork guide collection)
# =============================================================================

# Create ONE plot that still has the legend, to be collected:
p_with_legend <- plots[[1]] +
  theme(legend.position = "bottom")

# keep legend only on first plot
plots2 <- plots
plots2[[1]] <- plots2[[1]] + theme(legend.position = "bottom")

final <- wrap_plots(plots2, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))
ggsave("gamma_approximation_check.pdf", final, width = 10, height = 9, device = cairo_pdf)
ggsave("gamma_approximation_check.png", final, width = 10, height = 9, dpi = 300)

cat("Saved gamma_approximation_check.pdf / .png\n")

# =============================================================================
# Console table (Gamma INCLUDED)
# =============================================================================

approx_levels <- c("G","NB","LG","IG")

stat_df <- bind_rows(stat_rows) |>
  mutate(family = factor(family, levels = approx_levels)) |>
  arrange(omega, family)

fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))

cat("\n=== Goodness of approximation (vs True CDF from simulated data) ===\n")
cat(sprintf("%-6s  %-6s  %7s  %7s  %7s\n", "omega", "Family", "KS", "MAD", "CvM"))
cat(strrep("-", 46), "\n")
for (i in seq_len(nrow(stat_df))) {
  r <- stat_df[i, ]
  cat(sprintf("%-6.2f  %-6s  %7s  %7s  %7s\n",
              r$omega, as.character(r$family), fmt(r$KS), fmt(r$MAD), fmt(r$CvM)))
}
cat(strrep("-", 46), "\n")
cat("KS: max|True-F|  MAD: mean|True-F|  CvM: mean(True-F)^2\n")