# =============================================================================
# MDA temporal decay: worm burden vs prevalence
#
# Core message:
#   MDA reduces worm burden by the SAME relative amount regardless of baseline.
#   But the same relative worm burden reduction produces DIFFERENT relative
#   prevalence reductions depending on baseline endemicity.
#
# Two panels, time on x-axis:
#   Panel A: Relative worm burden  mu(t)/mu*  — one curve for all baselines
#   Panel B: Relative prevalence   p(t)/p*    — different curve per baseline
# =============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)

# ------ Model parameters ------------------------------------------------------
omega   <- 0.5     # NB aggregation
alpha   <- 1.0     # per-worm output rate
alpha_W <- 0.6     # immediate worm burden reduction from MDA
gamma_W <- 2.0     # recovery timescale (years)

# ------ Functions -------------------------------------------------------------
# MDA decay: proportion of baseline worm burden remaining at time t post-MDA
phi <- function(t) 1 - alpha_W * exp(-t / gamma_W)

# Prevalence from mean worm burden
prev <- function(mu) 1 - (omega / (omega + mu * (1 - exp(-alpha))))^omega

# Invert: mu from prevalence
mu_from_prev <- function(p0) {
  omega * ((1 - p0)^(-1 / omega) - 1) / (1 - exp(-alpha))
}

# ------ Three baseline levels -------------------------------------------------
baselines <- data.frame(
  label = c("Low (p* = 10%)", "Medium (p* = 40%)", "High (p* = 75%)"),
  p_star = c(0.10, 0.40, 0.75)
) |>
  mutate(mu_star = mu_from_prev(p_star))

cols <- c("Low (p* = 10%)"    = "#2166AC",
          "Medium (p* = 40%)" = "#F4A582",
          "High (p* = 75%)"   = "#B2182B")

# ------ Time grid (years since MDA) ------------------------------------------
t_seq <- seq(0, 8, length.out = 400)

# ------ Build data frame ------------------------------------------------------
df <- lapply(seq_len(nrow(baselines)), function(k) {
  mu_star <- baselines$mu_star[k]
  p_star  <- baselines$p_star[k]
  data.frame(
    t            = t_seq,
    label        = baselines$label[k],
    rel_mu       = phi(t_seq),                               # same for all k
    rel_prev     = prev(mu_star * phi(t_seq)) / p_star       # differs by k
  )
}) |> bind_rows() |>
  mutate(label = factor(label, levels = baselines$label))

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank())

# ------ Panel A: relative worm burden ----------------------------------------
# All three curves are identical — plot just one line with a note
df_a <- df |> filter(label == "Low (p* = 10%)")

p_a <- ggplot(df_a, aes(x = t, y = rel_mu * 100)) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 1.1, colour = "grey20") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(limits = c(38, 105), breaks = seq(0, 100, 20)) +
  labs(x     = "Years since MDA",
       y     = expression("Worm burden relative to baseline (%)"),
       title = "(a) Relative change in worm burden") +
  annotate("text", x = 5, y = 15,
           label = "Identical decay regardless\nof baseline endemicity",
           size = 3.2, colour = "grey30") +
  theme_paper

# ------ Panel B: relative prevalence -----------------------------------------
p_b <- ggplot(df, aes(x = t, y = rel_prev * 100, colour = label)) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(limits = c(38, 105), breaks = seq(0, 100, 20)) +
  labs(x     = "Years since MDA",
       y     = expression("Prevalence relative to baseline (%)"),
       title = "(b) Relative change in prevalence") +
  theme_paper

# ------ Combine ---------------------------------------------------------------
combined <- p_a | p_b

ggsave("mda_prevalence_illustration.pdf",
       plot = combined, width = 11, height = 5, device = cairo_pdf)
ggsave("mda_prevalence_illustration.png",
       plot = combined, width = 11, height = 5, dpi = 300)

# ------ Print reduction at peak effect (t -> 0+) and at t = 2 years ----------
cat("\n=== Worm burden and prevalence reductions at t = 0+ (immediate) ===\n")
cat(sprintf("Worm burden: %.1f%% of baseline (same for all)\n",
            phi(0.001) * 100))
cat(sprintf("\n%-22s  Prevalence at t=0+  Prevalence at t=2y\n", "Baseline"))
for (k in seq_len(nrow(baselines))) {
  mu_star <- baselines$mu_star[k]
  p_star  <- baselines$p_star[k]
  p0 <- prev(mu_star * phi(0.001)) / p_star * 100
  p2 <- prev(mu_star * phi(2))     / p_star * 100
  cat(sprintf("%-22s  %4.1f%%               %4.1f%%\n",
              baselines$label[k], p0, p2))
}
cat(sprintf("\nalpha_W=%.1f  gamma_W=%.1f  omega=%.1f  alpha=%.1f\n",
            alpha_W, gamma_W, omega, alpha))