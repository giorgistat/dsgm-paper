# =============================================================================
# Functional relationship between MF and antigen (ICT) prevalence
#
# p_mf  = 1 - [omega / (omega + mu*(1-exp(-alpha)))]^omega
# p_ict = gamma * {1 - [omega / (omega + mu)]^omega}
#
# The explicit curve p_ict = g(p_mf) is obtained by eliminating mu:
#   mu = omega * [(1 - p_mf)^(-1/omega) - 1] / (1 - exp(-alpha))
#
# Parameters varied:
#   alpha  — per-worm MF detection rate (controls MF sensitivity vs antigen)
#   omega  — NB aggregation (controls shape of the curve)
#   gamma  — ICT test sensitivity (fixed at 0.95)
# =============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)

gamma_sens <- 0.95   # fixed ICT sensitivity

# Vectorised curve function: given p_mf, omega, alpha, gamma -> p_ict
g_curve <- function(p_mf, omega, alpha, gamma = gamma_sens) {
  # invert eq (prob_gr_zero) to get mu
  mu <- omega * ((1 - p_mf)^(-1 / omega) - 1) / (1 - exp(-alpha))
  # antigen prevalence from eq (pict)
  gamma * (1 - (omega / (omega + mu))^omega)
}

# Grid of p_mf values
p_mf_seq <- seq(0.001, 0.999, length.out = 500)

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.size  = unit(0.45, "cm"))

# =============================================================================
# PANEL A — vary alpha, fix omega
# alpha controls how readily MF are detected per worm.
# Large alpha: MF test almost as sensitive as antigen -> curve close to diagonal
# Small alpha: MF test misses many infected hosts -> curve bows above diagonal
# =============================================================================
omega_fixed <- 0.5
alphas      <- c(0.1, 0.25, 0.5, 1.0, 2.0)
alpha_labels <- paste0("alpha == ", alphas)

df_a <- expand.grid(p_mf = p_mf_seq, alpha = alphas) |>
  mutate(p_ict  = g_curve(p_mf, omega = omega_fixed, alpha = alpha),
         alpha_f = factor(alpha, levels = alphas,
                          labels = paste0("\u03b1 = ", alphas)))

p_a <- ggplot(df_a, aes(x = p_mf * 100, y = p_ict * 100,
                        colour = alpha_f, linetype = alpha_f)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey60", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  scale_colour_brewer(palette = "Dark2", name = expression(alpha)) +
  scale_linetype_manual(values = c("solid","longdash","dashed","dotdash","dotted"),
                        name = expression(alpha)) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x     = "MF prevalence (%)",
       y     = "ICT prevalence (%)",
       title = bquote("(a) Varying " * alpha * "  (  " * omega == .(omega_fixed) * "  )")) +
  theme_paper +
  theme(legend.position = "right")

# =============================================================================
# PANEL B — vary omega, fix alpha
# omega controls worm burden aggregation.
# Small omega: highly aggregated -> a few heavily infected hosts drive both
#              diagnostics, strong nonlinearity in the curve
# Large omega: less aggregated (near Poisson) -> more linear relationship
# =============================================================================
alpha_fixed <- 0.1
omegas      <- c(0.1, 0.3, 0.5, 1.0, 5.0)

df_b <- expand.grid(p_mf = p_mf_seq, omega = omegas) |>
  mutate(p_ict  = g_curve(p_mf, omega = omega, alpha = alpha_fixed),
         omega_f = factor(omega, levels = omegas,
                          labels = paste0("\u03c9 = ", omegas)))

p_b <- ggplot(df_b, aes(x = p_mf * 100, y = p_ict * 100,
                        colour = omega_f, linetype = omega_f)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey60", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  scale_colour_brewer(palette = "Dark2", name = expression(omega)) +
  scale_linetype_manual(values = c("solid","longdash","dashed","dotdash","dotted"),
                        name = expression(omega)) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x     = "MF prevalence (%)",
       y     = "ICT prevalence (%)",
       title = bquote("(b) Varying " * omega * "  (  " * alpha == .(alpha_fixed) * "  )")) +
  theme_paper +
  theme(legend.position = "right")

# =============================================================================
# Combine and save
# =============================================================================
combined <- p_a | p_b

ggsave("diagnostic_curves.pdf",
       plot   = combined,
       width  = 10, height = 4.5,
       device = cairo_pdf)

ggsave("diagnostic_curves.png",
       plot   = combined,
       width  = 10, height = 4.5, dpi = 300)

cat("Saved diagnostic_curves.pdf / .png\n")
cat(sprintf("gamma (ICT sensitivity) fixed at %.2f\n", gamma_sens))