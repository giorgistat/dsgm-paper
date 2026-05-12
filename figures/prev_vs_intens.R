library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---------------------------------------------------------------------------
# Prevalence-intensity curve under the DSGM
#
#   p(lambda) = 1 - [ omega / (omega + lambda*(1-exp(-alpha))/alpha) ]^omega
#
# lambda = mean intensity (EPG), alpha = egg detection rate per worm,
# omega  = aggregation parameter
#
# Note: for fixed lambda, smaller alpha implies larger implied worm burden
# (mu = lambda/alpha), hence higher prevalence. Curves are ordered so that
# smaller alpha sits above larger alpha.
# ---------------------------------------------------------------------------

prev_intensity <- function(lambda, alpha, omega) {
  1 - (omega / (omega + (lambda / alpha) * (1 - exp(-alpha))))^omega
}

lambda_seq  <- seq(0, 500, length.out = 500)
omega_vals  <- c(0.01, 0.05, 0.1, 0.5, 1.0)
alpha_fixed <- c(10, 50, 100)
panel_labs  <- c("(a)", "(b)", "(c)")

colours     <- c("#1b4f72", "#2980b9", "#27ae60", "#e67e22", "#c0392b")
linetypes   <- c("solid", "dashed", "dotdash", "longdash", "twodash")

# ---------------------------------------------------------------------------
# Build one panel per fixed alpha value
# ---------------------------------------------------------------------------
make_panel <- function(alpha_val, panel_label) {

  df <- expand.grid(lambda = lambda_seq, omega = omega_vals) |>
    mutate(
      prevalence = prev_intensity(lambda, alpha = alpha_val, omega = omega),
      omega_lab  = factor(paste0("\u03c9 = ", omega),
                          levels = paste0("\u03c9 = ", omega_vals))
    )

  ggplot(df, aes(x = lambda, y = prevalence,
                 colour = omega_lab, linetype = omega_lab)) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_linetype_manual(values = linetypes, name = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(
      x     = "Mean intensity \u03bb (EPG)",
      y     = "Prevalence",
      title = bquote(.(panel_label)~"  "~alpha == .(alpha_val))
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key.width  = unit(1.2, "cm"),
      legend.text       = element_text(size = 10),
      plot.title        = element_text(size = 11, hjust = 0)
    )
}

panels <- mapply(make_panel, alpha_fixed, panel_labs, SIMPLIFY = FALSE)

# ---------------------------------------------------------------------------
# Combine: shared legend from first panel
# ---------------------------------------------------------------------------
fig <- wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
ggsave(
  filename = "figures/prev_intensity_curve.pdf",
  plot     = fig,
  width    = 13,
  height   = 5.2,
  device   = cairo_pdf
)

ggsave(
  filename = "prev_intensity_curve.png",
  plot     = fig,
  width    = 13,
  height   = 5.2,
  dpi      = 300
)

message("Saved: prev_intensity_curve.pdf and prev_intensity_curve.png")
