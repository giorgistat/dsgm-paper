library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(metR)
library(RiskMap)

# =============================================================================
# 1. LOAD SAVED FITS / PREDICTIONS
# =============================================================================
fit_dsgm <- readRDS("fit_sth_dsgm.rds")

pred_prev_dsgm      <- readRDS("pred_prev_dsgm.rds")       # DSGM: any (target$prevalence)
pred_intensity_dsgm <- readRDS("pred_intensity_dsgm.rds")  # DSGM: light/moderate/high
pred_prev_dast      <- readRDS("pred_prev_dast.rds")       # DAST: list(any, light, moderate, high)

grid_pred <- readRDS("ken_grid_pred.rds")

# =============================================================================
# 2. GRID COORDINATES AND HELPER TO BUILD A PLOTTING DATA FRAME
# =============================================================================
if (inherits(grid_pred, "sf")) {
  grid_xy <- st_coordinates(grid_pred)
} else {
  grid_xy <- as.matrix(grid_pred[, c("x", "y")])
}

make_df <- function(mean_vals) {
  r  <- rast(cbind(grid_xy, mean_vals), type = "xyz", crs = "EPSG:32736")
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "prevalence"
  df$prevalence <- df$prevalence * 100
  df
}

# DSGM: any / light / moderate / high
df_dsgm_any      <- make_df(pred_prev_dsgm$target$prevalence$mean)
df_dsgm_light    <- make_df(pred_intensity_dsgm$target$light_prev$mean)
df_dsgm_moderate <- make_df(pred_intensity_dsgm$target$moderate_prev$mean)
df_dsgm_high     <- make_df(pred_intensity_dsgm$target$high_prev$mean)

# DAST: any / light / moderate / high
df_dast_any      <- make_df(pred_prev_dast$any$target$prevalence$mean)
df_dast_light    <- make_df(pred_prev_dast$light$target$prevalence$mean)
df_dast_moderate <- make_df(pred_prev_dast$moderate$target$prevalence$mean)
df_dast_high     <- make_df(pred_prev_dast$high$target$prevalence$mean)

# =============================================================================
# 3. PER-BAND COLOUR SCALE (DSGM + DAST share a scale within each band,
#    but the scale range adapts separately per band)
# =============================================================================
zlim_any      <- range(c(df_dsgm_any$prevalence,      df_dast_any$prevalence),      na.rm = TRUE)
zlim_light    <- range(c(df_dsgm_light$prevalence,    df_dast_light$prevalence),    na.rm = TRUE)
zlim_moderate <- range(c(df_dsgm_moderate$prevalence, df_dast_moderate$prevalence), na.rm = TRUE)
zlim_high     <- range(c(df_dsgm_high$prevalence,     df_dast_high$prevalence),     na.rm = TRUE)

# =============================================================================
# 4. SHARED THEME
# =============================================================================
theme_map <- theme_bw(base_size = 11) +
  theme(panel.grid   = element_blank(),
        legend.position = "right",
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        plot.title   = element_text(size = 13, hjust = 0.5, margin = margin(b = 3)),
        plot.margin  = margin(2, 1, 2, 1),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11))

# =============================================================================
# 5. PANEL BUILDER (takes the band-specific colour limits)
# =============================================================================
make_panel <- function(df, title, zlim) {
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = prevalence)) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "YlOrRd"),
      limits  = zlim,
      name    = "Prev. (%)"
    ) +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_map
}

# =============================================================================
# 6. BUILD ALL 8 PANELS (rows: any / light / moderate / high;
#    left column: DSGM, right column: DAST)
# =============================================================================
p_dsgm_any      <- make_panel(df_dsgm_any,      "DSGM \u2014 any intensity",      zlim_any)
p_dast_any      <- make_panel(df_dast_any,      "DAST \u2014 any intensity",      zlim_any)
p_dsgm_light    <- make_panel(df_dsgm_light,    "DSGM \u2014 light intensity",    zlim_light)
p_dast_light    <- make_panel(df_dast_light,    "DAST \u2014 light intensity",    zlim_light)
p_dsgm_moderate <- make_panel(df_dsgm_moderate, "DSGM \u2014 moderate intensity", zlim_moderate)
p_dast_moderate <- make_panel(df_dast_moderate, "DAST \u2014 moderate intensity", zlim_moderate)
p_dsgm_high     <- make_panel(df_dsgm_high,     "DSGM \u2014 high intensity",     zlim_high)
p_dast_high     <- make_panel(df_dast_high,     "DAST \u2014 high intensity",     zlim_high)

# =============================================================================
# 7. BUILD EACH DSGM/DAST PAIR AS ITS OWN UNIT.
#    Each pair shares a legend (band-specific colour scale).
# =============================================================================
make_pair <- function(dsgm_panel, dast_panel) {
  pair <- (dsgm_panel | dast_panel) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  wrap_elements(full = pair)
}

pair_any      <- make_pair(p_dsgm_any,      p_dast_any)
pair_light    <- make_pair(p_dsgm_light,    p_dast_light)
pair_moderate <- make_pair(p_dsgm_moderate, p_dast_moderate)
pair_high     <- make_pair(p_dsgm_high,     p_dast_high)

# =============================================================================
# 8. COMBINE — 2 rows x 2 pair-columns (any & light on row 1;
#    moderate & high on row 2).
# =============================================================================
fig <- (pair_any | pair_light) / (pair_moderate | pair_high)

fig

# =============================================================================
# 9. SIZE OUTPUT TO MATCH THE DATA'S ASPECT RATIO (avoids blank margins left
#    by coord_equal() when the canvas doesn't match the true extent shape)
# =============================================================================
data_asp   <- diff(range(grid_xy[, 2])) / diff(range(grid_xy[, 1]))  # height / width
panel_w    <- 4    # inches per map panel (excluding legend)
legend_w   <- 1    # inches reserved for each pair's legend
fig_width  <- 2 * (2 * panel_w + legend_w)
fig_height <- 2 * panel_w * data_asp

ggsave("dsgm_dast_comparison_bands.png", fig, width = fig_width, height = fig_height, dpi = 900)
ggsave("dsgm_dast_comparison_bands.pdf", fig, width = fig_width, height = fig_height, device = cairo_pdf)
