library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(metR)

# =============================================================================
# 1. BUILD RASTERS FROM GRID PREDICTIONS
# =============================================================================
if (inherits(grid_pred, "sf")) {
  grid_xy <- st_coordinates(grid_pred)
} else {
  grid_xy <- as.matrix(grid_pred[, c("x", "y")])
}

r_dsgm <- rast(cbind(grid_xy, pred_prev_dsgm$target$prevalence$mean),
               type = "xyz", crs = "EPSG:32736")
r_dast <- rast(cbind(grid_xy, pred_prev_dast$target$prevalence$mean),
               type = "xyz", crs = "EPSG:32736")

df_dsgm <- as.data.frame(r_dsgm, xy = TRUE); names(df_dsgm)[3] <- "prevalence"
df_dast <- as.data.frame(r_dast, xy = TRUE); names(df_dast)[3] <- "prevalence"

df_dsgm$prevalence <- df_dsgm$prevalence * 100
df_dast$prevalence <- df_dast$prevalence * 100

# =============================================================================
# 2. COMMON COLOUR SCALE AND CONTOUR BREAKS
# =============================================================================
zlim <- range(c(df_dsgm$prevalence, df_dast$prevalence), na.rm = TRUE)
contour_breaks <- c(0, 2.5, 5, 10, 20, 30)

# =============================================================================
# 3. SURVEY LOCATIONS
# =============================================================================
survey_xy <- as.data.frame(st_coordinates(st_geometry(fit_dsgm$data_sf)))
names(survey_xy) <- c("x", "y")
survey_xy <- unique(survey_xy)

# =============================================================================
# 4. SHARED THEME / SCALE HELPERS
# =============================================================================
theme_map <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.position = "right")

fill_scale <- scale_fill_gradientn(
  colours = brewer.pal(9, "YlOrRd"),
  limits  = zlim,
  name    = "Prevalence (%)"
)

# =============================================================================
# 5. PANEL 1 — DSGM map
# =============================================================================
p1 <- ggplot() +
  geom_raster(data = df_dsgm, aes(x, y, fill = prevalence)) +
  geom_point(data = survey_xy, aes(x, y),
             shape = 21, fill = "black", colour = "black",
             size = 1.3, stroke = 0.5) +
  fill_scale +
  geom_contour(data = df_dsgm, aes(x, y, z = prevalence),
               breaks = contour_breaks, colour = "grey20", linewidth = 0.3) +
  geom_text_contour(data = df_dsgm,
                    aes(x, y, z = prevalence, label = after_stat(paste0(level, "%"))),
                    breaks = contour_breaks, colour = "grey20", size = 2.5,
                    stroke = 0.15, stroke.colour = "white",
                    skip = 1, min.size = 20) +
  coord_equal() +
  labs(title = "DSGM prevalence", x = NULL, y = NULL) +
  theme_map +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# =============================================================================
# 6. PANEL 2 — DAST map
# =============================================================================
p2 <- ggplot() +
  geom_raster(data = df_dast, aes(x, y, fill = prevalence)) +
  geom_point(data = survey_xy, aes(x, y),
             shape = 21, fill = "black", colour = "black",
             size = 1.3, stroke = 0.5) +
  fill_scale +
  geom_contour(data = df_dast, aes(x, y, z = prevalence),
               breaks = contour_breaks, colour = "grey20", linewidth = 0.3) +
  geom_text_contour(data = df_dast,
                    aes(x, y, z = prevalence, label = after_stat(paste0(level, "%"))),
                    breaks = contour_breaks, colour = "grey20", size = 2.5,
                    stroke = 0.15, stroke.colour = "white",
                    skip = 1, min.size = 20) +
  coord_equal() +
  labs(title = "DAST prevalence", x = NULL, y = NULL) +
  theme_map +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# =============================================================================
# 7. PANEL 3 — DSGM vs DAST scatter with identity line
# =============================================================================
scatter_df <- data.frame(
  dsgm = pred_prev_dsgm$target$prevalence$mean * 100,
  dast = pred_prev_dast$target$prevalence$mean * 100
)

p3 <- ggplot(scatter_df, aes(dsgm, dast)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.6) +
  geom_point(colour = "steelblue", size = 0.8, alpha = 0.6) +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 30)) +
  coord_equal() +
  labs(title = "", x = "DSGM prevalence (%)", y = "DAST prevalence (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

# =============================================================================
# 8. COMBINE — maps on top row, scatter below spanning full width
# =============================================================================
fig <- (p1 | p2) + plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")

fig

ggsave("dsgm_dast_comparison.png", fig, width = 9, height = 9, dpi = 300)
ggsave("dsgm_dast_comparison.pdf", fig, width = 9, height = 9, device = cairo_pdf)
