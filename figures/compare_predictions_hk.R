library(terra)
library(sf)

# =============================================================================
# 1. BUILD RASTERS FROM GRID PREDICTIONS
# =============================================================================
# Assumption: grid_pred is a regular grid (sf object or data.frame with x, y
# columns in the same UTM CRS as fit_dsgm$data_sf). If pred_prev_dsgm/pred_prev_dast
# are outputs of pred_target_grid(), their $target$prevalence$mean should align
# row-for-row with grid_pred.

if (inherits(grid_pred, "sf")) {
  grid_xy <- st_coordinates(grid_pred)
} else {
  grid_xy <- as.matrix(grid_pred[, c("x", "y")])
}

make_raster <- function(xy, z) {
  df <- data.frame(x = xy[, 1], y = xy[, 2], z = z)
  terra::rast(df, type = "xyz", crs = "EPSG:32736")  # UTM 36S, matches fit_dsgm$data_sf
}

r_dsgm <- rast(cbind(st_coordinates(grid_xy), pred_prev_dsgm$target$prevalence$mean),
               type = "xyz")
r_dast <- rast(cbind(st_coordinates(grid_xy), pred_prev_dast$target$prevalence$mean),
               type = "xyz")

# =============================================================================
# 2. COMMON COLOUR SCALE
# =============================================================================
zlim <- range(c(values(r_dsgm), values(r_dast)), na.rm = TRUE)
breaks <- seq(zlim[1], zlim[2], length.out = 100)
col_pal <- terrain.colors(99)  # swap for viridis::viridis(99) if preferred

# =============================================================================
# 3. SURVEY LOCATIONS (for overlay)
# =============================================================================
survey_xy <- st_coordinates(st_geometry(fit_dsgm$data_sf))
survey_xy_unique <- unique(survey_xy)  # avoid overplotting repeated individuals at same school

# =============================================================================
# 4. THREE-PANEL PLOT
# =============================================================================
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 1))

# --- Panel 1: DSGM prevalence map ---
par(mar = c(3, 3, 3, 5))
plot(r_dsgm, main = "DSGM prevalence", col = col_pal, breaks = breaks,
     range = zlim, axes = TRUE)
contour(r_dsgm, add = TRUE, nlevels = 6, col = "black", labcex = 0.6)
points(survey_xy_unique, pch = 20, cex = 0.3, col = "grey20")

# --- Panel 2: DAST prevalence map ---
par(mar = c(3, 3, 3, 5))
plot(r_dast, main = "DAST prevalence", col = col_pal, breaks = breaks,
     range = zlim, axes = TRUE)
contour(r_dast, add = TRUE, nlevels = 6, col = "black", labcex = 0.6)
points(survey_xy_unique, pch = 20, cex = 0.3, col = "grey20")

# --- Panel 3: DSGM vs DAST scatter with identity line ---
par(mar = c(4, 4, 3, 2))
plot(pred_prev_dsgm$target$prevalence$mean,
     pred_prev_dast$target$prevalence$mean,
     xlab = "DSGM prevalence", ylab = "DAST prevalence",
     main = "DSGM vs DAST", pch = 20, cex = 0.5, col = "steelblue",
     xlim = zlim, ylim = zlim, asp = 1)
abline(0, 1, col = "red", lwd = 1.5)
