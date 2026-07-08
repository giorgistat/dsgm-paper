rm(list = ls())

library(RiskMap)
library(sf)
fit_dsgm <- readRDS("fit_sth_dsgm.rds")
fit_dast <- readRDS("fit_sth_prev_dast.rds")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
sth <- read.csv("sth_ken_pi_sw.csv")
sth <- sth[complete.cases(sth[, c("longitude", "latitude", "hkepg")]), ]

sth <- st_as_sf(sth, coords = c("longitude", "latitude"), crs = 4326)
data_utm <- propose_utm(sth)
sth <- st_transform(sth, crs = data_utm)
sth$year <- as.numeric(substr(sth$surveydate, 7, 10))

# Dichotomise intensity outcome
sth$hkw_pos <- as.integer(sth$hkepg > 0)

# Add planar coordinates (needed for aggregation key and intervention matrix)
sth$x <- st_coordinates(sth)[, 1]
sth$y <- st_coordinates(sth)[, 2]

# Read MDA data
mda_data_sf <- st_read("KenyaCoverage_TotalMDA.shp", quiet = TRUE)
mda_data_sf <- st_transform(mda_data_sf, crs = data_utm)

# Extract MDA times
mda_columns <- grepl("^cov", names(mda_data_sf), ignore.case = TRUE)
cov_columns <- names(mda_data_sf)[mda_columns]
mda_times   <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000

# =============================================================================
# 2. AGGREGATE PREVALENCE DATA AT LOCATION LEVEL
# =============================================================================
sth_df <- st_drop_geometry(sth)

sth_agg <- aggregate(
  cbind(hkw_pos, hkw_den = 1) ~ x + y + year,
  data = sth_df,
  FUN  = sum
)
names(sth_agg)[names(sth_agg) == "hkw_pos"] <- "n_pos"
names(sth_agg)[names(sth_agg) == "hkw_den"] <- "n_tested"

sth_agg <- st_as_sf(sth_agg, coords = c("x", "y"), crs = data_utm)
sth_agg$x <- st_coordinates(sth_agg)[, 1]
sth_agg$y <- st_coordinates(sth_agg)[, 2]

# =============================================================================
# 3. BUILD INTERVENTION MATRIX ON AGGREGATED LOCATIONS
# =============================================================================
unique_coords  <- unique(data.frame(x = sth_agg$x, y = sth_agg$y))
unique_locs_sf <- st_as_sf(unique_coords, coords = c("x", "y"), crs = data_utm)
spatial_join   <- st_join(unique_locs_sf, mda_data_sf, join = st_intersects, left = TRUE)

spatial_join_df     <- st_drop_geometry(spatial_join)
cov_matrix          <- as.matrix(spatial_join_df[, cov_columns, drop = FALSE])
intervention_unique <- ifelse(cov_matrix > 0, 1L, 0L)
intervention_unique[is.na(intervention_unique)] <- 0L

unique_coords$key <- paste(unique_coords$x, unique_coords$y, sep = "_")
sth_agg$key       <- paste(sth_agg$x, sth_agg$y, sep = "_")
matched_idx       <- match(sth_agg$key, unique_coords$key)
intervention_prev <- intervention_unique[matched_idx, , drop = FALSE]
sth_agg$key       <- NULL

# Also build intervention matrix for the original intensity data
unique_coords_nb  <- unique(data.frame(x = sth$x, y = sth$y))
unique_locs_nb_sf <- st_as_sf(unique_coords_nb, coords = c("x", "y"), crs = data_utm)
spatial_join_nb   <- st_join(unique_locs_nb_sf, mda_data_sf, join = st_intersects, left = TRUE)

spatial_join_nb_df     <- st_drop_geometry(spatial_join_nb)
cov_matrix_nb          <- as.matrix(spatial_join_nb_df[, cov_columns, drop = FALSE])
intervention_unique_nb <- ifelse(cov_matrix_nb > 0, 1L, 0L)
intervention_unique_nb[is.na(intervention_unique_nb)] <- 0L

unique_coords_nb$key <- paste(unique_coords_nb$x, unique_coords_nb$y, sep = "_")
sth$key              <- paste(sth$x, sth$y, sep = "_")
matched_idx_nb       <- match(sth$key, unique_coords_nb$key)
intervention         <- intervention_unique_nb[matched_idx_nb, , drop = FALSE]
sth$key              <- NULL

grid_pred <- readRDS("ken_grid_pred.rds")

pred_S_dast <- pred_over_grid(fit_dast, grid_pred = grid_pred, type = "joint")
saveRDS(pred_S_dast, "pred_S_dast.rds")

pred_S_dsgm <- pred_over_grid(fit_dsgm, grid_pred = grid_pred)
saveRDS(pred_S_dsgm, "pred_S_dsgm.rds")

# =============================================================================
# 4. BUILD INTERVENTION MATRIX ON THE PREDICTION GRID
# =============================================================================
# grid_pred must contain coordinates in the same CRS/units as sth (UTM, scaled
# to km if scale_to_km = TRUE was used in fitting). Adjust the coordinate
# extraction below if grid_pred stores x/y differently (e.g. as an sf object
# vs a plain data.frame with x, y columns).

if (inherits(grid_pred, "sf")) {
  grid_coords_df <- st_drop_geometry(grid_pred)
  grid_sf        <- grid_pred
} else {
  grid_coords_df <- grid_pred
  grid_sf <- st_as_sf(grid_pred, coords = c("x", "y"), crs = data_utm)
}

grid_join    <- st_join(grid_sf, mda_data_sf, join = st_intersects, left = TRUE)
grid_join_df <- st_drop_geometry(grid_join)

grid_cov_matrix   <- as.matrix(grid_join_df[, cov_columns, drop = FALSE])
intervention_grid <- ifelse(grid_cov_matrix > 0, 1L, 0L)
intervention_grid[is.na(intervention_grid)] <- 0L

# =============================================================================
# 5. PREDICT TARGETS ON THE GRID FOR BOTH MODELS
# =============================================================================
# time_pred sets the survey time at which MDA impact (via alpha_W/gamma_W or
# alpha/gamma) is evaluated on the grid. Using the most recent survey year
# here so predictions reflect cumulative MDA impact up to the latest data.
time_pred <- 2013

pred_prev_dast <- pred_target_grid(
  pred_S_dast,
  f_target = list(prevalence = function(x) 1/(1+exp(-x))),
  mda_grid  = intervention_grid,
  time_pred = time_pred
)

pred_prev_dsgm <- pred_target_grid(
  pred_S_dsgm,
  mda_grid  = intervention_grid,
  time_pred = time_pred
)

saveRDS(pred_prev_dast, "pred_prev_dast.rds")
saveRDS(pred_prev_dsgm, "pred_prev_dsgm.rds")

