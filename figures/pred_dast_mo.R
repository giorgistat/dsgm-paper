rm(list = ls())

library(RiskMap)
library(sf)

# =============================================================================
# 1. LOAD FITTED DAST MODELS (any / light / moderate / high intensity)
# =============================================================================
fits <- list(
  any      = readRDS("fit_sth_prev_dast.rds"),
  light    = readRDS("fit_sth_light_dast.rds"),
  moderate = readRDS("fit_sth_mod_dast.rds"),
  high     = readRDS("fit_sth_high_dast.rds")
)

# =============================================================================
# 2. LOAD DATA (needed to rebuild CRS/UTM and the MDA intervention matrix
#    for the prediction grid; DAST models themselves were fitted on the
#    aggregated location-level data)
# =============================================================================
sth <- read.csv("sth_ken_pi_sw.csv")
sth <- sth[complete.cases(sth[, c("longitude", "latitude", "hkepg")]), ]

sth <- st_as_sf(sth, coords = c("longitude", "latitude"), crs = 4326)
data_utm <- propose_utm(sth)
sth <- st_transform(sth, crs = data_utm)

# Read MDA data
mda_data_sf <- st_read("KenyaCoverage_TotalMDA.shp", quiet = TRUE)
mda_data_sf <- st_transform(mda_data_sf, crs = data_utm)

mda_columns <- grepl("^cov", names(mda_data_sf), ignore.case = TRUE)
cov_columns <- names(mda_data_sf)[mda_columns]
mda_times   <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000

grid_pred <- readRDS("ken_grid_pred.rds")

# =============================================================================
# 3. BUILD INTERVENTION MATRIX ON THE PREDICTION GRID
# =============================================================================
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
# 4. PREDICT PREVALENCE ON THE GRID FOR EACH DAST OUTCOME
# =============================================================================
# time_pred sets the survey time at which MDA impact (via alpha/gamma) is
# evaluated on the grid. Using the most recent survey year here so
# predictions reflect cumulative MDA impact up to the latest data.
time_pred <- 2013

pred_S_dast  <- list()
pred_prev_dast <- list()

for (outcome_name in names(fits)) {

  cat(sprintf("\n=============================\n"))
  cat(sprintf("Predicting DAST: outcome = %s\n", outcome_name))
  cat(sprintf("=============================\n"))

  pred_S_dast[[outcome_name]] <- pred_over_grid(
    fits[[outcome_name]],
    grid_pred = grid_pred,
    type      = "joint"
  )

  pred_prev_dast[[outcome_name]] <- pred_target_grid(
    pred_S_dast[[outcome_name]],
    f_target  = list(prevalence = function(x) 1 / (1 + exp(-x))),
    mda_grid  = intervention_grid,
    time_pred = time_pred
  )
}
# pred_prev_dast$any / light / moderate / high each contain
# $target$prevalence with mean, median, sd, lower (2.5%), upper (97.5%)
# on the grid, from the default pd_summary in pred_target_grid().

saveRDS(pred_S_dast,    "pred_S_dast.rds")
saveRDS(pred_prev_dast, "pred_prev_dast.rds")
