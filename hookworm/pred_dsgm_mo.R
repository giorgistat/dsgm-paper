rm(list = ls())

library(RiskMap)
library(sf)
fit_dsgm <- readRDS("fit_sth_dsgm.rds")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
sth <- read.csv("sth_ken_pi_sw.csv")
sth <- sth[complete.cases(sth[, c("longitude", "latitude", "hkepg")]), ]

sth <- st_as_sf(sth, coords = c("longitude", "latitude"), crs = 4326)
data_utm <- propose_utm(sth)
sth <- st_transform(sth, crs = data_utm)
sth$year <- as.numeric(substr(sth$surveydate, 7, 10))

# Add planar coordinates (needed for intervention matrix)
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
# 2. BUILD INTERVENTION MATRIX FOR THE ORIGINAL INTENSITY DATA
# =============================================================================
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

#pred_S_dsgm <- pred_over_grid(fit_dsgm, grid_pred = grid_pred)
#saveRDS(pred_S_dsgm, "pred_S_dsgm.rds")
pred_S_dsgm <- readRDS("pred_S_dsgm.rds")
# =============================================================================
# 3. BUILD INTERVENTION MATRIX ON THE PREDICTION GRID
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
# 4. PREDICT PREVALENCE ON THE GRID (DSGM)
# =============================================================================
# time_pred sets the survey time at which MDA impact (via alpha_W/gamma_W) is
# evaluated on the grid. Using the most recent survey year here so predictions
# reflect cumulative MDA impact up to the latest data.
time_pred <- 2013

pred_prev_dsgm <- pred_target_grid(
  pred_S_dsgm,
  mda_grid  = intervention_grid,
  time_pred = time_pred
)

saveRDS(pred_prev_dsgm, "pred_prev_dsgm.rds")

# =============================================================================
# 5. LIGHT / MODERATE / HIGH INTENSITY PREVALENCE (DSGM)
# =============================================================================
# WHO hookworm intensity bands on the double-KK EPG scale (multiples of 12):
#   light    : 12   <= EPG <= 1992
#   moderate : 2004 <= EPG <= 3996
#   high     : EPG >= 4008
# Closed-form banded intensity probability P(a <= Y <= b | S(x)), derived from
# the shifted Negative Binomial approximation to the positive-count
# distribution (model formulation, Section on approximation of positive
# counts). Uses the same k(x)/rho parametrisation as the built-in
# 'prevalence'/'intensity' targets for the sth (intprev) family.

k_val   <- pred_S_dsgm$par_hat$k
rho_val <- pred_S_dsgm$par_hat$rho
c_rho   <- 1 - exp(-rho_val)
vary_k_flag <- isTRUE(pred_S_dsgm$vary_k) && !is.null(pred_S_dsgm$par_hat$omega1)
omega1_val  <- if (vary_k_flag) pred_S_dsgm$par_hat$omega1 else 0

k_fun <- function(mu_W) {
  if (vary_k_flag) exp(log(k_val) + omega1_val * log(pmax(mu_W, 1e-10)))
  else array(k_val, dim = dim(mu_W))
}

banded_intensity_prob <- function(a, b) {
  force(a); force(b)
  function(lp) {
    mu_W <- exp(lp)
    k_i  <- k_fun(mu_W)
    p    <- 1 - (k_i / (k_i + mu_W * c_rho))^k_i
    p    <- pmin(pmax(p, 1e-10), 1 - 1e-10)

    mu_C     <- rho_val * mu_W / p
    sigma2_C <- rho_val * mu_W * (1 + rho_val) / p +
      (rho_val^2 * mu_W^2 / p) * (1 / k_i + 1 - 1 / p)
    mu_Z <- mu_C - 1
    r    <- mu_Z^2 / (sigma2_C - mu_Z)
    p_z  <- r / (r + mu_Z)

    lower_arg <- max(a, 1) - 1
    Fc_lower  <- if (lower_arg <= 0) 0 else pnbinom(lower_arg - 1, size = r, prob = p_z)
    Fc_upper  <- if (is.infinite(b)) 1 else pnbinom(b - 1, size = r, prob = p_z)

    out <- p * (Fc_upper - Fc_lower)
    if (a <= 0) out <- out + (1 - p)
    pmin(pmax(out, 0), 1)
  }
}

f_target_intensity <- list(
  light_prev    = banded_intensity_prob(12,   1992),
  moderate_prev = banded_intensity_prob(2004, 3996),
  high_prev     = banded_intensity_prob(4008, Inf)
)

pred_intensity_dsgm <- pred_target_grid(
  pred_S_dsgm,
  f_target  = f_target_intensity,
  mda_grid  = intervention_grid,
  time_pred = time_pred
)
# pred_intensity_dsgm$target$light_prev / moderate_prev / high_prev each
# contain: mean, median, sd, lower (2.5%), upper (97.5%) on the grid,
# from the default pd_summary in pred_target_grid().

saveRDS(pred_intensity_dsgm, "pred_intensity_dsgm.rds")
