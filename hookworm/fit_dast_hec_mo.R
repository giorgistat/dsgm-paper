# =============================================================================
# FIT SCRIPT FOR dast (prevalence only) -- ANY / MODERATE / HIGH INTENSITY
# =============================================================================
rm(list = ls())
set.seed(12345)
library(sf)
library(RiskMap)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
sth <- read.csv("sth_ken_pi_sw.csv")
sth <- sth[complete.cases(sth[, c("longitude", "latitude", "hkepg")]), ]

sth <- st_as_sf(sth, coords = c("longitude", "latitude"), crs = 4326)
data_utm <- propose_utm(sth)
sth <- st_transform(sth, crs = data_utm)
sth$year <- as.numeric(substr(sth$surveydate, 7, 10))

# Dichotomize EPG (double Kato-Katz, values are multiples of 12) into
# four binary outcomes using WHO hookworm intensity cut-offs snapped
# to the nearest achievable multiple of 12:
#   any      : EPG > 0
#   light    : 12   <= EPG <= 1992  (bounded band; WHO 1-1999 cut-off)
#   moderate : 2004 <= EPG <= 3996  (bounded band; WHO 2000-3999 cut-off)
#   high     : EPG >= 4008          (open-ended; first multiple of 12 >= WHO 4000 cut-off)
sth$hkw_pos   <- as.integer(sth$hkepg > 0)
sth$hkw_light <- as.integer(sth$hkepg >= 12   & sth$hkepg <= 1992)
sth$hkw_mod   <- as.integer(sth$hkepg >= 2004 & sth$hkepg <= 3996)
sth$hkw_high  <- as.integer(sth$hkepg >= 4008)

sth$x <- st_coordinates(sth)[, 1]
sth$y <- st_coordinates(sth)[, 2]

# Read MDA data
mda_data_sf <- st_read("KenyaCoverage_TotalMDA.shp", quiet = TRUE)
mda_data_sf <- st_transform(mda_data_sf, crs = data_utm)

mda_columns <- grepl("^cov", names(mda_data_sf), ignore.case = TRUE)
cov_columns <- names(mda_data_sf)[mda_columns]
mda_times   <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000

# =============================================================================
# 2. AGGREGATE OUTCOMES AT LOCATION LEVEL (all three outcomes together,
#    since n_tested is shared across them)
# =============================================================================
sth_df <- st_drop_geometry(sth)

sth_agg <- aggregate(
  cbind(hkw_pos, hkw_light, hkw_mod, hkw_high, hkw_den = 1) ~ x + y + year,
  data = sth_df,
  FUN  = sum
)
names(sth_agg)[names(sth_agg) == "hkw_den"] <- "n_tested"

sth_agg <- st_as_sf(sth_agg, coords = c("x", "y"), crs = data_utm)
sth_agg$x <- st_coordinates(sth_agg)[, 1]
sth_agg$y <- st_coordinates(sth_agg)[, 2]

# =============================================================================
# 3. BUILD INTERVENTION MATRIX ON AGGREGATED LOCATIONS (shared across outcomes)
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

# =============================================================================
# 4. PENALTY (shared across the three intensity models)
# =============================================================================
var_beta <- function(a, b) a * b / ((a + b)^2 * (a + b + 1))

# DAST prevalence penalty
# alpha ~ Beta(16, 4): mean = 0.80, consistent with hookworm cure rates
#                      from Moser et al. 2017 (CR ~0.72-0.80 for albendazole)
# gamma ~ LogNormal(log(7), 0.2): recovery ~7 years
cat("DAST prevalence Beta(16,4): mean =", 16/20,
    "var =", round(var_beta(16, 4), 4), "\n")

pen_prev <- make_penalty(
  alpha_a    = 27*3,
  alpha_b    = 9*3,
  gamma_type = "lognormal",
  gamma_mean = log(10),
  gamma_sd   = 0.2
)

# =============================================================================
# 5. ITERATIVE MCML FIT -- LOOPED OVER OUTCOME (any / moderate / high)
# =============================================================================
par_hat_start <- list(
  beta   = -3.034348,
  sigma2 =  2.034045,
  phi    = 29.486620,
  gamma  =  6.777271,
  alpha  =  0.8
)

outcomes <- list(
  any      = list(col = "hkw_pos",   outfile = "fit_sth_prev_dast.rds"),
  light    = list(col = "hkw_light", outfile = "fit_sth_light_dast.rds"),
  moderate = list(col = "hkw_mod",   outfile = "fit_sth_mod_dast.rds"),
  high     = list(col = "hkw_high",  outfile = "fit_sth_high_dast.rds")
)

fits <- list()

for (outcome_name in names(outcomes)) {

  outcome_col <- outcomes[[outcome_name]]$col
  outfile     <- outcomes[[outcome_name]]$outfile

  cat(sprintf("\n=============================\n"))
  cat(sprintf("Fitting DAST: outcome = %s (%s)\n", outcome_name, outcome_col))
  cat(sprintf("=============================\n"))

  sth_agg$n_pos <- sth_agg[[outcome_col]]
  par_hat_prev  <- par_hat_start

  for (i in 1:3) {
    cat(sprintf("\n--- %s: iteration %d ---\n", outcome_name, i))

    fit_prev <- dast(
      n_pos ~ gp(),
      data        = sth_agg,
      den         = n_tested,
      time        = year,
      mda_times   = mda_times,
      int_mat     = intervention_prev,
      penalty     = penalty_to_dast(pen_prev),
      power_val   = 1,
      scale_to_km = TRUE,
      par0        = par_hat_prev,
      start_pars  = par_hat_prev
    )

    par_hat_prev <- coef(fit_prev)
    cat(sprintf("Current par_hat_prev (%s):\n", outcome_name)); print(par_hat_prev)
  }

  saveRDS(fit_prev, outfile)
  fits[[outcome_name]] <- fit_prev
  cat(sprintf("\nDone: %s. Model saved to %s\n", outcome_name, outfile))
}
