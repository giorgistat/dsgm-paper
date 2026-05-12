# =============================================================================
# TEST SCRIPT FOR dsgm + dast FUNCTIONS (intensity + prevalence)
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
# sth <- sth[sample(1:nrow(sth), 300),]
fit <- readRDS("fit_sth.rds")
sth <- st_as_sf(sth, coords = c("longitude", "latitude"), crs = 4326)
data_utm <- propose_utm(sth)
sth <- st_transform(sth, crs = data_utm)
sth$year <- as.numeric(substr(sth$surveydate, 7, 10))
sth$hkw_pos <- as.integer(sth$hkepg > 0)
sth$x <- st_coordinates(sth)[, 1]
sth$y <- st_coordinates(sth)[, 2]

# =============================================================================
# 2. READ MDA DATA
# =============================================================================
mda_data_sf <- st_read("KenyaCoverage_TotalMDA.shp", quiet = TRUE)
mda_data_sf <- st_transform(mda_data_sf, crs = data_utm)
mda_columns <- grepl("^cov", names(mda_data_sf), ignore.case = TRUE)
cov_columns <- names(mda_data_sf)[mda_columns]
mda_times   <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000

# =============================================================================
# 3. BUILD INTERVENTION MATRIX
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

# =============================================================================
# 4. WARM-START FROM EXISTING FIT
# =============================================================================
par0_vary        <- fit$model_params   # pull beta, k, rho, sigma2, phi, alpha_W, gamma_W
par0_vary$omega1 <- 0.5               # add starting value for the new slope parameter

# =============================================================================
# 5. FIT WITH VARYING OMEGA
# =============================================================================
pen <- penalty_to_dsgm(make_penalty(
  alpha_a    = 27*3, alpha_b = 9*3,
  gamma_type = "lognormal",
  gamma_mean = log(10), gamma_sd = 0.2
))

for(i in 1:10) {
  
  fit_vary_k <- dsgm(
    hkepg ~ gp(sf),
    data         = sth,
    time         = year,
    mda_times    = mda_times,
    int_mat      = intervention,
    penalty      = pen,
    intensity_family = "negbin",
    vary_k       = TRUE,
    omega1_start = 0.5,
    par0         = par0_vary,
    n_samples    = 1000,
    n_warmup     = 1000,
    n_chains     = 1,
    adapt_delta  = 0.9,
    messages     = TRUE
  )
  par0_vary <- coef(fit_vary_k)
  saveRDS(fit_vary_k, paste0("fit_sth_vary_k_",i,".rds"))
}
