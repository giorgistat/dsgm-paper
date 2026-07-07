# =============================================================================
# Build a prediction grid over the union of southwestern Kenya counties
# =============================================================================
library(sf)
library(rgeoboundaries)

# -----------------------------------------------------------------------------
# 1. Load county boundaries and filter to the study region
# -----------------------------------------------------------------------------
kenya <- geoboundaries(country = "Kenya", adm_lvl = "adm1")

sw_kenya <- c("Bungoma", "Kakamega", "Busia", "Vihiga", "Kericho",
              "Kisumu", "Bomet", "Kisii", "Nyamira",
              "Homa Bay", "Migori")

adm_hit <- kenya[kenya$shapeName %in% sw_kenya, ]

# -----------------------------------------------------------------------------
# 2. Reproject to the same UTM CRS used for model fitting
# -----------------------------------------------------------------------------
adm_hit <- st_transform(adm_hit, crs = data_utm)

# -----------------------------------------------------------------------------
# 3. Take the union of the study-area counties (single polygon)
# -----------------------------------------------------------------------------
study_area <- st_union(adm_hit)

# -----------------------------------------------------------------------------
# 4. Build a regular grid over the bounding box of the study area
# -----------------------------------------------------------------------------
grid_full <- st_make_grid(study_area, cellsize = 2500, what = "centers")
grid_full <- st_as_sf(grid_full)
st_geometry(grid_full) <- "geometry"   # rename away from default "x"

# -----------------------------------------------------------------------------
# 5. Keep only grid points falling inside the study area union
# -----------------------------------------------------------------------------
grid_pred <- st_filter(grid_full, study_area)

# -----------------------------------------------------------------------------
# 6. Add explicit x/y columns (needed by RiskMap's prediction functions)
# -----------------------------------------------------------------------------
coords <- st_coordinates(grid_pred)
grid_pred$x <- coords[, 1]
grid_pred$y <- coords[, 2]

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(grid_pred, "ken_grid_pred.rds")

cat("Grid built:", nrow(grid_pred), "prediction locations within the study area\n")
