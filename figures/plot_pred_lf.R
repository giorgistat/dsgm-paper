library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)

# Helper to extract the same data.frame that plot.RiskMap_pred_target_grid builds
extract_target_df <- function(pred_target, which_target, which_summary = "mean") {
  coords <- st_coordinates(pred_target$grid_pred)
  val <- pred_target$target[[which_target]][[which_summary]]
  data.frame(x = coords[, "X"], y = coords[, "Y"], value = val)
}

# Extract the three surfaces
df_mf      <- extract_target_df(pred_target, "mf_prevalence", "mean")
df_antigen <- extract_target_df(pred_target, "antigen_prevalence", "mean")
df_worm    <- extract_target_df(pred_target, "worm_burden", "mean")

# Sanity check -- remove once confirmed the two prevalence surfaces are genuinely distinct
stopifnot(!identical(df_mf$value, df_antigen$value))

# Difference surface: antigen minus MF prevalence
df_diff <- df_mf
df_diff$value <- df_antigen$value - df_mf$value

# Common scale for the two prevalence panels
prev_range <- range(c(df_mf$value, df_antigen$value), na.rm = TRUE)

# --- Row 1, panel 1: MF prevalence ---
p_mf <- ggplot(df_mf, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Prevalence", limits = prev_range) +
  coord_equal() +
  labs(title = "MF") +
  theme_minimal() +
  theme(axis.title = element_blank(), plot.title = element_text(face = "bold"))

# --- Row 1, panel 2: Antigen prevalence ---
p_antigen <- ggplot(df_antigen, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Prevalence", limits = prev_range) +
  coord_equal() +
  labs(title = "ICT") +
  theme_minimal() +
  theme(axis.title = element_blank(), plot.title = element_text(face = "bold"))

# --- Row 2, panel 1: ICT - MF ---
diff_lim <- max(abs(df_diff$value), na.rm = TRUE)
p_diff <- ggplot(df_diff, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_distiller(name = "Prev. difference", palette = "RdBu", direction = 1,
                       limits = c(-diff_lim, diff_lim)) +
  coord_equal() +
  labs(title = "ICT - MF") +
  theme_minimal() +
  theme(axis.title = element_blank(), plot.title = element_text(face = "bold"))

# --- Row 2, panel 2: Mean worm burden ---
p_worm <- ggplot(df_worm, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Worm burden", option = "magma") +
  coord_equal() +
  labs(title = "Mean worm burden") +
  theme_minimal() +
  theme(axis.title = element_blank(), plot.title = element_text(face = "bold"))

# --- Combine: 2x2 grid ---
final_plot <- (p_mf | p_antigen) / (p_diff | p_worm)

pdf("pred_maps_lf.pdf", width = 8, height = 8)
final_plot
dev.off()

