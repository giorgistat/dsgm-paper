# =============================================================================
# Motivating examples — two separate figures, each with 3 panels
#
# Figure 1 (STH):  (a) prevalence map  |  (b) intensity map  |  (c) mean-variance
# Figure 2 (LF):   (d) MF map          |  (e) ICT map        |  (f) MF vs ICT scatter
# =============================================================================
library(sf)
library(ggplot2)
library(dplyr)
library(rgeoboundaries)
library(rnaturalearth)
library(patchwork)

fit_sth <- readRDS("fit_sth.rds")
fit_lf  <- readRDS("lf_fit.rds")
sth     <- fit_sth$data_sf
lf      <- fit_lf$data_sf

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_blank(),
        legend.key.size   = unit(0.4, "cm"))

prev_scale_pct <- function(name = "Prevalence\n(%)") {
  scale_colour_gradientn(
    colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
    name    = name
  )
}

# =============================================================================
# Boundaries
# =============================================================================
kenya   <- geoboundaries(country = "Kenya", adm_lvl = "adm1")
kenya   <- st_transform(kenya, st_crs(sth))
adm_hit <- st_filter(kenya, sth)
sw_kenya <- c("Bungoma", "Kakamega", "Busia", "Vihiga", "Kericho", 
              "Kisumu", "Bomet", "Kisii", "Nyamira", 
              "Homa Bay", "Migori")

adm_hit <- adm_hit[adm_hit$shapeName %in% sw_kenya, ]

sf::sf_use_s2(FALSE)
wa_sf <- ne_countries(country = c("Togo", "Benin", "Ghana", "Burkina Faso"),
                      scale = "medium", returnclass = "sf")
wa_sf <- st_transform(wa_sf, st_crs(lf))
sf::sf_use_s2(TRUE)

# =============================================================================
# School-level summaries
# =============================================================================
sth_schools <- sth |>
  st_drop_geometry() |>
  group_by(schoolcode) |>
  summarise(
    prev      = mean(hkepg > 0) * 100,
    intensity = mean(hkepg, na.rm = TRUE),
    lon       = mean(x),
    lat       = mean(y),
    .groups   = "drop"
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sth)) |>
  st_filter(adm_hit)                          # <-- add this line

# =============================================================================
# FIGURE 1 — STH
# =============================================================================

# Panel a — prevalence map
p_a <- ggplot() +
  geom_sf(data = adm_hit,
          fill = "grey92", colour = "grey60", linewidth = 0.25) +
  geom_sf(data = sth_schools,
          aes(colour = prev), size = 1.5, alpha = 0.85) +
  scale_colour_gradientn(
    colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
    name    = "Prevalence\n(%)",
    breaks  = c(0, 10, 25, 50, 75)
  ) +
  labs(title = "(a) Hookworm prevalence") +
  theme_paper +
  theme(axis.title = element_blank(), legend.position = "right")

# Panel b — intensity map
p_b <- ggplot() +
  geom_sf(data = adm_hit,
          fill = "grey92", colour = "grey60", linewidth = 0.25) +
  geom_sf(data = sth_schools,
          aes(colour = intensity + 1), size = 1.5, alpha = 0.85) +
  scale_colour_gradientn(
    colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
    name    = "Mean EPG\n(log scale)",
    trans   = "log10",
    breaks  = c(1, 5, 20, 100, 500),
    labels  = c("0", "4", "19", "99", "499")
  ) +
  labs(title = "(b) Hookworm intensity (mean EPG)") +
  theme_paper +
  theme(axis.title = element_blank(), legend.position = "right")

# Panel c — mean vs variance
mv_df <- sth |>
  st_drop_geometry() |>
  group_by(schoolcode, year) |>
  summarise(mu = mean(hkepg),
            v  = var(hkepg),
            n  = n(),
            .groups = "drop") |>
  filter(n >= 5, mu > 0, v > 0)

p_c <- ggplot(mv_df, aes(x = mu, y = v)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_point(colour = "#2166AC", size = 1.1, alpha = 0.45) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x     = "Mean EPG (log scale)",
       y     = "Variance EPG (log scale)",
       title = "(c) Overdispersion in egg counts") +
  theme_paper

fig_sth <- p_a | p_b | p_c

ggsave("motivating_sth.png",
       plot  = fig_sth,
       width = 14, height = 5, dpi = 300)

ggsave("motivating_sth.pdf",
       plot   = fig_sth,
       width  = 14, height = 5,
       device = cairo_pdf)

cat("STH figure saved\n")

# =============================================================================
# FIGURE 2 — LF
# =============================================================================

# Panel d — MF locations
lf_mf <- lf |> filter(diagnostic == "par")

p_d <- ggplot() +
  geom_sf(data = wa_sf,
          fill = "grey92", colour = "grey60", linewidth = 0.25) +
  geom_sf(data  = lf_mf,
          aes(colour = Prevalence * 100),
          shape = 16, size = 1.8, alpha = 0.85) +
  prev_scale_pct() +
  labs(title = "(d) MF locations") +
  theme_paper +
  theme(axis.title = element_blank(), legend.position = "right")

# Panel e — ICT locations
lf_ict <- lf |> filter(diagnostic == "ser")

p_e <- ggplot() +
  geom_sf(data = wa_sf,
          fill = "grey92", colour = "grey60", linewidth = 0.25) +
  geom_sf(data  = lf_ict,
          aes(colour = Prevalence * 100),
          shape = 17, size = 1.8, alpha = 0.85) +
  prev_scale_pct() +
  labs(title = "(e) ICT locations") +
  theme_paper +
  theme(axis.title = element_blank(), legend.position = "right")

# Panel f — MF vs ICT scatter
lf_wide <- lf |>
  st_drop_geometry() |>
  select(IU_ID, LocationName, Country, diagnostic, Prevalence) |>
  tidyr::pivot_wider(names_from  = diagnostic,
                     values_from = Prevalence,
                     values_fn   = mean) |>
  filter(!is.na(par), !is.na(ser))

country_shapes <- c("Benin"        = 21,
                    "Ghana"        = 24,
                    "Togo"         = 22,
                    "Burkina Faso" = 23)

country_fills  <- c("Benin"        = "#1B7837",
                    "Ghana"        = "#762A83",
                    "Togo"         = "#E66101",
                    "Burkina Faso" = "#2166AC")

p_f <- ggplot(lf_wide,
              aes(x = par * 100, y = ser * 100,
                  fill = Country, shape = Country)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_point(colour = "white", size = 2.8, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values  = country_fills,  name = "Country") +
  scale_shape_manual(values = country_shapes, name = "Country") +
  labs(x     = "MF prevalence (%)",
       y     = "ICT prevalence (%)",
       title = "(f) MF vs ICT prevalence") +
  theme_paper +
  theme(legend.position = "right")

fig_lf <- p_d | p_e | p_f

ggsave("motivating_lf.png",
       plot  = fig_lf,
       width = 14, height = 5, dpi = 300)

ggsave("motivating_lf.pdf",
       plot   = fig_lf,
       width  = 14, height = 5,
       device = cairo_pdf)

cat("LF figure saved\n")