library(MASS)
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

sth$hkw_pos <- as.integer(sth$hkepg > 0)

sth$x <- st_coordinates(sth)[, 1]
sth$y <- st_coordinates(sth)[, 2]

sth_df <- st_drop_geometry(sth)

# Compute group-level stats and fit NB to each school-year
school_year_groups <- split(sth_df, list(sth_df$x, sth_df$y, sth_df$year), drop = TRUE)

results <- lapply(school_year_groups, function(g) {
  epg <- g$hkepg
  n   <- length(epg)
  mn  <- mean(epg)
  pv  <- mean(epg > 0)

  # Need enough observations and at least some positives and some zeros
  # to identify both mean and overdispersion
  if (n < 10 || sum(epg > 0) < 2 || sum(epg == 0) < 2) return(NULL)

  fit <- tryCatch(
    fitdistr(epg, "negative binomial"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  data.frame(
    n       = n,
    mean    = mn,
    prev    = pv,
    omega   = fit$estimate["size"],
    omega_se = fit$sd["size"]
  )
})

results_df <- do.call(rbind, Filter(Negate(is.null), results))

# Plot
pdf("agg_vs_mean.pdf", width = 8, height = 8)
plot(log(results_df$mean + 1), log(results_df$omega),
     xlab = "mean eggs per gram (log scale)",
     ylab = expression("Aggregation parameter" ~omega~" (log scale)"),
     pch  = 20,
     cex  = sqrt(results_df$n) / max(sqrt(results_df$n)) * 2)
abline(lm(log(omega) ~ log(mean + 1), data = results_df), col = "red")
dev.off()

# How many groups in total
length(school_year_groups)

# Names of the groups (location_year combinations)
head(names(school_year_groups))

# Look at one specific group
school_year_groups[[1]]

# Or by name if you know the location-year
school_year_groups[["123456.78_987654.32_2015"]]

# See the size of each group
group_sizes <- sapply(school_year_groups, nrow)
table(group_sizes)  # distribution of group sizes
sort(group_sizes)   # ordered from smallest to largest
