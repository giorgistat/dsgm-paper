rm(list = ls())

library(RiskMap)

fit_dsgm <- readRDS("fit_sth_dsgm.rds")

anpit_hk <-
  assess_pp(list(DSGM_hk = fit_dsgm),
            method = "regularized",
            n_size = 29,
            min_dist = 5,
            iter = 1)

saveRDS(anpit_hk, file  = "anpit_hk_one.rds")

anpit_hk_mult <-
  assess_pp(list(DSGM_hk = fit_dsgm),
            method = "regularized",
            n_size = 29,
            min_dist = 5,
            iter = 10)

saveRDS(anpit_hk_mult, file  = "anpit_hk_mult.rds")
