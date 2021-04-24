# clearing all loaded packages to prevent conflicts
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T, unload = T,force = T))

library(blavaan)
library(tidyverse)

bsem.summary <- function(fit) {
  out <- trimws(summary(fit, neff = TRUE, fit.measures = TRUE))
  return(out)
}

set.seed(2019)
options(mc.cores = parallel::detectCores()) ## Run chains in parallel
blavaan_data <- readRDS("data/blavaan_data.rds")

# H3 ---------------------------
## Does dynamic norm (versus static or no norm) information increase participants’ positive attitude, intentions, and expectations to reduce their meat consumption?
simple.mod <- '
INTEREST       ~ conditionbi
attitude_mean  ~ conditionbi
expintent_avg  ~ conditionbi'

### Fitting models
simpleuninf.fit <- bsem(model = simple.mod, data = blavaan_data) # uninformative
simplelowinf.fit <- bsem(model = simple.mod, data = blavaan_data, dp = dpriors(beta = "normal(0.5,0.75)")) # informative
simplehighinf.fit <- bsem(model = simple.mod, data = blavaan_data, dp = dpriors(beta = "normal(0.5, 0.35)")) # informative

### Fit measures
h3.fit <- list(simpleuninf.fit, simplelowinf.fit, simplehighinf.fit)
h3.global <- lapply(h3.fit, fitmeasures)

### Summary outputs
h3.out <- list(lowinf = simplelowinf.fit,
               highinf = simplehighinf.fit,
               uninf = simpleuninf.fit) %>%
  lapply(., bsem.summary)

### Estimate bias
h3.estimates <- lapply(h3.out, function(x) as.numeric(x[, 2][1:3]))
h3.bias <- c(lowinf_uninf = round(100 * ((h3.estimates$lowinf - h3.estimates$uninf) / h3.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h3.estimates$highinf - h3.estimates$uninf) / h3.estimates$uninf), 2))

### Summary table
h3.table <- cbind("Parameter" = rbind("Interest", "Attitude", "Intention/Expectation"),
                  "Mean (SD)" = paste0(h3.out$uninf[, 2], " (", h3.out$uninf[, 3], ")"),
                  "95% PPI" = paste(h3.out$uninf[, 4], h3.out$uninf[, 5], sep = ", "),
                  "neff" = h3.out$uninf[, 7],
                  "PSRF" = round(blavInspect(simpleuninf.fit, "psrf"), 3),
                  "Prior" = h3.out$uninf[, 8],
                  "Mean (SD)" = paste0(h3.out$lowinf[, 2], " (", h3.out$lowinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$lowinf[, 4], h3.out$lowinf[, 5], sep = ", "),
                  "neff" = h3.out$lowinf[, 7],
                  "PSRF" = round(blavInspect(simplelowinf.fit, "psrf"), 3),
                  "Bias (%)" = h3.bias[1:3],
                  "Prior" = h3.out$lowinf[, 8],
                  "Mean (SD)" = paste0(h3.out$highinf[, 2], " (", h3.out$highinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$highinf[, 4], h3.out$highinf[, 5], sep = ", "),
                  "neff" = h3.out$highinf[, 7],
                  "PSRF" = round(blavInspect(simplehighinf.fit, "psrf"), 3),
                  "Bias (%)" = h3.bias[4:6],
                  "Prior" = h3.out$highinf[, 8])

# H4 ---------------------------
full.mod <- '
INTEREST       ~ conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ conditionbi + AGE + genderbi + POLITICS
expintent_avg ~ conditionbi + AGE + genderbi + POLITICS'

full_low.mod <- '
INTEREST       ~ prior("normal(0.5,0.75)")*conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ prior("normal(0.5,0.75)")*conditionbi + AGE + genderbi + POLITICS
expintent_avg ~ prior("normal(0.5,0.75)")*conditionbi + AGE + genderbi + POLITICS'

full_high.mod <- '
INTEREST       ~ prior("normal(0.5,0.35)")*conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ prior("normal(0.5,0.35)")*conditionbi + AGE + genderbi + POLITICS
expintent_avg ~ prior("normal(0.5,0.35)")*conditionbi + AGE + genderbi + POLITICS'

### Fitting models
full.fit <- bsem(full.mod, data = blavaan_data, target = "stan", seed = 2019) # default
full_lowinf.fit <- bsem(model = full_low.mod, data = blavaan_data, target = "stan", seed = 2019) # informative
full_highinf.fit <- bsem(model = full_high.mod, data = blavaan_data, target = "stan", seed = 2019) # informative
h4.fit <- list(full.fit, full_lowinf.fit, full_highinf.fit)
h4.global <- lapply(h4.fit, fitmeasures)


### Summary outputs
h4.out <- list(lowinf = full_lowinf.fit,
               highinf = full_highinf.fit,
               uninf = full.fit) %>%
  lapply(., bsem.summary)

### Estimate bias
h4.estimates <- lapply(h4.out, function(x) as.numeric(x[, 2][1:12]))
h4.bias <- c(lowinf_uninf = round(100 * ((h4.estimates$lowinf - h4.estimates$uninf) / h4.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h4.estimates$highinf - h4.estimates$uninf) / h4.estimates$uninf), 2))

### Summary table
h4.table <- cbind("Parameter" = rep(c("Condition", "Age", "Gender", "Politics"), 4),
                  "Mean (SD)" = paste0(h4.out$uninf[1:12, 2], " (", h4.out$uninf[1:12, 3], ")"),
                  "95% PPI" = paste(h4.out$uninf[1:12, 4], h4.out$uninf[1:12, 5], sep = ", "),
                  "neff" = h4.out$uninf[1:12, 7],
                  "PSRF" = round(blavInspect(full.fit, "psrf"), 3) %>% .[1:12],
                  "Prior" = h4.out$uninf[1:12, 8],
                  "Mean (SD)" = paste0(h4.out$lowinf[1:12, 2], " (", h4.out$lowinf[1:12, 3], ")"),
                  "95% PPI" = paste(h4.out$lowinf[1:12, 4], h4.out$lowinf[1:12, 5], sep = ", "),
                  "neff" = h4.out$lowinf[1:12, 7],
                  "PSRF" = round(blavInspect(full_lowinf.fit, "psrf"), 3) %>% .[1:12],
                  "Bias (%)" = h4.bias[1:12],
                  "Prior" = h4.out$lowinf[1:12, 8],
                  "Mean (SD)" = paste0(h4.out$highinf[1:16, 2], " (", h4.out$highinf[1:12, 3], ")"),
                  "95% PPI" = paste(h4.out$highinf[1:12, 4], h4.out$highinf[1:12, 5], sep = ", "),
                  "neff" = h4.out$highinf[1:12, 7],
                  "PSRF" = round(blavInspect(full_highinf.fit, "psrf"), 3) %>% .[1:12],
                  "Bias (%)" = h4.bias[13:24],
                  "Prior" = h4.out$highinf[1:12, 8])


# H5 ---------------------------
interact.fit <- bsem(simple.mod, data = blavaan_data, target = "stan", group = "agebi", seed = 2019) # uninformative
interact_lowinf.fit <- bsem(model = simple.mod, data = blavaan_data, target = "stan", group = "agebi", seed = 2019, dp = dpriors(beta = "normal(0.5,0.75)")) # informative
interact_highinf.fit <- bsem(model = simple.mod, data = blavaan_data, target = "stan", group = "agebi", seed = 2019, dp = dpriors(beta = "normal(0.5,0.35)")) # highly informative
h5.fit <- list(interact.fit, interact_lowinf.fit, interact_highinf.fit)
h5.global <- lapply(h5.fit, fitmeasures)

### Summary outputs
h5.out <- list(uninf = interact.fit,
               lowinf = interact_lowinf.fit,
               highinf = interact_highinf.fit) %>%
  lapply(., bsem.summary)

### Estimate bias
h5.estimates <- lapply(h5.out, function(x) as.numeric(x[, 2][c(1:3, 15:17)]))
h5.bias <- c(lowinf_uninf = round(100 * ((h5.estimates$lowinf - h5.estimates$uninf) / h5.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h5.estimates$highinf - h5.estimates$uninf) / h5.estimates$uninf), 2))

interact.table <- cbind("Parameter" = rep(c("Interest", "Attitude", "Intention", "Expectation"),2),
                        "Mean (SD)" = paste0(h5.out$uninf[c(1:3, 15:17), 2], " (", h5.out$uninf[c(1:3, 15:17), 3], ")"),
                        "95% PPI" = paste(h5.out$uninf[c(1:3, 15:17), 4], h5.out$uninf[c(1:3, 15:17), 5], sep = ", "),
                        "neff" = h5.out$uninf[c(1:3, 15:17), 7],
                        "PSRF" = round(blavInspect(interact.fit, "psrf"), 3) %>% .[c(1:3,13:15)],
                        "Prior" = h5.out$uninf[c(1:3, 15:17), 8],
                        "Mean (SD)" = paste0(h5.out$lowinf[c(1:3, 15:17), 2], " (", h5.out$lowinf[c(1:3, 15:17), 3], ")"),
                        "95% PPI" = paste(h5.out$lowinf[c(1:3, 15:17), 4], h5.out$lowinf[c(1:3, 15:17), 5], sep = ", "),
                        "neff" = h5.out$lowinf[c(1:3, 15:17), 7],
                        "PSRF" = round(blavInspect(interact_lowinf.fit, "psrf"), 3) %>% .[c(1:3,13:15)],
                        "Bias (%)" = h5.bias[1:6],
                        "Prior" = h5.out$lowinf[c(1:3, 15:17), 8],
                        "Mean (SD)" = paste0(h5.out$highinf[c(1:3, 15:17), 2], " (", h5.out$highinf[c(1:3, 15:17), 3], ")"),
                        "95% PPI" = paste(h5.out$highinf[c(1:3, 15:17), 4], h5.out$highinf[c(1:3, 15:17), 5], sep = ", "),
                        "neff" = h5.out$highinf[c(1:3, 15:17), 7],
                        "PSRF" = round(blavInspect(interact_highinf.fit, "psrf"), 3) %>% .[c(1:3,13:15)],
                        "Bias (%)" = h5.bias[7:12],
                        "Prior" = h5.out$highinf[c(1:3, 15:17), 8])

save(h3.fit, file =  "data/blav_fit/h3_fit.RData")
save(h4.fit, file =  "data/blav_fit/h4_fit.RData")
save(h5.fit, file =  "data/blav_fit/h5_fit.RData")

save(h3.global, h4.global, h5.global, h3.table, h4.table, interact.table, file =  "data/exp_results.RData")

