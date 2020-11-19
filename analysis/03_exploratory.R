# Load packages ---------------------------
list.of.packages <- c("ggplot2", "Rcpp", "papaja", "tidyverse", "haven", "ggpubr",
                      "psy", "knitr", "kableExtra", "rlang", "blavaan")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

options(mc.cores = parallel::detectCores()) ## Run chains in parallel

clean <- haven::read_spss("../data/02_clean_data.sav")

clean[, c("GENDER", "condition", "conditionbi", "genderbi", "agebi")] <- haven::as_factor(clean[, c("GENDER", "condition", "conditionbi", "genderbi", "agebi")])

# Create functions ---------------------------
bsem.plot <- function(data, column) {
  ggplot(data = data[[column]],
         mapping = aes(x        = value,
                       fill     = id1,
                       colour   = id2,
                       linetype = id2,
                       alpha    = id2)) +
    geom_density(size = 1) +
    scale_x_continuous(limits = c(0, 2)) +
    scale_colour_manual(name = 'Posterior/Prior', values = c("black", "red"), labels = c("posterior", "prior")) +
    scale_linetype_manual(name = 'Posterior/Prior', values = c("solid", "dotted"), labels = c("posterior", "prior")) +
    scale_alpha_discrete(name = 'Posterior/Prior', range = c(.7, .3), labels = c("posterior", "prior")) +
    scale_fill_manual(name = "Densities", values = c("Yellow", "darkred", "blue", "green", "pink")) +
    labs(title = bquote(~ beta[.(as.name(gsub("^(.)", "\\U\\1", column, perl = T)))]))
}
bsem.summary <- function(fit) {
  out <- trimws(summary(fit, neff = TRUE, fit.measures = TRUE))
  return(out)
}

# H3 ---------------------------
## Does dynamic norm (versus static or no norm) information increase participants’ positive attitude, intentions, and expectations to reduce their meat consumption?
simple.mod <- '
INTEREST       ~ conditionbi
attitude_mean  ~ conditionbi
intention_mean ~ conditionbi
expect_mean    ~ conditionbi'

### Fitting models
simpleuninf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019) # uninformative
simplelowinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019, dp = dpriors(beta = "normal(0.5,0.447)")) # informative
simplehighinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", seed = 2019, dp = dpriors(beta = "normal(0.5, 0.316)")) # informative

### Summary outputs
h3.out <- list(lowinf = simplelowinf.fit,
               highinf = simplehighinf.fit,
               uninf = simpleuninf.fit) %>%
  lapply(., bsem.summary)

h3.table <- cbind("Parameter" = rbind("Interest", "Attitude", "Intention", "Expectation"),
                  "Mean (SD)" = paste0(h3.out$uninf[, 2], " (", h3.out$uninf[, 3], ")"),
                  "95% PPI" = paste(h3.out$uninf[, 4], h3.out$uninf[, 5], sep = ", "),
                  "neff" = h3.out$uninf[, 7],
                  "PSRF" = round(blavInspect(simpleuninf.fit, "psrf"), 3),
                  "Prior" = h3.out$uninf[, 8],
                  "Mean (SD)" = paste0(h3.out$lowinf[, 2], " (", h3.out$lowinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$lowinf[, 4], h3.out$lowinf[, 5], sep = ", "),
                  "neff" = h3.out$lowinf[, 7],
                  "PSRF" = round(blavInspect(simplelowinf.fit, "psrf"), 3),
                  "Prior" = h3.out$lowinf[, 8],
                  "Mean (SD)" = paste0(h3.out$highinf[, 2], " (", h3.out$highinf[, 3], ")"),
                  "95% PPI" = paste(h3.out$highinf[, 4], h3.out$highinf[, 5], sep = ", "),
                  "neff" = h3.out$highinf[, 7],
                  "PSRF" = round(blavInspect(simplehighinf.fit, "psrf"), 3),
                  "Prior" = h3.out$highinf[, 8])

### Estimate bias
h3.estimates <- lapply(h3.out, function(x) as.numeric(x[, 2][1:4]))
h3.bias <- c(lowinf_uninf = round(100 * ((h3.estimates$lowinf - h3.estimates$uninf) / h3.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h3.estimates$highinf - h3.estimates$uninf) / h3.estimates$uninf), 2))

posterior1.2.3 <- list(interest = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[1]']),
                                            .id = "id1"),
                       attitude = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[2]']),
                                            .id = "id1"),
                       intention = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[3]']),
                                             .id = "id1"),
                       expectation = bind_rows("uninformative prior"      = enframe(as.matrix(blavInspect(simpleuninf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               "informative prior (low)"  = enframe(as.matrix(blavInspect(simplelowinf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               "informative prior (high)" = enframe(as.matrix(blavInspect(simplehighinf.fit, what = "mcmc"))[, 'bet_sign[4]']),
                                               .id = "id1"))

prior1.2.3 <- bind_rows("uninformative prior" = enframe(rnorm(10000, mean = 0, sd = sqrt(1 / 1e-2))),
                        "informative prior (low)" = enframe(rnorm(10000, mean = 0.5, sd = sqrt(0.2))),
                        "informative prior (high)" = enframe(rnorm(10000, mean = 0.5, sd = sqrt(0.01))),
                        .id = "id1") # here we sample a large number of values from the prior distributions to be able to plot them.

h3priors.posterior <- list(interest = bind_rows("posterior" = posterior1.2.3$interest, "prior" = prior1.2.3, .id = "id2"),
                           attitude = bind_rows("posterior" = posterior1.2.3$attitude, "prior" = prior1.2.3, .id = "id2"),
                           intention = bind_rows("posterior" = posterior1.2.3$intention, "prior" = prior1.2.3, .id = "id2"),
                           expectation = bind_rows("posterior" = posterior1.2.3$expectation, "prior" = prior1.2.3, .id = "id2"))

h3.plots <- lapply(seq_along(h3priors.posterior), function(x) bsem.plot(h3priors.posterior, names(h3priors.posterior)[x]))
ggarrange(plotlist = h3.plots, common.legend = TRUE, legend = "bottom", labels = c("A", "B", "C", "D"))

# H4 ---------------------------
interact.fit <- bsem(simple.mod, data = clean, target = "stan", group = "agebi", seed = 2019) # uninformative
interact_lowinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", group = "agebi", seed = 2019, dp = dpriors(beta = "normal(0.5,0.447)")) # informative
interact_highinf.fit <- bsem(model = simple.mod, data = clean, target = "stan", group = "agebi", seed = 2019, dp = dpriors(beta = "normal(0.5,0.316)")) # highly informative

h4.out <- list(uninf = interact.fit,
               lowinf = interact_lowinf.fit,
               highinf = interact_highinf.fit) %>%
  lapply(., bsem.summary)

### Estimate bias
h4.estimates <- lapply(h4.out, function(x) as.numeric(x[, 2][c(1:4, 21:24)]))
h4.bias <- c(lowinf_uninf = round(100 * ((h4.estimates$lowinf - h4.estimates$uninf) / h4.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h4.estimates$highinf - h4.estimates$uninf) / h4.estimates$uninf), 2))

interact.table <- cbind("Parameter" = rep(c("Interest", "Attitude", "Intention", "Expectation"),2),
                        "Mean (SD)" = paste0(h4.out$uninf[c(1:4,21:24), 2], " (", h4.out$uninf[c(1:4,21:24), 3], ")"),
                        "95% PPI" = paste(h4.out$uninf[c(1:4,21:24), 4], h4.out$uninf[c(1:4,21:24), 5], sep = ", "),
                        "neff" = h4.out$uninf[c(1:4,21:24), 7],
                        "PSRF" = round(blavInspect(interact.fit, "psrf"), 3) %>% .[c(1:4,19:22)],
                        "Prior" = h4.out$uninf[c(1:4,21:24), 8],
                        "Mean (SD)" = paste0(h4.out$lowinf[c(1:4,21:24), 2], " (", h4.out$lowinf[c(1:4,21:24), 3], ")"),
                        "95% PPI" = paste(h4.out$lowinf[c(1:4,21:24), 4], h4.out$lowinf[c(1:4,21:24), 5], sep = ", "),
                        "neff" = h4.out$lowinf[c(1:4,21:24), 7],
                        "PSRF" = round(blavInspect(interact_lowinf.fit, "psrf"), 3) %>% .[c(1:4,19:22)],
                        "Prior" = h4.out$lowinf[c(1:4,21:24), 8],
                        "Mean (SD)" = paste0(h4.out$highinf[c(1:4,21:24), 2], " (", h4.out$highinf[c(1:4,21:24), 3], ")"),
                        "95% PPI" = paste(h4.out$highinf[c(1:4,21:24), 4], h4.out$highinf[c(1:4,21:24), 5], sep = ", "),
                        "neff" = h4.out$highinf[c(1:4,21:24), 7],
                        "PSRF" = round(blavInspect(interact_highinf.fit, "psrf"), 3) %>% .[c(1:4,19:22)],
                        "Prior" = h4.out$highinf[c(1:4,21:24), 8])

# H5 ---------------------------
full.mod <- '
INTEREST       ~ conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ conditionbi + AGE + genderbi + POLITICS
intention_mean ~ conditionbi + AGE + genderbi + POLITICS
expect_mean    ~ conditionbi + AGE + genderbi + POLITICS'

full_low.mod <- '
INTEREST       ~ prior("normal(0.5,0.447)")*conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ prior("normal(0.5,0.447)")*conditionbi + AGE + genderbi + POLITICS
intention_mean ~ prior("normal(0.5,0.447)")*conditionbi + AGE + genderbi + POLITICS
expect_mean    ~ prior("normal(0.5,0.447)")*conditionbi + AGE + genderbi + POLITICS'

full_high.mod <- '
INTEREST       ~ prior("normal(0.5,0.316)")*conditionbi + AGE + genderbi + POLITICS
attitude_mean  ~ prior("normal(0.5,0.316)")*conditionbi + AGE + genderbi + POLITICS
intention_mean ~ prior("normal(0.5,0.316)")*conditionbi + AGE + genderbi + POLITICS
expect_mean    ~ prior("normal(0.5,0.316)")*conditionbi + AGE + genderbi + POLITICS'

### Fitting models
full.fit <- bsem(full.mod, data = clean, target = "stan", seed = 2019) # default
full_lowinf.fit <- bsem(model = full_low.mod, data = clean, target = "stan", seed = 2019) # informative
full_highinf.fit <- bsem(model = full_high.mod, data = clean, target = "stan", seed = 2019) # informative

### Summary outputs
h5.out <- list(lowinf = full_lowinf.fit,
               highinf = full_highinf.fit,
               uninf = full.fit) %>%
  lapply(., bsem.summary)

h5.table <- cbind("Parameter" = rep(c("Condition", "Age", "Gender", "Politics"), 4),
                  "Mean (SD)" = paste0(h5.out$uninf[1:16, 2], " (", h5.out$uninf[1:16, 3], ")"),
                  "95% PPI" = paste(h5.out$uninf[1:16, 4], h5.out$uninf[1:16, 5], sep = ", "),
                  "neff" = h5.out$uninf[1:16, 7],
                  "PSRF" = round(blavInspect(full.fit, "psrf"), 3) %>% .[1:16],
                  "Prior" = h5.out$uninf[1:16, 8],
                  "Mean (SD)" = paste0(h5.out$lowinf[1:16, 2], " (", h5.out$lowinf[1:16, 3], ")"),
                  "95% PPI" = paste(h5.out$lowinf[1:16, 4], h5.out$lowinf[1:16, 5], sep = ", "),
                  "neff" = h5.out$lowinf[1:16, 7],
                  "PSRF" = round(blavInspect(full_lowinf.fit, "psrf"), 3) %>% .[1:16],
                  "Prior" = h5.out$lowinf[1:16, 8],
                  "Mean (SD)" = paste0(h5.out$highinf[1:16, 2], " (", h5.out$highinf[1:16, 3], ")"),
                  "95% PPI" = paste(h5.out$highinf[1:16, 4], h5.out$highinf[1:16, 5], sep = ", "),
                  "neff" = h5.out$highinf[1:16, 7],
                  "PSRF" = round(blavInspect(full_highinf.fit, "psrf"), 3) %>% .[1:16],
                  "Prior" = h5.out$highinf[1:16, 8])

### Estimate bias
h5.estimates <- lapply(h5.out, function(x) as.numeric(x[, 2][1:16]))
h5.bias <- c(lowinf_uninf = round(100 * ((h5.estimates$lowinf - h5.estimates$uninf) / h5.estimates$uninf), 2),
             highinf_uninf = round(100 * ((h5.estimates$highinf - h5.estimates$uninf) / h5.estimates$uninf), 2))

