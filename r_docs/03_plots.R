# Setup ---------------------------
library(blavaan)
library(ggplot2)
library(bayesplot)
library(papaja)
library(tidyverse)

bintervals.plot <- function(data, discrete_labs = NULL, discrete_lims = NULL, color_group = NULL) {
  ggplot(data = data, aes(x = m, y = parameter, color = color_group)) +
    geom_linerange(aes(xmin = l, xmax = h), position = position_dodge(width=0.8), size=2)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = position_dodge(width=0.8)) +
    geom_point(size = 3, position = position_dodge(width=0.8)) +
    geom_text(aes(label = paste0("M = ", round(m, 2), " [", round(ll, 2), ", ", round(hh, 2), "]")),
              position = position_dodge(width=0.8), vjust = -1, size = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    labs(x = "Parameter Value", y = "Parameter", color = "Model",
         title = "Posterior distributions", subtitle = "with means and credibility intervals") +
    scale_y_discrete(labels = discrete_labs, limits = discrete_lims) +
    xlim(-0.5,0.5) +
    papaja::theme_apa(box = TRUE) +
    theme(legend.position = "right")
}
bdensity.plot <- function(data, discrete_labs) {
  ggplot(data = data, aes_(x = ~x, y = ~parameter, color = Model)) +
    geom_linerange(aes(xmin = l, xmax = h), position = position_dodge(width=0.8), size=2)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = position_dodge(width=0.8)) +
    geom_point(position = position_dodge(width=0.8)) +
    geom_text(aes(label = paste0("M = ", round(m, 2), " [", round(ll, 2), ", ", round(hh, 2), "]")),
              position = position_dodge(width=0.8), vjust = -1, size = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    labs(x = "Parameter Value", y = "Parameter (~ Condition)", title = "Posterior distributions", subtitle = "with means and 95% credibility intervals") +
    scale_y_discrete(labels = discrete_labs) +
    papaja::theme_apa(box = TRUE) +
    theme(legend.position = "right")
}

load("data/exp_results.RData")
set.seed(2019)

# H3 ---------------------------
## extracting MCMC chains and intervals
h3.mcmc <- lapply(h3.fit, function(x) as.matrix(blavInspect(x, 'mcmc')))

h3.mcmc_intervals <- bind_rows(lapply(h3.mcmc, function(x)
  mcmc_intervals_data(x, regex_pars = "^bet", point_est = "mean", prob = 0.8, prob_outer = 0.95))) %>%
  mutate(Model = rep(c("1: Uninformative", "2: Weakly informative", "3: Moderately informative"), each = 3))

h3.density <- bind_rows(lapply(h3.mcmc, function(x)
  mcmc_areas_data(x, regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean"))) %>%
  mutate(Model = rep(c("1: Uninformative", "2: Weakly informative", "3: Moderately informative"), each = 6211) %>%
           as.factor())

## traceplots
h3.traceplots <- lapply(h3.fit, function(x) plot(x, type = "trace"))

### intervals
h3_intervals.plot <- bintervals.plot(h3.mcmc_intervals,
                                     color_group = h3.mcmc_intervals$Model,
                                     discrete_labs = c("bet_sign[3]" = "Interest ~ Condition", "bet_sign[2]" = "Attitude ~ Condition", "bet_sign[1]" = "Intent/Exp ~ Condition"))

## density plots
h3_density.plot <- lapply(h3.mcmc, function(x) mcmc_areas(x, regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean") +
                       scale_y_discrete(labels = c("bet_sign[3]" = "Interest", "bet_sign[2]" = "Attitude", "bet_sign[1]" = "Intent/Exp")) +
                       xlim(-0.5,0.5) +
                       theme_apa(box = TRUE))

test<- bayesplot_grid(plots = h3_density.plot, subtitles = c("1. Uninformative", "2. Weakly informative", "3. Moderately informative"),
                      grid_args = list(ncol = 1))

# H4 ---------------------------
## extracting MCMC chains and intervals
h4.mcmc <- lapply(h4.fit, function(x) as.matrix(blavInspect(x, 'mcmc')))

h4.mcmc_intervals <- bind_rows(lapply(h4.mcmc, function(x)
  mcmc_intervals_data(x, regex_pars = "^bet", point_est = "mean", prob = 0.8, prob_outer = 0.95))) %>%
  mutate(Model = rep(c("1: Uninformative", "2: Weakly informative", "3: Moderately informative"), each = 12))

h4.density <- lapply(h4.mcmc, function(x)
  mcmc_areas_data(x, regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean"))
h4.density[[1]]$Model <- rep(c("1: Uninformative"))
h4.density[[2]]$Model <- rep(c("2: Weakly informative"))
h4.density[[3]]$Model <- rep(c("2: Moderately informative"))
h4.density <- bind_rows(h4.density) %>% mutate(Model = as_factor(Model))

## plots
h4.traceplots <- lapply(h4.fit, function(x) plot(x, type = "trace"))

h4_intervals.plot <- bintervals.plot(h4.mcmc_intervals[h4.mcmc_intervals$Model == "1: Uninformative",],
                                     discrete_labs = c("bet_sign[1]" = "Intent/Exp ~ Condition", "bet_sign[2]" = "Attitude ~ Condition", "bet_sign[3]" = "Interest ~ Condition",
                                                       "bet_sign[4]" = "Interest ~ Age", "bet_sign[5]" = "Attitude ~ Age", "bet_sign[6]" = "Intent/Exp ~ Age",
                                                       "bet_sign[7]" = "Interest ~ Gender", "bet_sign[8]" = "Attitude ~ Gender", "bet_sign[9]" = "Intent/Exp ~ Gender",
                                                       "bet_sign[10]" = "Interest ~ Politics", "bet_sign[11]" = "Attitude ~ Politics", "bet_sign[12]" = "Intent/Exp ~ Politics"),
                                     discrete_lims = c("bet_sign[1]", "bet_sign[2]", "bet_sign[3]", "bet_sign[4]",
                                                       "bet_sign[5]", "bet_sign[6]", "bet_sign[7]", "bet_sign[8]",
                                                       "bet_sign[9]", "bet_sign[10]", "bet_sign[11]", "bet_sign[12]"))

# H5 ---------------------------
## Trace plots
h5.traceplots <- lapply(h5.fit, function(x) plot(x, type = "trace"))
