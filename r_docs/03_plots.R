# Setup ---------------------------
set.seed(2019)
library(blavaan)
library(ggplot2)
library(bayesplot)
library(papaja)
library(tidyverse)

bintervals.plot <- function(data, discrete_labs = NULL, discrete_lims = NULL, color_group = NULL, xlims = limits(c(...), "x")) {
  ggplot(data = data, aes(x = m, y = parameter, color = color_group)) +
    geom_linerange(aes(xmin = l, xmax = h), position = position_dodge(width=0.8), size=1.5)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = position_dodge(width=0.8)) +
    geom_point(size = 2.5, position = position_dodge(width=0.8)) +
    geom_text(aes(label = paste0("M = ", round(m, 2), " [", round(ll, 2), ", ", round(hh, 2), "]")),
              position = position_dodge(width=0.8), vjust = -1, size = 3) +
    geom_vline(xintercept = 0, linetype = 3, alpha = 0.4) +
    labs(x = "Parameter Value", y = "Regression Parameter", color = "Model",
         title = "Posterior distributions", subtitle = "with means and credibility intervals") +
    scale_y_discrete(labels = discrete_labs, limits = discrete_lims) +
    xlim(xlims) +
    papaja::theme_apa(box = TRUE) +
    theme(legend.position = "right", legend.key.height = unit(1.3, "cm"))
}

load("data/blav_fit/h3_fit.RData")
load("data/blav_fit/h4_fit.RData")
load("data/blav_fit/h5_fit.RData")

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

## intervals
h3_intervals.plot <- bintervals.plot(h3.mcmc_intervals, xlims = c(-0.5,0.5),
                                     color_group = h3.mcmc_intervals$Model,
                                     discrete_labs = c("bet_sign[3]" = "Interest ~ Condition", "bet_sign[2]" = "Attitude ~ Condition", "bet_sign[1]" = "Intent/Exp ~ Condition")) +
  scale_color_manual(name="Beta prior",
                     limits = c("3: Moderately informative", "2: Weakly informative", "1: Uninformative"),
                     labels=c("1. Uninformative\nN(0,10)", "2. Weakly informative\nN(0.5,0.75)",
                              "3. Moderately informative\nN(0.5,0.35)"),
                     values=c("#33BBEE", "#EE7733", "#CC3311"))
## density plots
h3_density.plot <- lapply(h3.mcmc, function(x) mcmc_areas(x, regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean") +
                            scale_y_discrete(labels = c("bet_sign[3]" = "Interest", "bet_sign[2]" = "Attitude", "bet_sign[1]" = "Intent/Exp")) +
                            xlim(-0.5,0.5) +
                            theme_apa(box = TRUE))

h3_density.grid <- bayesplot_grid(plots = h3_density.plot, subtitles = c("1: Uninformative prior, N(0, 10)", "2: Weakly informative prior, N(0.5, 0.75)", "3: Moderately informative prior, N(0.5, 0.35)"),
                      grid_args = list(ncol = 1, left = "Parameter (~ Condition)", bottom = "Parameter Value", top="Posterior distributions with means 80% credibility intervals"))

# H4 ---------------------------
## extracting MCMC chains and intervals
h4.mcmc <- lapply(h4.fit, function(x) as.matrix(blavInspect(x, 'mcmc')))

h4.mcmc_intervals <- bind_rows(lapply(h4.mcmc, function(x)
  mcmc_intervals_data(x, regex_pars = "^bet", point_est = "mean", prob = 0.8, prob_outer = 0.95))) %>%
  mutate(Model = rep(c("1: Uninformative", "2: Weakly informative", "3: Moderately informative"), each = 12))

h4.density <- lapply(h4.mcmc, function(x)
  mcmc_areas_data(x, regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean"))

## traceplots
h4.traceplots <- lapply(h4.fit, function(x) plot(x, type = "trace"))

## intervals
h4_intervals.plot <- bintervals.plot(h4.mcmc_intervals[h4.mcmc_intervals$Model == "1: Uninformative",],
                                     xlims = c(-1.1,1.1),
                                     discrete_labs = c("bet_sign[1]" = "Intent/Exp ~ Condition", "bet_sign[2]" = "Attitude ~ Condition", "bet_sign[3]" = "Interest ~ Condition",
                                                       "bet_sign[4]" = "Interest ~ Age", "bet_sign[5]" = "Attitude ~ Age", "bet_sign[6]" = "Intent/Exp ~ Age",
                                                       "bet_sign[7]" = "Interest ~ Gender", "bet_sign[8]" = "Attitude ~ Gender", "bet_sign[9]" = "Intent/Exp ~ Gender",
                                                       "bet_sign[10]" = "Interest ~ Politics", "bet_sign[11]" = "Attitude ~ Politics", "bet_sign[12]" = "Intent/Exp ~ Politics"),
                                     discrete_lims = c("bet_sign[1]", "bet_sign[2]", "bet_sign[3]", "bet_sign[4]",
                                                       "bet_sign[5]", "bet_sign[6]", "bet_sign[7]", "bet_sign[8]",
                                                       "bet_sign[9]", "bet_sign[10]", "bet_sign[11]", "bet_sign[12]"))

## density plots
h4_density.plot <- mcmc_areas(h4.mcmc[[1]], regex_pars = "^bet", prob = 0.8, prob_outer = 0.99, point_est = "mean") +
  scale_y_discrete(labels = c("bet_sign[1]" = "Intent/Exp ~ Condition", "bet_sign[2]" = "Attitude ~ Condition", "bet_sign[3]" = "Interest ~ Condition",
                              "bet_sign[4]" = "Interest ~ Age", "bet_sign[5]" = "Attitude ~ Age", "bet_sign[6]" = "Intent/Exp ~ Age",
                              "bet_sign[7]" = "Interest ~ Gender", "bet_sign[8]" = "Attitude ~ Gender", "bet_sign[9]" = "Intent/Exp ~ Gender",
                              "bet_sign[10]" = "Interest ~ Politics", "bet_sign[11]" = "Attitude ~ Politics", "bet_sign[12]" = "Intent/Exp ~ Politics"),
                   limits = c("bet_sign[1]", "bet_sign[2]", "bet_sign[3]", "bet_sign[4]",
                              "bet_sign[5]", "bet_sign[6]", "bet_sign[7]", "bet_sign[8]",
                              "bet_sign[9]", "bet_sign[10]", "bet_sign[11]", "bet_sign[12]")) +
  labs(x = "Parameter Value", y = "Regression Parameter", color = "Model",
       title = "Posterior distributions", subtitle = "with means and 80% credibility intervals") +
  xlim(-0.5,0.5) +
  theme_apa(box = TRUE)

# H5 ---------------------------
## Trace plots
h5.traceplots <- lapply(h5.fit, function(x) plot(x, type = "trace"))

save(h3.traceplots, h4.traceplots, h5.traceplots,
     h3.mcmc_intervals, h4.mcmc_intervals,
     h3_intervals.plot, h4_intervals.plot,
     h3_density.grid, h4_density.plot,
     file =  "figures/exp_plots.RData")

