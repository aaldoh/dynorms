# Load packages ---------------------------
## devtools::install_github("crsh/papaja")
## devtools::install_github("debruine/bfrr")
req_packages <- c("MOTE", "MASS", "ggplot2", "Rcpp", "tm", "SnowballC", "wordcloud",
                  "RColorBrewer", "RCurl", "XML", "papaja", "tidyverse", "knitr",
                  "here", "kableExtra", "codebook", "psych", "rlang", "bfrr",
                  "qualtRics", "emmeans", "TOSTER", "psy")
new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) install.packages(new_packages)
lapply(req_packages, require, character.only = TRUE)

set.seed(2019)
options(mc.cores = parallel::detectCores()) ## Run chains in parallel
source('http://www.sthda.com/upload/rquery_wordcloud.r')

export <- qualtRics::read_survey(here("data/01_exported.csv")) %>%
  .[ , -c(1:11)]

# Create functions ---------------------------
sci.num <- function(x, zero.rm = F) {
  if (length(x) > 1) {
    sapply(x, sci.num, zero.rm)
  } else {
    if (abs(x) < .0005) {
      x <- sub("^(-?\\d+\\.?\\d{0,2}).*?e(-\\d+)$", "$\\1\\\\times10^{\\2}$", format(x, scientific = T))
      x <- sub("^\\$-?1\\\\times", "$", x)
    } else {
      x <- round(x, 3)
      if (x < 0) x <- sub("(.*)", "$\\1$", x)
      if (zero.rm) x <- sub("^(\\$-)?0", "\\1", x)
    }
    return(x)
  }
}
Bf <- function(sd, obtained, dfdata, meanoftheory, sdtheory, dftheory, tail = 2) {
  area <- 0
  normarea <- 0
  theta <- meanoftheory - 10 * sdtheory
  incr <- sdtheory / 200

  for (A in -2000:2000) {
    theta <- theta + incr
    dist_theta <- dt((theta - meanoftheory) / sdtheory, df = dftheory)
    if (identical(tail, 1)) {
      if (theta <= 0) {
        dist_theta <- 0

      } else {
        dist_theta <- dist_theta * 2
      }

    }
    height <- dist_theta * dt((obtained - theta) / sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta * incr

  }
  LikelihoodTheory <- area / normarea
  Likelihoodnull <- dt(obtained / sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  BayesFactor
}
glrstab<- function(x, mat.rows = c(...), export=FALSE) {

  r <-corr.test(x)$r	#taking just the correlation matrix; no N, or p
  p <-corr.test(x)$p	#taking the p*s

  #define notions for significance levels
  mystars <- ifelse(p < .001, "**"
                    , ifelse(p < .01, "**"
                             , ifelse(p < .05, "*"
                                      , ifelse(p < .10, "+", " "))))

  #round r, define new matrix Rnew with the correlations from rnd and paste mystars
  rnd  <- papaja::printnum(r, gt1 = FALSE, digits = 2)  #round, drop leading 0 - Thanks CRSH!
  Rnew <- matrix(paste(rnd, mystars, sep=""), ncol=ncol(rnd))

  #remove 1.0 correlations from diagonal  and set the strings
  diag(Rnew) <- ''
  Rnew[upper.tri(Rnew)] <- ''

  rownames(Rnew) <- paste(1:ncol(rnd), colnames(rnd), sep=" ")         #define number and name
  colnames(Rnew) <- paste(1:ncol(rnd), "", sep="") 			       #define number

  #fun-part: we trim the top half
  Rnew[upper.tri(Rnew)] <- ''
  Rnew

  Rnew <- cbind(round(data.frame(psych::describe(x))[,3:4],2), Rnew)		     #describe x, M sD - put them in the matrix
  colnames(Rnew)[1:2] <- c("M","SD")					      		#Beschriftung der neuen Spalten
  Rnew <- Rnew[,1:(ncol(Rnew)-1)]							        	#delete the last column (ugly)


  rownames(Rnew) <- mat.rows
  #export to clipboard

  if (export==TRUE){
    result<-write.table(Rnew
                        , "clipboard"
                        , sep=";"
                        , row.names=FALSE)
  }
  else result <- Rnew
  return(result)

}
TOSTtwo.sci <- function(m1, m2, sd1, sd2, n1, n2, low_eqbound, high_eqbound, alpha, var.equal, plot = TRUE, verbose = TRUE) {

  if (missing(alpha)) {
    alpha <- 0.05
  }
  if (missing(var.equal)) {
    var.equal <- FALSE
  }
  if (low_eqbound >= high_eqbound)
    warning("The lower bound is equal to or larger than the upper bound. Check the plot and output to see if the bounds are specified as you intended.")
  if (n1 < 2 | n2 < 2)
    stop("The sample size should be larger than 1.")
  if (1 <= alpha | alpha <= 0)
    stop("The alpha level should be a positive value between 0 and 1.")
  if (sd1 <= 0 | sd2 <= 0)
    stop("The standard deviation should be a positive value.")
  if (var.equal == TRUE) {
    sdpooled <- sqrt((((n1 - 1) * (sd1^2)) + (n2 - 1) *
                        (sd2^2)) / ((n1 + n2) - 2))
    t1 <- ((m1 - m2) - low_eqbound) / (sdpooled * sqrt(1 / n1 +
                                                         1 / n2))
    degree_f <- n1 + n2 - 2
    p1 <- pt(t1, degree_f, lower.tail = FALSE)
    t2 <- ((m1 - m2) - high_eqbound) / (sdpooled * sqrt(1 / n1 +
                                                          1 / n2))
    p2 <- pt(t2, degree_f, lower.tail = TRUE)
    LL90 <- (m1 - m2) - qt(1 - alpha, degree_f) * (sdpooled *
                                                     sqrt(1 / n1 + 1 / n2))
    UL90 <- (m1 - m2) + qt(1 - alpha, degree_f) * (sdpooled *
                                                     sqrt(1 / n1 + 1 / n2))
    LL95 <- (m1 - m2) - qt(1 - (alpha / 2), degree_f) * (sdpooled *
                                                           sqrt(1 / n1 + 1 / n2))
    UL95 <- (m1 - m2) + qt(1 - (alpha / 2), degree_f) * (sdpooled *
                                                           sqrt(1 / n1 + 1 / n2))
    t <- (m1 - m2) / (sdpooled * sqrt(1 / n1 + 1 / n2))
    pttest <- 2 * pt(-abs(t), df = degree_f)
  }
  else {
    sdpooled <- sqrt((sd1^2 + sd2^2) / 2)
    t1 <- ((m1 - m2) - low_eqbound) / sqrt(sd1^2 / n1 + sd2^2 / n2)
    degree_f <- (sd1^2 / n1 + sd2^2 / n2)^2 / (((sd1^2 / n1)^2 / (n1 -
                                                                    1)) + ((sd2^2 / n2)^2 / (n2 - 1)))
    p1 <- pt(t1, degree_f, lower.tail = FALSE)
    t2 <- ((m1 - m2) - high_eqbound) / sqrt(sd1^2 / n1 + sd2^2 / n2)
    p2 <- pt(t2, degree_f, lower.tail = TRUE)
    t <- (m1 - m2) / sqrt(sd1^2 / n1 + sd2^2 / n2)
    pttest <- 2 * pt(-abs(t), df = degree_f)
    LL90 <- (m1 - m2) - qt(1 - alpha, degree_f) * sqrt(sd1^2 / n1 +
                                                         sd2^2 / n2)
    UL90 <- (m1 - m2) + qt(1 - alpha, degree_f) * sqrt(sd1^2 / n1 +
                                                         sd2^2 / n2)
    LL95 <- (m1 - m2) - qt(1 - (alpha / 2), degree_f) * sqrt(sd1^2 / n1 +
                                                               sd2^2 / n2)
    UL95 <- (m1 - m2) + qt(1 - (alpha / 2), degree_f) * sqrt(sd1^2 / n1 +
                                                               sd2^2 / n2)
  }
  ptost <- max(p1, p2)
  ttost <- ifelse(abs(t1) < abs(t2), t1, t2)
  dif <- (m1 - m2)
  testoutcome <- ifelse(pttest < alpha, "significant", "non-significant")
  TOSToutcome <- ifelse(ptost < alpha, "significant", "non-significant")
  if (plot == TRUE) {
    plot(NA, ylim = c(0, 1), xlim = c(min(LL90, low_eqbound) -
                                        max(UL90 - LL90, high_eqbound - low_eqbound) / 10,
                                      max(UL90, high_eqbound) + max(UL90 - LL90, high_eqbound -
                                                                      low_eqbound) / 10), bty = "l", yaxt = "n", ylab = "",
         xlab = "Mean Difference")
    points(x = dif, y = 0.5, pch = 15, cex = 2)
    abline(v = high_eqbound, lty = 2)
    abline(v = low_eqbound, lty = 2)
    abline(v = 0, lty = 2, col = "grey")
    segments(LL90, 0.5, UL90, 0.5, lwd = 3)
    segments(LL95, 0.5, UL95, 0.5, lwd = 1)
    title(main = paste("Equivalence bounds ", round(low_eqbound,
                                                    digits = 3), " and ", round(high_eqbound, digits = 3),
                       "\nMean difference = ", round(dif, digits = 3),
                       " \n TOST: ", 100 * (1 - alpha * 2), "% CI [", round(LL90,
                                                                            digits = 3), ";", round(UL90, digits = 3), "] ",
                       TOSToutcome, " \n NHST: ", 100 * (1 - alpha), "% CI [",
                       round(LL95, digits = 3), ";", round(UL95, digits = 3),
                       "] ", testoutcome, sep = ""), cex.main = 1)
  }
  if (missing(verbose)) {
    verbose <- TRUE
  }
  if (verbose == TRUE) {
    cat("TOST results:\n")
    cat("*t*-value lower bound:", sci.num(t1), "\t*p*-value lower bound:", sci.num(p1, T))
    cat("\n")
    cat("*t*-value upper bound:", sci.num(t2), "\t*p*-value upper bound:", sci.num(p2, T))
    cat("\n")
    cat("degrees of freedom :", round(degree_f, digits = 2))
    cat("\n\n")
    cat("Equivalence bounds (raw scores):")
    cat("\n")
    cat("low eqbound:", sci.num(low_eqbound),
        "\nhigh eqbound:", sci.num(high_eqbound))
    cat("\n\n")
    cat("TOST confidence interval:")
    cat("\n")
    cat("lower bound ", 100 * (1 - alpha * 2), "% CI: ",
        sci.num(LL90), "\nupper bound ",
        100 * (1 - alpha * 2), "% CI:  ", sci.num(UL90), sep = "")
    cat("\n\n")
    cat("NHST confidence interval:")
    cat("\n")
    cat("lower bound ", 100 * (1 - alpha), "% CI: ", sci.num(LL95), "\nupper bound ", 100 * (1 - alpha),
        "% CI:  ", sci.num(UL95), sep = "")
    cat("\n\n")
    cat("Equivalence Test Result:\n")
    message(cat("The equivalence test was ", TOSToutcome,
                ", *t*(", round(degree_f, digits = 2), ") = ", sci.num(ttost),
                ", *p* = ", sci.num(ptost, T), ", given equivalence bounds of ",
                sci.num(low_eqbound),
                " and ", sci.num(high_eqbound), " (on a raw scale) and an alpha of ",
                alpha, ".", sep = ""))
    cat("\n")
    cat("Null Hypothesis Test Result:\n")
    message(cat("The null hypothesis test was ", testoutcome,
                ", *t*(", round(degree_f, digits = 2), ") = ", sci.num(t),
                ", *p* = ", sci.num(pttest, T), ", given an alpha of ",
                alpha, ".", sep = ""))
    if (pttest <= alpha && ptost <= alpha) {
      combined_outcome <- "statistically different from zero and statistically equivalent to zero"
    }
    if (pttest < alpha && ptost > alpha) {
      combined_outcome <- "statistically different from zero and statistically not equivalent to zero"
    }
    if (pttest > alpha && ptost <= alpha) {
      combined_outcome <- "statistically not different from zero and statistically equivalent to zero"
    }
    if (pttest > alpha && ptost > alpha) {
      combined_outcome <- "statistically not different from zero and statistically not equivalent to zero"
    }
    cat("\n")
    message(cat("Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is ",
                combined_outcome, ".", sep = ""))
  }
  invisible(list(diff = dif, TOST_t1 = t1, TOST_p1 = p1, TOST_t2 = t2,
                 TOST_p2 = p2, TOST_df = degree_f, alpha = alpha, low_eqbound = low_eqbound,
                 high_eqbound = high_eqbound, LL_CI_TOST = LL90, UL_CI_TOST = UL90,
                 LL_CI_TTEST = LL95, UL_CI_TTEST = UL95))
}
# Data cleaning ---------------------------
raw <- export %>%
  mutate(across(c(6:23), as.numeric),
         across(c(where(is.character), UKNATION, RESIDENT,VEG, GENDER, DQInclude, -DYUK, -STUK), as_factor),
         condition = ifelse(!is.na(DYUK), 1, ifelse(!is.na(STUK), 2, ifelse(is.na(c(DYUK, STUK)), 3, NA))) %>%
           factor(., labels = c("Dynamic", "Static", "No norm")),
         attitude_mean = rowMeans(across(ATT1:ATT3), na.rm = T),
         intention_mean = rowMeans(across(INTENT1:INTENT3), na.rm = T),
         expect_mean = rowMeans(across(EXPECT1:EXPECT3), na.rm = T),
         expintent_avg = rowMeans(across(intention_mean:expect_mean), na.rm = T),
         conditionbi = na_if(condition, "No norm"),
         genderbi = na_if(GENDER, 3),
         agebi = cut(export$AGE, c(0, median(export$AGE, na.rm = T), max(export$AGE, na.rm = T)), labels = c("younger", "older")),
         age_cent = AGE - mean(AGE, na.rm = TRUE),
         PERCEPTCHANGE_r = 8 - PERCEPTCHANGE,
         comb_future = rowMeans(across(c("PERCEPTCHANGE_r", "PRECONFORMITY")), na.rm = T)) %>%
  droplevels() %>%
  filter(VEG == 2, !is.na(DQInclude)) %>%
  dplyr::select(where(is.character), where(is.factor), where(is.numeric), -VEG)

raw <- cbind(raw, dummy.code(raw$condition))

# excluding data quality fails
clean <- raw %>%
  filter(INTENT_QUALITY == 6 & DQInclude == 1) %>%
  dplyr::select(-INTENT_QUALITY, -DQInclude)

# outliers
mcd     <- cov.mcd(clean[ , c(11, 28, 31, 21:25)], quantile.used = nrow(clean[ , c(11, 28, 31, 21:25)])*.75)
mcd_mah <- mahalanobis(clean[ , c(11, 28, 31, 21:25)], mcd$center,mcd$cov)
cutoff  <- qchisq(p = 0.99, df = ncol(clean[ , c(11, 28, 31, 21:25)]))
no_out  <- clean[mcd_mah <= cutoff, ]

# word cloud
comb_text <- toString(clean[1:2]) %>%
  gsub('environmental', ' environment ', .) %>%
  gsub('animal | animals | animal welfare', ' animals ', .) %>%
  gsub('healthier', ' health ', .) %>%
  gsub('climate change | climate-change', ' climate ', .) %>%
  gsub('aware | awareness', ' awareness ', .) %>%
  gsub('red meat', ' red-meat ', .)

# Data Overview ---------------------------
age_desc <- apply(clean[27], MARGIN = 2, function(x) c(min = min(x), max = max(x), avg = mean(x), sd = sd(x))) ## age
gender_freq <- 100 * prop.table(table(clean$gender)) ## gender
cron <- apply(matrix(12:20, ncol = 3), 2, function(x) printnum(cronbach(clean[x])$alpha))
measured_vars <- clean[, c(11, 28:31, 33, 25)]
measure.tib <- glrstab(measured_vars, mat.rows = c("Interest", "Attitude", "Intention", "Expectation", "Intent/expectation composite", "Perception of change", "Preconformity")) %>%
  cbind(Alpha = c("-", cron, "-", "-", "-"), .)

# Randomization check ---------------------------
demobycond <- group_by(clean, condition) %>%
  summarise(across(c("AGE", "POLITICS"), list(mean = mean, sd = sd), .names = "{.fn}_{.col}"), .groups = "rowwise") %>%
  printnum()
age_check <- paste(demobycond$mean_AGE, demobycond$sd_AGE, sep = ' $\\pm$ ')
age_stat <- apa_print(aov(AGE ~ condition, clean))
pol_check <- paste(demobycond$mean_POLITICS, demobycond$sd_POLITICS, sep = ' $\\pm$ ')
pol_stat <- apa_print(aov(POLITICS ~ condition, clean))

# gender
gender_check <- clean %>%
  count(condition, GENDER) %>%
  group_by(condition) %>%
  mutate(rel.freq = paste0(round(100 * n / sum(n), 2), "\\%")) %>%
  summarise(Gender = paste0(GENDER, " (", rel.freq, ")", collapse = "\n"), .groups = "rowwise") %>%
  pull(Gender)

gender_stat <- apa_print(chisq.test(clean$condition, clean$GENDER), n = nrow(clean))

# nation
nation_check <- clean %>%
  count(condition, UKNATION) %>%
  group_by(condition) %>%
  mutate(rel.freq = round(100 * n / sum(n), 2)) %>%
  summarise(Nationality = paste0(UKNATION, " (", rel.freq, "\\%)", collapse = "\n"), .groups = "rowwise") %>%
  pull(Nationality)

nation_stat <- apa_print(chisq.test(clean$condition, clean$UKNATION), n = nrow(clean))

# table
random_check <- rbind(age_check, gender_check, pol_check, nation_check) %>%
  cbind(c("Age (years)", "Gender (\\%)", "Political position", "Nationality (\\%)"), .,
        c(age_stat$statistic$condition, gender_stat$statistic, pol_stat$statistic$condition, nation_stat$statistic))

# H1 interest ---------------------------
outcomes_desc <- group_by(clean, condition) %>%
  summarise(n = n(),
            across(c("INTEREST", "attitude_mean", "intention_mean", "expect_mean"), list(mean = mean, sd = sd),
                   .names = "{.fn}_{.col}"), .groups = "rowwise")
# contrasts
contrasts = list(DYST = c(1, -1, 0),
                 DYNO = c(1, 0, -1),
                 STNO = c(0, 1, -1),
                 DYCONT = c(1, -0.5, -0.5),
                 EXPNO = c(-0.5, -0.5, 1))

ls_interest <- lm(INTEREST ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_interest$table <- ls_interest$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

# effect sizes
H1.effect <- list(
  DYST = printnum(d.ind.t(outcomes_desc[[1, 3]], outcomes_desc[[2, 3]], outcomes_desc[[1, 4]], outcomes_desc[[2, 4]], outcomes_desc[[1, 2]], outcomes_desc[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(outcomes_desc[[1, 3]], outcomes_desc[[3, 3]], outcomes_desc[[1, 4]], outcomes_desc[[3, 4]], outcomes_desc[[1, 2]], outcomes_desc[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(outcomes_desc[[2, 3]], outcomes_desc[[3, 3]], outcomes_desc[[2, 4]], outcomes_desc[[3, 4]], outcomes_desc[[2, 2]], outcomes_desc[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_interest$SE[1], obtained = ls_interest$estimate[1], dfdata = ls_interest$df[1], meanoftheory = 0,
                 sdtheory = 0.69, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_interest$estimate[1], sample_se = ls_interest$SE[1], sample_df = ls_interest$df[1], model = "normal", mean = 0, sd = 0.69, tail = 1, criterion = 5,
                rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# full regression
pol.interest_out <- apa_print(summary(lm(INTEREST ~ POLITICS, clean)))
gen.interest_out <- apa_print(t.test(INTEREST ~ genderbi, clean))
demoreg.out <- summary(lm(INTEREST ~ conditionbi + genderbi + POLITICS + conditionbi:genderbi + conditionbi:POLITICS, clean)) %>%
  apa_print()

demoreg.out$table <- demoreg.out$table %>%
  mutate(predictor = c("Intercept", "Condition", "Gender", "Political position", "Condition $\\times$ Gender", "Condition $\\times$ Political position"),
         rownames = NULL)

# H2 ---------------------------
future <- group_by(clean, condition) %>%
  summarise(n = n(),
            across(c("PERCEPTCHANGE_r", "PRECONFORMITY", "comb_future"), list(mean = mean, sd = sd),
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

# Future consumption
ls_change <- lm(PERCEPTCHANGE_r ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_change$table <- ls_change$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2change.effect <- list(
  DYST = printnum(d.ind.t(future[[1, 3]], future[[2, 3]], future[[1, 4]], future[[2, 4]], future[[1, 2]], future[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(future[[1, 3]], future[[3, 3]], future[[1, 4]], future[[3, 4]], future[[1, 2]], future[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(future[[2, 3]], future[[3, 3]], future[[2, 4]], future[[3, 4]], future[[2, 2]], future[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_change$SE[1], obtained = ls_change$estimate[1], dfdata = ls_change$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_change$estimate[1], sample_se = ls_change$SE[1], sample_df = ls_change$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# Preconformity
# contrasts
ls_preconformity <- lm(PRECONFORMITY ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_preconformity$table <- ls_preconformity$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2preconformity.effect <- list(
  DYST = printnum(d.ind.t(future[[1, 5]], future[[2, 5]], future[[1, 6]], future[[2, 6]], future[[1, 2]], future[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(future[[1, 5]], future[[3, 5]], future[[1, 6]], future[[3, 6]], future[[1, 2]], future[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(future[[2, 5]], future[[3, 5]], future[[2, 6]], future[[3, 6]], future[[2, 2]], future[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_preconformity$SE[1], obtained = ls_preconformity$estimate[1], dfdata = ls_preconformity$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_preconformity$estimate[1], sample_se = ls_preconformity$SE[1], sample_df = ls_preconformity$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# combined
# contrasts
ls_comb_future <- lm(comb_future ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_comb_future$table <- ls_comb_future$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2comb_future.effect <- list(
  DYST = printnum(d.ind.t(future[[1, 7]], future[[2, 7]], future[[1, 8]], future[[2, 8]], future[[1, 2]], future[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(future[[1, 7]], future[[3, 7]], future[[1, 8]], future[[3, 8]], future[[1, 2]], future[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(future[[2, 7]], future[[3, 7]], future[[2, 8]], future[[3, 8]], future[[2, 2]], future[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_comb_future$SE[1], obtained = ls_comb_future$estimate[1], dfdata = ls_comb_future$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_comb_future$estimate[1], sample_se = ls_comb_future$SE[1], sample_df = ls_comb_future$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# Secondary analyses ---------------------------
secondary <- group_by(clean, condition) %>%
  summarise(n = n(),
            across(c("PERCEPTNUM", "PERCEPTSCALE", "CONSTRUAL_1"),
                   list(mean = mean, sd = sd),
                   .names = "{.fn}_{.col}"), .groups = "rowwise")

clean %>%
  filter(!is.na(conditionbi)) %>%
  saveRDS(., file = "../data/blavaan_data.rds")
