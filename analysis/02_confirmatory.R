# Load packages ---------------------------
list.of.packages <- c("ggplot2", "Rcpp", "tm", "SnowballC", "wordcloud", "RColorBrewer", "RCurl", "XML", "papaja", "tidyverse", "haven",
                      "emmeans", "TOSTER", "psy", "corrr", "MOTE", "knitr", "kableExtra", "codebook", "psych", "rlang", "blavaan",
                      "ggpubr", "bfrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source('http://www.sthda.com/upload/rquery_wordcloud.r')

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
corstars <- function(x, method = c("pearson", "spearman"), removeTriangle = c("upper", "lower"), result = c("none", "html", "latex")) {
  # Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix <- rcorr(x, type = method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value

  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****",
                    ifelse(p < .001, "*** ",
                           ifelse(p < .01, "**  ",
                                  ifelse(p < .05, "*   ", "    "))))

  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[, -1]

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep = "")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep = "")

  ## remove upper triangle of correlation matrix
  if (removeTriangle[1] == "upper") {
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }

  ## remove lower triangle of correlation matrix
  else if (removeTriangle[1] == "lower") {
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }

  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew) - 1])
  if (result[1] == "none") return(Rnew)
  else {
    if (result[1] == "html") print(xtable(Rnew), type = "html")
    else print(xtable(Rnew), type = "latex")
  }
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

# Load data ---------------------------
clean <- haven::read_spss("../data/02_clean_data.sav")

clean[, c("GENDER", "condition", "conditionbi", "genderbi", "agebi")] <- haven::as_factor(clean[, c("GENDER", "condition", "conditionbi", "genderbi", "agebi")])

# Data Overview ---------------------------
describeBy(clean, clean$condition)

# participants
data_desc <- c(total_n = nrow(raw),
               clean_n = nrow(clean),
               veg_n = nrow(raw) - nrow(noveg),
               comp_n = nrow(noveg) - nrow(clean))

# age
age_desc <- clean %>%
  summarise(min_age = min(AGE),
            max_age = max(AGE),
            m_age = printnum(mean(AGE)),
            sd_age = printnum(sd(AGE)))

# gender
gender_freq <- round(100 * prop.table(table(clean$GENDER)), digits = 2)

# word cloud
comb_text <- toString(clean[3:4]) %>%
  gsub('environmental', ' environment ', .) %>%
  gsub('animal | animals | animal welfare', ' animals ', .) %>%
  gsub('healthier', ' health ', .) %>%
  gsub('climate change | climate-change', ' climate ', .) %>%
  gsub('aware | awareness', ' awareness ', .) %>%
  gsub('red meat', ' red-meat ', .)

comb_cloud <- rquery.wordcloud(comb_text, type ="text", lang = "english", min.freq = 5,
                               excludeWords = c("reasons", "eating", "also", "meat", "people", "think", "much", "become", "due", "lot",
                                                "less", "eat", "consumption"))
head(comb_cloud$freqTable, 20)

# Randomization check ---------------------------

# age
agebycond <- clean %>%
  group_by(condition) %>%
  summarise(age_m = mean(AGE),
            age_sd = sd(AGE), .groups = "rowwise") %>%
  printnum()

age_check <- paste(agebycond$age_m, agebycond$age_sd, sep = ' $\\pm$ ')
age_stat <- apa_print(aov(AGE ~ condition, clean))

# political position
polbycond <- clean %>%
  group_by(condition) %>%
  summarise(pol_m = mean(POLITICS), pol_sd = sd(POLITICS), .groups = "rowwise") %>%
  printnum()

pol_check <- paste(polbycond$pol_m, polbycond$pol_sd, sep = ' $\\pm$ ')
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
random_check <- tibble(Item = c("Age (years)", "Gender (\\%)", "Political position", "Home country (\\%)"),
                       Dynamic = c(age_check[1], gender_check[1], pol_check[1], nation_check[1]),
                       Static = c(age_check[2], gender_check[2], pol_check[2], nation_check[2]),
                       None = c(age_check[3], gender_check[3], pol_check[3], nation_check[3]),
                       Test = c(age_stat$statistic$condition, gender_stat$statistic, pol_stat$statistic$condition, nation_stat$statistic)) %>%
  mutate_all(linebreak)

# Correlation matrix ---------------------------
# reliability
cron <- apply(matrix(6:14, ncol = 3), 2, function(x) printnum(cronbach(clean[x])$alpha))

# correlation matrix
measure_sum <- clean %>%
  summarise(m_interest = mean(INTEREST),
            m_attitude = mean(attitude_mean),
            m_expect = mean(expect_mean),
            m_intent = mean(intention_mean),
            m_change = mean(PERCEPTCHANGE),
            m_conformity = mean(PRECONFORMITY),
            sd_interest = sd(INTEREST),
            sd_attitude = sd(attitude_mean),
            sd_expect = sd(expect_mean),
            sd_intent = sd(intention_mean),
            sd_change = sd(PERCEPTCHANGE),
            sd_conformity = sd(PRECONFORMITY)) %>%
  printnum()

mcor <- clean %>% select("INTEREST", "attitude_mean", "expect_mean", "intention_mean", "PERCEPTCHANGE", "PRECONFORMITY") %>% corstars()

measure_sum <- tibble(Measure = c("1. Interest", "2. Attitudes", "3. Expectations", "4. Intentions", "5. Future norm", "6. Preconformity"),
                   Mean = unlist(measure_sum[, 1:6]),
                   SD = unlist(measure_sum[, 7:12]),
                   Alpha = c("-", cron, "-", "-"))

cortable <- cbind(measure_sum, mcor) %>% as_tibble()

# H1 interest ---------------------------
interest_desc <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(INTEREST),
            SD = sd(INTEREST), .groups = "rowwise")

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
  DYST = printnum(d.ind.t(interest_desc[[1, 3]], interest_desc[[2, 3]], interest_desc[[1, 4]], interest_desc[[2, 4]], interest_desc[[1, 2]], interest_desc[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(interest_desc[[1, 3]], interest_desc[[3, 3]], interest_desc[[1, 4]], interest_desc[[3, 4]], interest_desc[[1, 2]], interest_desc[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(interest_desc[[2, 3]], interest_desc[[3, 3]], interest_desc[[2, 4]], interest_desc[[3, 4]], interest_desc[[2, 2]], interest_desc[[3, 2]], a = .05)),
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

# Future consumption
future_cons <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTCHANGE, na.rm = TRUE),
            SD = sd(PERCEPTCHANGE), .groups = "rowwise")

ls_change <- lm(PERCEPTCHANGE ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_change$table <- ls_change$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2change.effect <- list(
  DYST = printnum(d.ind.t(future_cons[[1, 3]], future_cons[[2, 3]], future_cons[[1, 4]], future_cons[[2, 4]], future_cons[[1, 2]], future_cons[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(future_cons[[1, 3]], future_cons[[3, 3]], future_cons[[1, 4]], future_cons[[3, 4]], future_cons[[1, 2]], future_cons[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(future_cons[[2, 3]], future_cons[[3, 3]], future_cons[[2, 4]], future_cons[[3, 4]], future_cons[[2, 2]], future_cons[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_change$SE[1], obtained = ls_change$estimate[1], dfdata = ls_change$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_change$estimate[1], sample_se = ls_change$SE[1], sample_df = ls_change$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# Preconformity

# summary table
preconformity <- clean %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n = n(),
                   Mean = mean(PRECONFORMITY),
                   SD = sd(PRECONFORMITY), .groups = "rowwise")

# contrasts
ls_preconformity <- lm(PRECONFORMITY ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_preconformity$table <- ls_preconformity$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2preconformity.effect <- list(
  DYST = printnum(d.ind.t(preconformity[[1, 3]], preconformity[[2, 3]], preconformity[[1, 4]], preconformity[[2, 4]], preconformity[[1, 2]], preconformity[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(preconformity[[1, 3]], preconformity[[3, 3]], preconformity[[1, 4]], preconformity[[3, 4]], preconformity[[1, 2]], preconformity[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(preconformity[[2, 3]], preconformity[[3, 3]], preconformity[[2, 4]], preconformity[[3, 4]], preconformity[[2, 2]], preconformity[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_preconformity$SE[1], obtained = ls_preconformity$estimate[1], dfdata = ls_preconformity$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_preconformity$estimate[1], sample_se = ls_preconformity$SE[1], sample_df = ls_preconformity$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# combined

# summary table
comb_future <- clean %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n = n(),
                   Mean = mean(comb_future),
                   SD = sd(comb_future), .groups = "rowwise")

# contrasts
ls_comb_future <- lm(comb_future ~ condition, clean) %>%
  emmeans(., "condition", contr = contrasts) %>%
  .$contrasts %>%
  c(summary(.), apa_print(.))

ls_comb_future$table <- ls_comb_future$table %>%
  mutate(contrast = c("Dynamic, static", "Dynamic, control", "Static, control", "Dynamic, both", "Norms, control"),
         rownames = NULL)

H2comb_future.effect <- list(
  DYST = printnum(d.ind.t(comb_future[[1, 3]], comb_future[[2, 3]], comb_future[[1, 4]], comb_future[[2, 4]], comb_future[[1, 2]], comb_future[[2, 2]], a = .05)),
  DYNO = printnum(d.ind.t(comb_future[[1, 3]], comb_future[[3, 3]], comb_future[[1, 4]], comb_future[[3, 4]], comb_future[[1, 2]], comb_future[[3, 2]], a = .05)),
  STNO = printnum(d.ind.t(comb_future[[2, 3]], comb_future[[3, 3]], comb_future[[2, 4]], comb_future[[3, 4]], comb_future[[2, 2]], comb_future[[3, 2]], a = .05)),
  DYST.Bf = Bf(sd = ls_comb_future$SE[1], obtained = ls_comb_future$estimate[1], dfdata = ls_comb_future$df[1], meanoftheory = 0, sdtheory = 0.40, dftheory = 10^10, tail = 1),
  rr      = bfrr(sample_mean = ls_comb_future$estimate[1], sample_se = ls_comb_future$SE[1], sample_df = ls_comb_future$df[1], model = "normal", mean = 0, sd = 0.40, tail = 1, criterion = 5,
                 rr_interval = list(mean = c(-2, 2), sd = c(0, 2)), precision = 0.05))

# Secondary analyses ---------------------------

# Perception change by condition on numerican 1-100 scale
PERCEPTNUM_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTNUM),
            SD = sd(PERCEPTNUM), .groups = "rowwise")

TOSTtwo.sci(m1 = PERCEPTNUM_cond[[1, 3]], m2 = PERCEPTNUM_cond[[2, 3]], sd1 = PERCEPTNUM_cond[[1, 4]], sd2 = PERCEPTNUM_cond[[2, 4]], n1 = PERCEPTNUM_cond[[1, 2]], n2 = PERCEPTNUM_cond[[2, 2]], low_eqbound = -5, high_eqbound = 5)

# Perception change by condition on Likert scale
PERCEPTSCALE_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(PERCEPTSCALE),
            SD = sd(PERCEPTSCALE), .groups = "rowwise")

# construal of consumption
construal_cond <- clean %>%
  group_by(condition) %>%
  summarise(n = n(),
            Mean = mean(CONSTRUAL_1),
            SD = sd(CONSTRUAL_1), .groups = "rowwise") %>%
  mutate_if(is.numeric, round, 2)

