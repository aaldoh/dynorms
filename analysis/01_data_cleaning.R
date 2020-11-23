# Title: Data cleaning

# Load packages ---------------------------
list.of.packages <- c("papaja", "tidyverse", "haven", "rlang", "psych")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# Load data ---------------------------
export <- haven::read_spss("../data/01_exported.sav")
complete <- filter(export, !is.na(DQInclude))
raw <- select(complete, -1:-13, -39)

# Data cleaning ---------------------------
raw[, c(1, 2, 23:25)] <- haven::as_factor(raw[, c(1, 2, 23:25)])

raw <- raw %>%
  mutate(
    condition = ifelse(DYUK != "", 1, ifelse(STUK != "", 2, ifelse(DYUK == "" & STUK == "", 3, NA))) %>%
      factor(., labels = c("Dynamic", "Static", "No norm")),
    attitude_mean = select(., ATT1, ATT2, ATT3) %>% rowMeans(., na.rm = TRUE),
    intention_mean = select(., INTENT1, INTENT2, INTENT3) %>% rowMeans(., na.rm = TRUE),
    expect_mean = select(., EXPECT1, EXPECT2, EXPECT3) %>% rowMeans(., na.rm = TRUE),
    conditionbi = na_if(condition, "No norm"),
    genderbi = na_if(GENDER, "Other"),
    agebi = cut(raw$AGE, c(0, median(raw$AGE, na.rm = T), max(raw$AGE, na.rm = T)), labels = c("younger", "older")),
    age_cent = AGE - mean(AGE, na.rm = TRUE),
    PERCEPTCHANGE = 8 - PERCEPTCHANGE %>% #recoded  - higher scores indicate increasing limiting meat cons
      set_attrs(., labels = c("Increase significantly" = 1, "Increase moderately" = 2, "Increase slightly" = 3, "Stay the same" = 4,
                              "Decrease slightly" = 5, "Decrease moderately" = 6, "Decrease significantly" = 7)),
    comb_future = select(., PERCEPTCHANGE, PRECONFORMITY) %>% rowMeans(., na.rm = T)) %>%
  droplevels()

raw <- cbind(raw, dummy.code(raw$conditionbi))

# excluding vegetarians/vegan
noveg <- filter(raw, VEG == "No")

# excluding data quality fails
clean <- noveg %>%
  filter(INTENT_QUALITY == 6 | DQInclude != 1) %>%
  select(-INTENT_QUALITY, -VEG, -DQInclude)

# outliers
mahalfiltered = mahalanobis(clean[ , 5:21],
                            colMeans(clean[ , 5:21], na.rm = T),
                            cov(clean[ , 5:21]))
cutoff = qchisq(1-.001, ncol(clean[ , 5:21]))
ncol(clean[ , 5:21]) #df
noout = subset(clean, mahalfiltered < cutoff)

# write datasets
write_sav(clean, "../data/02_clean_data.sav")
write_sav(noout, "../data/03_no_outliers.sav")

data_desc <- c(total_n = nrow(raw),
               clean_n = nrow(clean),
               veg_n = nrow(raw) - nrow(noveg),
               comp_n = nrow(noveg) - nrow(clean))


