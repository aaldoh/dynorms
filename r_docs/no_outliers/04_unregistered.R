set.seed(2019)
options(mc.cores = parallel::detectCores()) ## Run chains in parallel

lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T, unload = T,force = T))

library(tidyverse)
library(papaja)
library(blavaan)
library(pwr)

blav_data_out <- readRDS("data/blav_data_no_out.rds")
no_out <- readRDS("data/no_out.rds")

# Power analysis
h1.pwr <- pwr.t.test(d = 0.36, sig.level = 0.05, power = .80, type = "two.sample", alternative= "two.sided")

# Mediation ---------------------------
##Is attitude a mediator of the effect of trending minority norms vs. minority only on intentions/expectations?

med.mod <- '
INTEREST       ~ a*conditionbi
attitude_mean  ~ b*conditionbi
expintent_avg  ~ c*conditionbi + d*attitude_mean

bd:= b*d
totmed:= bd + c'

med.fit <- bsem(med.mod, data = blav_data_out, target = "stan", seed = 2019)
med.out <- summary(med.fit, fit.measures = T)
med.pam <- parameterEstimates(med.fit)

med.table <- med.out[c(1:4, 14:15),] %>%
  as.tibble() %>%
  rename("Parameter" = `        `)

med.table <- rbind("", med.table[1:4,], "", med.table[5:6,]) %>%
         mutate(Parameter = c("Direct effects", "Interest ~ Condition", "Attitude ~ Condition", "Intent/Expectation ~ Condition", "Intent/Expectation ~ Attitude",
                 "Mediation", "Indirect effects/nCondition > Attitude > Intention/Expectation", "Total effect"))

##Is perception of future norm a mediator of the effect of trending minority norms vs. minority only on interest?
medpercept.mod <- '
expintent_avg       ~ a*conditionbi
attitude_mean  ~ b*conditionbi
PERCEPTCHANGE  ~ c*conditionbi
INTEREST  ~ d*conditionbi + e*PERCEPTCHANGE

ce:= c*e
totmed:= ce + d'

medpercept.fit <- bsem(medpercept.mod, data = blav_data_out, target = "stan", seed = 2019)
medpercept.out <- summary(medpercept.fit, fit.measures = T)
medpercept.pam <- parameterEstimates(medpercept.fit)

# Moderation ---------------------------
attitude_mod.out <- summary(lm(attitude_mean ~ conditionbi + genderbi + POLITICS + conditionbi:genderbi + conditionbi:POLITICS, no_out)) %>%
  apa_print()

intent_mod.out <- summary(lm(expintent_avg ~ conditionbi + genderbi + POLITICS + conditionbi:genderbi + conditionbi:POLITICS, no_out)) %>%
  apa_print()


save(h1.pwr, med.table, medpercept.out, intent_mod.out, attitude_mod.out, file =  "data/unregistered_out.RData")
