#### INITIAL CODE ####
rm(list=ls())

#load in packages
library(here)
library(tidyverse)
library(MASS)
library(emmeans)
library(DHARMa)
library(lme4)

#read in data
susc <- read.csv(here("Data", "susceptibility.csv"), header = TRUE)
inf <- read.csv(here("Data", "infectiousness.csv"), header = TRUE)

#### SUSCEPTIBILITY - DATA RE-STRUCTURING ####
#remove dead, pupated, bad dissections
susc <- susc %>%
  filter(!gamont %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(trophont, gamont, gut_gametocyst, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(no_gut_parasites = trophont + (2*gamont) + (2*gut_gametocyst)) %>% #calculate total number of gut parasites
  mutate(no_gut_parasites = as.integer(no_gut_parasites)) %>%
  mutate(across(c(colony, group, experimental_repeat, day_removed), as.factor))

#get only experimental_repeat individuals
susc.exp <- susc %>%
  filter(group == "experimental")

#get prevalence and intensity data sets
susc_prev <- susc.exp
susc_int <- susc.exp %>%
  filter(no_gut_parasites > 0)

#### SUSCEPTIBILITY - PREVALENCE GLMM ####
susc_prev_mod <- glmer(par_presence ~ colony + experimental_repeat + (1|day_removed), data = susc_prev, family = "binomial")
summary(susc_prev_mod)
#AIC 724.6

#check fit
susc_prev_output <- simulateResiduals(fittedModel = susc_prev_mod)
plot(susc_prev_output)

#posthoc test for overall difference between colonies
susc_prev_posthoc <- emmeans(susc_prev_mod, specs = pairwise ~ colony, weights = "proportional")
summary(susc_prev_posthoc)
#significant difference between colonies

#get estimates for colony mean prevalence
susc_prev_mean <- as.data.frame(summary(susc_prev_posthoc$emmeans, type = "response"))

#test for differences between experimental repeats
susc_prev_posthoc2 <- emmeans(susc_prev_mod, specs = pairwise ~ experimental_repeat, weights = "proportional")
summary(susc_prev_posthoc2)
#no differences

#get estimates for colony mean prevalence per experimental repeat
susc_prev_posthoc3 <- emmeans(susc_prev_mod, specs = pairwise ~ colony + experimental_repeat, weights = "proportional")
susc_prev_mean_repeat <- as.data.frame(summary(susc_prev_posthoc3$emmeans, type = "response"))

#### SUSCEPTIBILITY - INTENSITY GLMM ####
susc_int_mod <- glmer.nb(no_gut_parasites ~ colony + experimental_repeat + (1|day_removed), data = susc_int)
summary(susc_int_mod)
#AIC 1934.2

#check fit
susc_int_output <- simulateResiduals(fittedModel = susc_int_mod)
plot(susc_int_output)

#posthoc test for overall difference between colonies
susc_int_posthoc <- emmeans(susc_int_mod, specs = pairwise ~ colony, weights = "proportional")
summary(susc_int_posthoc)
#significant difference between colonies

#get estimates for colony mean intensity
susc_int_mean <- as.data.frame(summary(susc_int_posthoc$emmeans, type = "response"))

#test for differences between experimental repeats
susc_int_posthoc2 <- emmeans(susc_int_mod, specs = pairwise ~ experimental_repeat, weights = "proportional")
summary(susc_int_posthoc2)
#significant difference between 1 and 2

#### INFECTIOUSNESS - DATA RE-STRUCTURING ####
inf <- inf %>%
  filter(!oocysts %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(oocysts, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(across(c(colony, experimental_repeat, day_removed), as.factor)) %>% #other columns into factors for GLMM
  mutate(oocysts = as.integer(400*oocysts)) %>% #multiply oocyst count by 400 to get number of oocysts in entire sample
  filter(oocysts > 0) #only keep individuals who produced oocysts
  
inf <- inf %>% 
  mutate(OLRE = 1:nrow(inf)) #add an observation-level random effect column

#### INFECTIOUSNESS GLMM ####
inf_mod <- glmer(oocysts ~ colony + experimental_repeat + (1|day_removed) + (1|OLRE),
                 data = inf, family = "poisson",
                 control = glmerControl(optimizer = "bobyqa"))
summary(inf_mod)
#AIC 7497.9

#check fit
inf_output <- simulateResiduals(fittedModel = inf_mod)
plot(inf_output)

#posthoc test for overall difference between colonies
inf_posthoc <- emmeans(inf_mod, specs = pairwise ~ colony, weights = "proportional")
summary(inf_posthoc)
#significant difference between colonies

#get estimates for colony mean infectiousness
inf_mean <- as.data.frame(summary(inf_posthoc$emmeans, type = "response"))

#post hoc test for differences in treatment
inf_posthoc2 <- emmeans(inf_mod, specs = pairwise ~ experimental_repeat, weights = "proportional")
summary(inf_posthoc2)
#there are significant differences
