#### INITIAL CODE ####
rm(list=ls())

#load in packages
library(here)
library(tidyverse)
library(reshape)
library(MASS)
library(emmeans)
library(DHARMa)
library(lme4)
library(glmmTMB)
library(car)

#read in data
gam_ooc <- read.csv(here("Data", "gametocyst_oocyst.csv"), header = TRUE)
dev <- read.csv(here("Data", "development.csv"), header = TRUE)

full <- read.csv(here("Data", "transmission_experiment.csv"), header = TRUE)
#for the treatment column, the letters correspond to treatments as follows:
#A -- 100% Colony A
#B -- 75% Colony A
#C -- 50% Colony A
#D -- 25% Colony A

#### TRANSMISSION EXPERIMENT DATA RE-STRUCTURING ####
#adjust data structure
full <- full %>%
  mutate(treat_rep = paste0(treatment, replicate_no)) %>% #add column combining treatment number and replicate number
  mutate(across(c(treatment, replicate_no, dissection_no, treat_rep), as.factor))

#remove dead, pupated, bad dissections
full <- full %>%
  filter(!trophont %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(trophont, gamont, gut_gametocyst, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(no_gut_parasites = trophont + (2*gamont) + (2*gut_gametocyst)) %>% #calculate total number of gut parasites
  mutate(no_gut_parasites = as.integer(no_gut_parasites))

#get parasite intensity data (only individuals with > 0 parasites)
int <- full %>%
  filter(no_gut_parasites > 0)

#### GAMETOCYST/OOCYST - DATA RE-STRUCTURING ####
gam_ooc <- gam_ooc %>%
  filter(!oocysts %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(oocysts, gametocysts, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(across(c(colony, experimental_repeat, day_removed), as.factor)) %>% #other columns into factors for GLMM
  mutate(oocysts = as.integer(400*oocysts)) %>% #multiply oocyst count by 400 to get number of oocysts in entire sample
  filter(oocysts > 0) #only keep individuals who produced oocysts

#### GAMETOCYST/OOCYST - GLMM ####
gam_mod <- glmer(oocysts ~ gametocysts + colony + experimental_repeat + (1|day_removed) + colony*gametocysts, data = gam_ooc,
                 family = Gamma(link = "log"))
summary(gam_mod)
#AIC 1002

#check fit
gam_output <- simulateResiduals(fittedModel = gam_mod)
plot(gam_output)

#posthoc test for overall difference between colonies
gam_posthoc <- emmeans(gam_mod, specs = pairwise ~ colony, weights = "proportional")
summary(gam_posthoc)
#no significant difference between colonies

#### DEVELOPMENT - DATA RE-STRUCTURING####
dev <- dev %>%
  mutate(across(c(colony, group), as.factor))

#### DEVELOPMENT - GLM ####
dev_mod <- glm(development ~ day + day*group + colony + colony*group + experimental_repeat, data = dev, family = "binomial")
summary(dev_mod)
#AIC 1589 

#check fit
dev_output <- simulateResiduals(fittedModel = dev_mod)
plot(dev_output)

#posthoc test for differences in development between exposed and unexposed larvae
dev_posthoc <- emmeans(dev_mod, specs = pairwise ~ group, weights = "proportional")
summary(dev_posthoc)
#development rate faster in unexposed than exposed

#### TEST WHETHER VARIANCE IN SUSCEPTIBILITY & INFECTIOUSNESS DIFFERED AMONG COLONIES ####
#read in data
susc <- read.csv(here("Data", "susceptibility.csv"), header = TRUE)
inf <- read.csv(here("Data", "infectiousness.csv"), header = TRUE)

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

#get variance and mean data
data_agg <- melt(susc.exp, id.vars = c("id", "day_dissected", "colony", "experimental_repeat"), measure.vars = c("no_gut_parasites"))
data_agg <- data_agg[which(data_agg$value > 0),]
var_mean <- as.data.frame(aggregate(value ~ colony + experimental_repeat, data_agg, function(x) c(mean = mean(x), var = var(x))))
var_mean <- do.call(data.frame, var_mean)
names(var_mean)[3:4] <- c("mean", "var")
var_mean$ratio <- var_mean$var/var_mean$mean

#test for differences in variance between colonies
leveneTest(no_gut_parasites ~ colony, data = susc_int)
#significant, so Colony A has a lower variance than Colony B

#data manipulation for infectiousness data
inf <- inf %>%
  filter(!oocysts %in% c("dead", "bad_dissection", "pupated")) %>%
  mutate(across(c(oocysts, par_presence), as.integer)) %>% #turn columns into integers, rather than characters
  mutate(across(c(colony, experimental_repeat, day_removed), as.factor)) %>% #other columns into factors for GLMM
  mutate(oocysts = as.integer(400*oocysts)) %>% #multiply oocyst count by 400 to get number of oocysts in entire sample
  filter(oocysts > 0) #only keep individuals who produced oocysts

#get variance and mean data
data_agg_inf <- melt(inf, id.vars = c("id", "day_removed", "colony", "experimental_repeat"), measure.vars = c("oocysts"))
var_mean_inf <- as.data.frame(aggregate(value ~ colony + experimental_repeat, data_agg_inf, function(x) c(mean = mean(x), var = var(x))))
var_mean_inf <- do.call(data.frame, var_mean_inf)
names(var_mean_inf)[3:4] <- c("mean", "var")
var_mean_inf$ratio <- var_mean_inf$var/var_mean_inf$mean

#test for differences in variance between colonies
leveneTest(oocysts ~ colony, data = inf)
#significant, so Colony A has a lower variance than Colony B

#### TEST FOR WEEK*TREATMENT INTERACTION TERM ####

##PREVALENCE
#first the original prevalence model, with no week*treatment interaction term
m0.prev <- glmmTMB(par_presence ~ dissection_no + treatment + (1 | treat_rep),
                   family = binomial,
                   data = full)

#then the prevalence model with a week*treatment interaction term
m1.prev <- glmmTMB(par_presence ~ dissection_no + treatment + (1 | treat_rep) + dissection_no*treatment,
                   family = binomial,
                   data = full)

#then compare model fits using a likelihood ratio test
anova(m0.prev, m1.prev)
#so the model without the interaction term is significantly better than the one with it

##INTENSITY
#first the original intensity model, with no week*treatment interaction term
m0.int <- glmmTMB(no_gut_parasites ~ dissection_no + treatment + (1 | treat_rep),
              family = nbinom2,
              data = int)

#then the intensity model with a week*treatment interaction term
m1.int <- glmmTMB(no_gut_parasites ~ dissection_no + treatment + (1 | treat_rep) + dissection_no*treatment,
              family = nbinom2,
              data = int)

#then compare model fits using a likelihood ratio test
anova(m0.int, m1.int)
#so the model without the interaction term is significantly better than the one with it

#### TEST WHETHER VARIANCE IN PREVALENCE IN TRANSMISSION EXPERIMENT IS DIFFERENT AMONG TREATMENTS####
#first aggregate data
prev_agg <- full %>%
  group_by(treatment, dissection_no, treat_rep) %>%
  summarise(
    prevalence = mean(par_presence, na.rm = TRUE),
    .groups = "drop"
  )

#then test for overall homogeneity of variance
leveneTest(prevalence ~ treatment, data = prev_agg, center = median)

#then test for homogeneity of variance when split by week
leveneTest(prevalence ~ treatment * dissection_no, data = prev_agg, center = median)

#### TEST WHETHER VARIANCE IN INTENSITY IN TRANSMISSION EXPERIMENT IS DIFFERENT AMONG TREATMENTS####
#first aggregate data
int_agg <- full %>%
  filter(no_gut_parasites > 0) %>%
  group_by(treatment, dissection_no, treat_rep) %>%
  summarise(mean_int = mean(no_gut_parasites, na.rm = TRUE),
            n_inf = n(),
            .groups = "drop")

#then test for overall homogeneity of variance
leveneTest(mean_int ~ treatment, data = int_agg, center = median)

#then test for homogeneity of variance when split by week
leveneTest(mean_int ~ treatment * dissection_no, data = int_agg, center = median)


##Additional method (not used in SI)
# #first test the original intensity model, where there is the same variance for all treatments
# m0.int <- glmmTMB(no_gut_parasites ~ dissection_no + treatment + (1 | treat_rep),
#                   family = nbinom2, #same mean-variance relationship as glmer.nb
#                   data = int)
# 
# #then test a treatment-specific dispersion model (variance changes with treatment)
# m1.int <- glmmTMB(no_gut_parasites ~ dissection_no + treatment + (1 | treat_rep),
#                   family = nbinom2,
#                   dispformula = ~ 0 + treatment,  #allows a separate dispersion per treatment
#                   data = int)
# 
# #then compare model fits
# anova(m0.int, m1.int)


#### SIMULATION BASED POWER CALCULATION - PREVALENCE ####
#set seed for replicability
set.seed(2025)

#Experimental setup
n_treat <- 5 #number of treatments
n_rep_per_treat <- 5 #replicates per treatment
N_pop <- 40 #individuals per replicate/population
nsim <- 2000  #number simulations
alpha <- 0.05 #singificance threshold

#treatment proportions of Colony A
propA <- c(1.0, 0.75, 0.5, 0.25, 0.0)

#colony-level prevalence
p_A <- 0.26 #colony A
p_B <- 0.42 #colony B

#Simulation function
sim_one_prevalence <- function(){
  dat_list <- list()
  id <- 1
  for(t in seq_along(propA)){
    for(rep in 1:n_rep_per_treat){
      #assign colony to N_pop individuals
      nA <- round(N_pop * propA[t])
      nB <- N_pop - nA
      colony <- c(rep("A", nA), rep("B", nB))
      #infection draws
      prob <- ifelse(colony == "A", p_A, p_B)
      infected <- rbinom(length(prob), size = 1, prob = prob)
      df <- data.frame(treatment = factor(paste0("T", t), levels = paste0("T", 1:n_treat)),
                       rep = factor(paste0("T", t, "_R", rep)),
                       colony = colony,
                       infected = infected)
      dat_list[[id]] <- df
      id <- id + 1
    }
  }
  data_sim <- bind_rows(dat_list)
  #fit GLMM: infected ~ treatment + (1|rep)
  #(use binomial family, individual-level)
  m <- tryCatch(glmmTMB(infected ~ treatment + (1|rep),
                        data = data_sim, family = binomial),
                error = function(e) NULL)
  if(is.null(m)) return(NA)
  #test overall treatment effect (Type II or likleihood ratio test): do anova comparing to null model
  m0 <- update(m, . ~ . - treatment)
  lrt <- anova(m0, m)
  pval <- lrt$`Pr(>Chisq)`[2]   #p for treatment
  return(pval < alpha)
}

#Run simulations
system.time({
  results <- replicate(nsim, sim_one_prevalence())
})

power_est <- mean(results, na.rm = TRUE)

cat("Estimated power for prevalence test (GLMM):", round(power_est,3), "\n")


#### SIMULATION BASED POWER CALCULATION - INTENSITY ####
#can use many of the same variables as the prevalence power estimate
#colony-level mean intensity among infected (mu) and NB size parameter (dispersion)
mu_A <- 29.6
mu_B <- 76.3

#size parameter
size_A <- 0.5573
size_B <- 0.5573

#Simulation function
sim_one_intensity <- function(test_dispersion = FALSE){
  dat_list <- list()
  id <- 1
  for(t in seq_along(propA)){
    for(rep in 1:n_rep_per_treat){
      nA <- round(N_pop * propA[t])
      nB <- N_pop - nA
      colony <- c(rep("A", nA), rep("B", nB))
      pcol <- ifelse(colony == "A", p_A, p_B)
      #infection status
      infected <- rbinom(length(pcol), 1, pcol)
      #counts: infected have NB with mean mu and size param
      mu <- ifelse(colony == "A", mu_A, mu_B)
      size <- ifelse(colony == "A", size_A, size_B)
      counts <- ifelse(infected == 1, rnbinom(sum(infected==1), mu = mu[infected==1], size = size[infected==1]), 0)
      counts_full <- numeric(length(colony))
      for(i in seq_along(colony)){
        if(infected[i]==1){
          mm <- ifelse(colony[i]=="A", mu_A, mu_B)
          ss <- ifelse(colony[i]=="A", size_A, size_B)
          counts_full[i] <- rnbinom(1, mu = mm, size = ss)
        } else counts_full[i] <- 0
      }
      df <- data.frame(treatment = factor(paste0("T", t), levels = paste0("T", 1:n_treat)),
                       rep = factor(paste0("T", t, "_R", rep)),
                       colony = colony,
                       count = counts_full)
      dat_list[[id]] <- df
      id <- id + 1
    }
  }
  data_sim <- bind_rows(dat_list)
  #fit NB GLMM: count ~ treatment + (1|rep)
  m <- tryCatch(glmmTMB(count ~ treatment + (1|rep), data = data_sim, family = nbinom2),
                error = function(e) NULL)
  if(is.null(m)) return(NA)
  #test overall treatment effect (LRT)
  m0 <- update(m, . ~ . - treatment)
  lrt <- anova(m0, m)
  p_treat <- lrt$`Pr(>Chisq)`[2]
  sig_treat <- (p_treat < alpha)
  sig_disp <- NA
  if(test_dispersion){
    #compare dispformula
    m_disp <- tryCatch(glmmTMB(count ~ treatment + (1|rep), dispformula = ~ treatment, data = data_sim, family = nbinom2),
                       error = function(e) NULL)
    if(!is.null(m_disp)){
      cmp <- anova(m, m_disp)
      #the anova ordering may differ; get p-value for disp model
      p_disp <- cmp$`Pr(>Chisq)`[2]
      sig_disp <- (p_disp < alpha)
    }
  }
  return(c(sig_treat = sig_treat, sig_disp = sig_disp))
}

#Run simulations
res <- replicate(nsim, sim_one_intensity(test_dispersion = TRUE))
res_mat <- t(res)
power_treat <- mean(res_mat[,"sig_treat"], na.rm = TRUE)
power_disp <- mean(res_mat[,"sig_disp"], na.rm = TRUE)

cat("Power (detect treatment effect on intensity):", round(power_treat,3), "\n")
cat("Power (detect dispersion differences):", round(power_disp,3), "\n")
