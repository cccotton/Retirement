cat("\014")
#Dplyr approach - Cary

#Functions
summaryFailureYear <- function (x) {
  balance <- portfolio
  
  runningbalance <- vector(mode = "numeric", length = 0)
  
  for (k in 1:50) {
    balance <- (balance - spending) * x[k]
    runningbalance[k] <- balance
    if (balance < 0) break
  }
  if(k == 50) k <- 99
  return(k)
}

##Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(survival)
library(cmprsk)

##Files and formatting
marketReturns <- read.csv("Random Market Returns Log-N 10Kx50.csv",header = FALSE)
lifeSpanMale <- read.csv("Male Random Lifetimes from Age 65.csv",header=FALSE)
lifeSpanFemale <- read.csv("Female Random Lifetimes from Age 65.csv",header=FALSE)
lifeSpanJoint <- read.csv("Joint Random Lifetimes from Age 65.csv",header=FALSE)

#marketReturns
marketReturns.tdf <- tbl_df(marketReturns) %>%
  mutate(series = row_number())
cleanedMarketReturns <- marketReturns.tdf %>%
  gather(year, return.annual, -series, convert = T) %>%
  mutate(year = as.numeric(sub("V", "", year))) %>%
  arrange(series, year) %>%
  group_by(series)
cleanedMarketReturns

#lifeSpanMale
lifeSpanMale.tdf <- tbl_df(lifeSpanMale)

cleanedlifeSpanMale <- lifeSpanMale.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanMale

#lifeSpanFemale
lifeSpanFemale.tdf <- tbl_df(lifeSpanFemale)

cleanedlifeSpanFemale <- lifeSpanFemale.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanFemale

#lifeSpanJoint
lifeSpanJoint.tdf <- tbl_df(lifeSpanJoint)

cleanedlifeSpanJoint <- lifeSpanJoint.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanJoint

##Plotting distributions
weight.inverse <- 1 / (max(cleanedMarketReturns$series) * max(cleanedMarketReturns$year))
ggplot(data = cleanedMarketReturns, aes(x = return.annual)) +
  geom_bar(aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Rate of annual returns") + ylab("Proportion of annual returns") +
  theme_classic()

weight.inverse <- 1 / 10000
ggplot(data = cleanedlifeSpanMale, aes(x = death.years)) + ggtitle("Male") +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") +
  theme_classic()

ggplot(data = cleanedlifeSpanFemale, aes(x = death.years)) +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") + ggtitle("Female") +
  theme_classic()

ggplot(data = cleanedlifeSpanJoint, aes(x = death.years)) +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") + ggtitle("Joint") +
  theme_classic()

##Parameters
portfolio <- 1000000
wr = 0.05 # percentage of initial portfolio value to spend annually
spending <- wr * portfolio

##Failure years (conditional on above parameters)
solvedMarketScenarios <- cleanedMarketReturns %>%
  summarise(failure.year = summaryFailureYear(return.annual) + 65)

weight.inverse <- 1 / (max(solvedMarketScenarios$series))
ggplot(data = solvedMarketScenarios, aes(x = failure.year)) +
  geom_bar(aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Year of Portfolio Failure") + ylab("Proportion of cases (%)") +
  theme_classic()

######Run an indexed set of n iterations of the simulation########
n <- 5000
#######Query the CI estimations from year.start to year.end#######
year.start <- 66
year.end <- 100
CISpan <- seq(year.start, year.end)

simulation <- tbl_df(data.frame(sim.index = as.numeric(rep(1:n, each=1000)), cohort.member = rep(as.numeric(seq(1, 1000)), n))) %>%
  bind_cols(sample_n(cleanedlifeSpanJoint, n * 1000, replace = T)) %>%
  mutate(series = sample(1:10000, n * 1000, replace=T)) %>%
  left_join(solvedMarketScenarios, by = "series") %>%
  rowwise() %>%
  mutate(event.year = pmin(death.years, failure.year),
         event.type = ifelse(failure.year < death.years, 1, 0)) %>%
  group_by(sim.index)

times <- simulation %>%
  do(time = as.vector(survfit(Surv(as.numeric(.$event.year), .$event.type) ~ 1)$time)) %>%
  unnest(time)

survivals <- simulation %>%
  group_by(sim.index) %>%
  do(surv = as.vector(survfit(Surv(as.numeric(.$event.year), .$event.type) ~ 1)$surv)) %>%
  unnest(surv) %>%
  rename(index2 = sim.index)

combinedSurvResults <- times %>%
  bind_cols(survivals) %>%
  filter(index2 == sim.index)

simwisemeans.surv <- combinedSurvResults %>%
  group_by(time) %>%
  summarise(smean = mean(surv, na.rm = T),
            svar = var(surv, na.rm = T)) %>%
  ungroup()

cumulative.inc <- simulation %>%
  do(cifailure = timepoints(cuminc(as.numeric(.$event.year), .$event.type, cencode = 2), CISpan)) %>%
  mutate(cifunpack = cifailure[1]) %>%
  select(-cifailure) %>%
  do(sim.index = .$sim.index,
     failure = .$cifunpack[1, ], 
     death = .$cifunpack[2, ]) %>%
  mutate(sim.index = as.integer(sim.index))

failure <- cumulative.inc %>%
  select(sim.index, failure) %>%
  unnest(failure)

death <- cumulative.inc %>%
  select(sim.index, death) %>%
  unnest(death) %>%
  rename(index2 = sim.index)

combinedCIResults <- failure %>%
  bind_cols(death) %>%
  filter(index2 == sim.index)

simwisemeans.ci <- combinedCIResults %>%
  mutate(time = rep(CISpan, n)) %>%
  group_by(time) %>%
  summarise(cifmean = mean(failure),
            cidmean = mean(death),
            cifvar = var(failure),
            cidvar = var(death))

##Graph the mean curves with simulation standard error 95% CL
ggplot(data = simwisemeans.surv, aes(x = time, y = smean)) +
  geom_ribbon(aes(ymin = (smean - 1.96 * sqrt(svar)), ymax = (smean + 1.96 * sqrt(svar))), alpha = 0.5) +
  geom_path() +
  xlab("Age") + ylab("Proportion 1000-Member Cohort Without Portfolio Failure") + xlim(65, 100) + ylim(0, 1) +
  theme_classic()

crsim <- cbind(gather(select(simwisemeans.ci, time, cidmean, cifmean), curve, mean, -time),
               select(gather(select(simwisemeans.ci, time, cidvar, cifvar), curve, var, -time), var))

crsim$curve <- factor(crsim$curve, labels = c("Portfolio failure", "Death"))

ggplot(data = crsim, aes(x = time, y = mean)) +
  geom_path(aes(linetype = curve)) +
  geom_ribbon(aes(ymin = (mean - 1.96 * sqrt(var)),
                  ymax = (mean + 1.96 * sqrt(var)),
                  linetype = curve),
              alpha = 0.05) +
  xlab("Age") + ylab("Proportion 1000-Member Cohort Specified Outcome") + xlim(65, 100) +
  scale_linetype_discrete(name = "") +
  theme_classic()
