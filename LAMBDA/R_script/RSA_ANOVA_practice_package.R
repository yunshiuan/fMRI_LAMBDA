#RSA: ANOVA for relatedness test results===========================================
library(dplyr)
library(ez)
library(WRS2)
#ezANOVA---------------------------------------------------------------------------
#Problem: This is a traditional parametric ANOVA.
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)
ezPrecis(ANT)

#Run an ANOVA on the mean correct RT data.
rt_anova = ezANOVA(
  data = ANT[ANT$error==0,]
  , dv = .(rt)
  , wid = .(subnum)
  , within = .(cue,flank)
  , between = .(group)
)

#Show the ANOVA & assumption tests.
print(rt_anova)

#Run an ANOVA on the mean_rt data, ignoring group.
rt_anova2 = ezANOVA(
  data = ANT[ANT$error==0,]
  , dv = .(rt)
  , wid = .(subnum)
  , within = .(cue,flank)
)

#Show the ANOVA & assumption tests.
print(rt_anova2)

#Run a purely between-Ss ANOVA on the mean_rt data.
rt_anova3 = ezANOVA(
  data = ANT[ANT$error==0,]
  , dv = .(rt)
  , wid = .(subnum)
  , between = .(group)
)

#Show the ANOVA & assumption tests.
print(rt_anova3)

#Greenhouse-Geisser and Huynh-Feldt are given Multiply dfn and dfd by the correction estimate
#Post-hoc tests:
#In this case, we should use the Bonferroni correction (because we have a lack of sphericity!)
#pairwise.t.test(bushLong$retch, bushLong$animal,
#                paired=T, p.adjust.method="bonferroni") 
#Robust version of ANOVA----------------------------------------------------------------------------
#Problem: Only support up to 2-way mixed ANOVA
#The WRS2 Package: https://cran.r-project.org/web/packages/WRS2/vignettes/WRS2.pdf
#Robust measures of location: (p.16 for mixed design)
#(1)trimmed mean which discards a certain percentage at both ends of the distribution. 
#(2)Winsorized mean: giving less weight to observations in the tails of the distribution and higher weight to the ones in the center
#(3)M-estimators (the "M" stands for "maximum likelihood-type", use bootstrap based functions)

#(Trimmed mean)fit the between-within subjects ANOVA on the 20% trimmed means (defaults is 20%)
bwtrim(symptoms ~ group*time, id = id, data = hangover)

#we fit a standard between-within subjects ANOVA through bwtrim by 
#setting the trimming level to 0 and check whether we get the same results as with ezANOVA
bwtrim(symptoms ~ group*time, id = id, data = hangover, tr = 0)
hangover$id <- as.factor(hangover$id)
fitF <- ezANOVA(hangover, symptoms, between = group, within = time, wid = id)
fitF$ANOVA

#(M-estimators)base our comparisons on Huber's M-estimator (bootstrap based functions)
#(for which we have to apply three separate functions, one for each effect.)
sppba(symptoms ~ group*time, id, data = hangover)
sppbb(symptoms ~ group*time, id, data = hangover)
sppbi(symptoms ~ group*time, id, data = hangover)




#Nonparametric Analysis of Longitudinal Data in Factorial Experiments---------------------------------------------
#Problem: It assumes a "time" variable
#nonparametric rank-based methods: enables to analyze ordered categorical, dichotomous, and heavily skewed data
#https://www.jstatsoft.org/article/view/v050i12
#MCMCglmm---------------------------------------------------------------------------------------------------------