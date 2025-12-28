library(cTOST)
data(skin)
diff <- skin$Generic - skin$Reference
hist(diff, probability = TRUE, breaks = 5)
#set probability = TRUE when we want to see the distribution density 
#set prob = FALSE when we want to see the raw count 
qqnorm(diff)
qqline(diff)
#qq plot compares quantile of observed data to quantiles of theoretical distribution to see if sample data follows a distribution

library(ggplot2)
ggplot(skin, aes(x = diff)) +
  geom_histogram(aes(y = ..density..), bins = 5) +
  geom_density(color = "red")
#prettier histogram

t.test(skin$Generic-skin$Reference, conf.level=0.9)
#t test to find what's the interval that will contain the true mean in 90% of repeated experiments
#H0 is true mean is 0; H1 is true mean is NOT 0.
# 90% CI (-0.205, 0.250), contains true mean H0
#p-value is 0.864: assuming true mean = 0, there is 86.4% chance of observing a value with |T| <= 0.174
#p-value is > 0.05, so fail to reject H0 and the result is consistent with 90% CI


theta_hat = diff(apply(skin,2,mean))
#the diff of means of two columns. 2 means apply the func across column, 1 means across rows
nu = nrow(skin) - 1
#find the degree of freedom
sig_hat = var(apply(skin,1,diff))/nu
#find the sample var of paired difference (summing squares of diff between each pair - mean diff), applying the func across rows; divide by nrows to get var of sample mean diff
stost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
              delta = log(1.25), method = "unadjusted")
stost
#standard TOST, test to see if the true mean (diff) is within margin log(1.25). 
#we calculate set alpha as 0.05 (FIXED alpha), and find the corresponding CI (1-2alpha), we check if CI lies within the margin--CI =  (-0.61707 ; 0.66248)--it DOES NOT, Can't accept bioequivalence

opt_tost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
                 delta = log(1.25), method = "optimal")
opt_tost
#optimal TOST, down-adjusted alpha (0.0385) so that the calculated CI is tighter and might fit within the acceptable margin 
#use opt_tost because CI sometimes might be too wide for small sample, so equivalence fail. So adjust CI to make it tighter around the mean 

###TOST
#two one-sided tests: traditional t-test looks for difference between two groups
  #H0 is u = 0, H1 is u NOT = 0
#tost looks if two groups are similar enough to be considered equivalence 
  #H1 is μ≤−Δ or μ≥Δ (not equivalent); H1 is -Δ < μ < Δ, delta is equivalence margin

###Q4-1
#A scientist knows that the mean height of females in England is 165cm and wants to know whether her patients with a certain disease “X” have heights that differ significantly from the population mean - we will use a one-sample t-test to test this.
#H0:sample mu = 165
#H1:sample mu NOT = 165
library(readr)
DiseaseX <- read_csv("DiseaseX.csv")

ggplot(DiseaseX, aes(Height)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, colour = "black", fill = "white") +
  stat_function(fun = dnorm, color = "red", args = list(mean = mean(DiseaseX$Height), sd = sd(DiseaseX$Height)))
qqnorm(DiseaseX$Height)
qqline(DiseaseX$Height)
#check if it's normally distributed
#stat_function(): adds a curve defined by a math function to a ggplot
#fun = dnorm: the function to draw → the normal (Gaussian) density
#color = "red": sets the curve color
#args = list(mean = ..., sd = ...): define parameters to draw the norm

t.test(DiseaseX$Height, mu = 165, alternative ="two.sided")
#this test calculates the sample mean and variance. Then it uses true mean u to calculate a T value, using this critical t value and df, we can use software to compute out a 95% CI assuming the test statistics follow t-distribution. If true mean is within the CI, we fail to reject H0. 
#alpha--using to define the CI--is the prob that "true mean" is outside CI when H0 is true. p-value is calculated from t value.
#p-value < 0.05 (true mean has very little chance to be outside) or CI coverage are equivalent, having either one is sufficient to say mu is likely not the true mean, hence reject H0
#in all other cases, we will conclude we fail to reject H0. 

#CI (169.34, 172,88) it does not contain true mean
#p-value = 7.293e-07
#these two consistently suggest that H0 is rejected, mean height differs from 165cm. 


###Q4-2
#reduction in bone marrow vessel density after 7 days of treatment
#H0 diff = 0 
#H1 diff < 0, one-tailed 
BloodV <- read_csv("BloodVesselFormation1.csv")
ggplot(BloodV, aes(Difference)) +
  geom_histogram(aes(y = ..density..), binwidth= 50, colour = "black", fill = "white") +
  stat_function(fun = dnorm, color = "red", args = list(mean = mean(BloodV$Difference), sd = sd(BloodV$Difference)))
qqnorm(BloodV$Difference)
qqline(BloodV$Difference)

t.test(BloodV$Difference, mu = 0, alternative ="less")
#mean diff = -86.143, t-value = -1.843, p-value = 0.057, CI (-inf, 4.709)
#cannot reject H0, because CI contains true mean, and p-value > alpha (prob that true mean is outside 95% CI)

###Q5-1
#diff in duration of a biological process between WT and geneKO cell
#H0 diff = 0; H1 diff is NOt 0
BioProcess <- read_csv("BiologicalProcessDurations.csv")
ggplot(BioProcess, aes(Time)) +
  geom_histogram(aes(y = ..density..), binwidth= 30, colour = "black", fill = "white") +
  facet_wrap(vars(Group)) + #split plots according to varieties of "Group"
  stat_function(fun = dnorm, color = "red", args = list(mean = mean(BioProcess$Time), sd = sd(BioProcess$Time)))
ggplot(BioProcess, aes(sample = Time)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(vars(Group))
t.test(skin$Generic-skin$Reference, conf.level=0.9)