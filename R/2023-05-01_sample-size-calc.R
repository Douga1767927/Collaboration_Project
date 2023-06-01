## Stats Consulting - Collaboration Project (Task 3)
## Douglas Dally - a1767927
## 04/04/2023

## Libs
pacman::p_load(pwr)

## Sample Size Calculation

## known parameters
power = 0.9  # required power
beta = 1 - power
alpha = 0.05 # significance level
r2 <- 0.1 # r squared value
u <- 5 # number of indpendent predictors

## Calculate effect size
f2 <- r2/(1-r2)

## Get power calculation
pwr_calc <- pwr.f2.test(u = k, f2 = f2, power = power)

## sample size required
n <- round(k + pwr_calc[[2]] + 1)
n
