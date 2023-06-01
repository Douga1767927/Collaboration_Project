## Stats Consulting - Collaboration Project (Task 4) - IMRaD Code
## Douglas Dally - a1767927
## 29/05/2023

## Load Libs
pacman::p_load(tidyverse, lmerTest, MuMIn, gt, gglm)

## load data
data <- read.csv(here::here("data", "2023-05-29_cleaned-data-final.csv"))

# clean data
data <- data |>
  mutate(cell_type = as.factor(exp_name),
         cell_line = as.factor(cell_line),
         treatment = as.factor(group))

# obtain necesarry columns
data <- data |>
  select(cell_line, gene_exp, concentration, treatment, cell_type)
data

# Full fixed effect model
M_full <- lmer(gene_exp ~ concentration*treatment*cell_line + (1|cell_type), data = data)
summary(M_full)$coefficients
R2_full <- r.squaredGLMM(M_full)

# remove 3 way interaction term
M2 <- update(M_full, .~. - concentration:treatment:cell_line)
summary(M2)$coefficients
R2_2 <- r.squaredGLMM(M2)

# remove 2 way interaction term with insignificant P-value
M3 <- update(M2, .~. - treatment:cell_line)
summary(M3)$coefficients
R2_3 <- r.squaredGLMM(M3)

# still have insignificant terms as single predictors. hence to remove them we must remove all interaction terms
M4 <- update(M3, .~. - concentration:treatment - concentration:cell_line)
summary(M4)$coefficients
R2_4 <- r.squaredGLMM(M4)

# cell_line wild is insignificant
M5 <- update(M4, .~. - cell_line)
summary(M5)$coefficients
R2_5 <- r.squaredGLMM(M5)

# get summary table for 5 models
AIC <- anova(M5, M4, M3, M2, M_full)[2]
R_squared <- rbind(R2_5, R2_4, R2_3, R2_2, R2_full)
n_pred <- c(2,3,5,6,7)

sum_tab <- cbind(n_pred, AIC, R_squared)
sum_tab

# Predicting data
new_data <- data.frame(
  cell_line = "wild",
  cell_type = "GL-bNo",
  concentration = 8,
  treatment = "placebo")
new_data

prediction <- predict(M_full, new_data)
actual <- 8.23
tab2 <- cbind(new_data, prediction, actual)
tab2 |> gt()


