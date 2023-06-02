## Stats Consulting - Collaboration Project (Task 4) - IMRaD Code
## Douglas Dally - a1767927
## 29/05/2023

## Load Libs
install.packages("pacman")
pacman::p_load(tidyverse, lmerTest, MuMIn, gt, gglm, reshape)

## load data
data <- read.csv(here::here("data", "2023-06-01_cleaned-data-final.csv"))
data

# Full fixed effect model
M_full <- lmer(gene_exp ~ concentration*treatment*cell_line + (1|cell_type), data = data)
anova(M_full)
R2_full <- r.squaredGLMM(M_full)

# remove 3 way interaction term
M2 <- update(M_full, .~. - concentration:treatment:cell_line)
anova(M2)
R2_2 <- r.squaredGLMM(M2)

# remove 2 way interaction term with insignificant P-value
M3 <- update(M2, .~. - treatment:cell_line)
anova(M3)
R2_3 <- r.squaredGLMM(M3)

# still have insignificant terms as single predictors. hence to remove them we must remove all interaction terms
M4 <- update(M3, .~. - concentration:treatment - concentration:cell_line)
anova(M4)
R2_4 <- r.squaredGLMM(M4)

# cell_line wild is insignificant
M5 <- update(M4, .~. - cell_line)
anova(M5)
R2_5 <- r.squaredGLMM(M5)

M_anova <- anova(M_full, M2, M3, M4, M5)

# get summary table for 5 models
AIC <- M_anova$AIC
R_squared <- rbind(R2_5, R2_4, R2_3, R2_2, R2_full)
n_fixed_effects <- c(2,3,5,6,7)
model <- c("step 5", "step 4", "step 3", "step 2", "full")

sum_tab <- cbind(model, n_fixed_effects, AIC, R_squared)
sum_tab <- as.data.frame(sum_tab)

# save sum_tab as dataframe to data folder
write.csv(sum_tab, file = here::here("data","2023-06-01_models-summary.csv"), row.names=FALSE)

sum_tab <- sum_tab |>
  gt() |>
  tab_header(
    title = "Summary of Backward Selection from Full Model")

#save table as .docx file
sum_tab |> gtsave(
      here::here("tabs", "2023-06-01_model-selection-summary.docx")
    )

## Predicting data for wild type GL-bNo cells with placebo treatment

# getting test data
test_data <- data |>
  filter(cell_type == "GL-bNo") |>
  select(-gene_exp)
test_data

# getting the actual gene expression values
test_actual_gene_exp <- data |>
  filter(cell_type == "GL-bNo") |>
  select(gene_exp)

# predict values for gene expression based on the test data
prediction <- predict(M_full, test_data)

# get a prediction-actual dataset so we can plot predicted and actual gene expression
pred_df <- cbind(test_data, prediction, test_actual_gene_exp) |>
  #melt dataframe to combine predicted and actual values into one column
  melt(measure.vars = c("prediction", "gene_exp"))
pred_df

# plot gene expression vs concentration for predicted and actual values
pred_plot <- pred_df |>
  ggplot(aes(x = concentration, y = value)) +
  geom_point(aes(col = variable), size = 3) +
  labs(title = "Prediction vs Actual for Wild GL-bNo cell with Placebo",
       col = "Variable",
       x = expression(paste("Growth Factor Concentration (", mu,"g/ml)")),
       y = "Gene Expression") +
  scale_color_manual(values = c("orange", "darkblue"), labels = c("Predicted", "Actual")) +
  theme_bw()

# save the prediction-actual plot
ggsave(
  filename = here::here("figs", "2023-06-01_prediction-actual.png"),
  plot = pred_plot,
  width = 12,
  height = 9
)


