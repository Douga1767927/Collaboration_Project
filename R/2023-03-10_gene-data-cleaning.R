# Load Libs
pacman::p_load(tidyverse, here)

# Read in data
data <- read_csv(
  here::here(
    "data", "2023-03-22_gene-exp_cleaned.csv"
  )
)
data

## Cleaning ##

# Remove missing value
data <- na.omit(data)

# Getting Cell type from block Id
data$cell_type <- as.character(data$block)
data$cell_type[data$cell_type == "1"] <- "GL-CsE"
data$cell_type[data$cell_type == "2"] <- "GL-bNo"
data$cell_type[data$cell_type == "3"] <- "GL-JZC"
data$cell_type[data$cell_type == "4"] <- "GL-fUg"
data$cell_type[data$cell_type == "5"] <- "GL-jEK"
data$cell_type[data$cell_type == "6"] <- "GL-Hoe"
data$cell_type[data$cell_type == "7"] <- "GL-Rza"
data$cell_type[data$cell_type == "8"] <- "GL-xpo"

# Creating new variables "treatment" and "cell_type"
data <- data |>
  mutate(
    treatment = str_replace(group, "treatment", "AF42"),
    treatment = as.factor(treatment),
    cell_line = as.factor(cell_line),
    cell_type = as.factor(cell_type)
  )

# obtaining only necessary columns
data <- data |>
  select(gene_exp, concentration, cell_line, cell_type, treatment)

# Save cleaned dataset
data
write.csv(data, file = here::here("data","2023-06-01_cleaned-data-final.csv"), row.names=FALSE)
