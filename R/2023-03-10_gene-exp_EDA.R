## Stats Consulting - Collaboration Project
## Douglas Dally - a1767927
## 10/03/2023

# Load Libs
pacman::p_load(tidyverse, here, patchwork, viridis, gt, vtable)

# Read in data
data <- read_csv(
  here::here(
    "data", "2023-03-22_gene-exp_cleaned.csv"
  )
)
data

## EDA

# CLEANING ----
# Skim data
skimr::skim_without_charts(data)

# Remove missing value
data <- na.omit(data)

# EDIT: Changing group level treatment to AF42 & new variable exp_name
data <- data |>
  mutate(
    group = str_replace(group, "treatment", "AF42")
    # exp_name = str_replace(as.character(block),
    #                        "[12345678]",
    #                        c("GL-CsE",
    #                          "GL-bNo",
    #                          "GL-JZC",
    #                          "GL-fug",
    #                          "GL-jEK",
    #                          "GL-Hoe",
    #                          "GL-Rza",
    #                          "GL-xpo"))
  )

# Converting character variables to factor
data <- data |>
  mutate(cell_line = as.factor(cell_line),
         group = as.factor(group),
         block = as.factor(block))
data

# Getting experiment names from block Id
data$exp_name <- as.character(data$block)
data$exp_name[data$exp_name == "1"] <- "GL-CsE"
data$exp_name[data$exp_name == "2"] <- "GL-bNo"
data$exp_name[data$exp_name == "3"] <- "GL-JZC"
data$exp_name[data$exp_name == "4"] <- "GL-fUg"
data$exp_name[data$exp_name == "5"] <- "GL-jEK"
data$exp_name[data$exp_name == "6"] <- "GL-Hoe"
data$exp_name[data$exp_name == "7"] <- "GL-Rza"
data$exp_name[data$exp_name == "8"] <- "GL-xpo"

# Save cleaned dataset
data
write.csv(data, file = here::here("data","2023-05-29_cleaned-data-final.csv"), row.names=FALSE)

# FIGURES ----
### Histograms
# Histogram and scatterplot of response (gene_exp)
p1 <- ggplot(data, aes(x = gene_exp)) +
  geom_histogram(fill = "orange", colour = "black") +
  labs(title = "Histogram of Gene Expression") +
       xlab("Gene Expression") +
       ylab("Count") +
  theme_bw()
p1

# Histogram of Gene Expression for Treatment and Placebo Groups
p2 <- ggplot(data, aes(x = gene_exp)) +
  geom_histogram(fill = "green4", colour = "black") +
  labs(title = "Histogram of Gene Expression for Treatment and Placebo Groups") +
  xlab("Gene Expression") +
  ylab("Count") +
  theme_bw()+
  facet_wrap(vars(group))
p2

# Histogram of Gene Expression for Wild and Type-101 Cell Lines
p3 <- ggplot(data, aes(x = gene_exp)) +
  geom_histogram(fill = "cyan4", colour = "black") +
  labs(title = "Histogram of Gene Expression for Wild and Type-101 Cell Lines") +
  xlab("Gene Expression") +
  ylab("Count") +
  theme_bw()+
  facet_wrap(vars(cell_line))
p3

### Scatterplots
# Scatterplot of Gene Expression vs. Concentration
p4 <- ggplot(data, aes(x = concentration,
                       y = gene_exp,
                       col = group)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = group), method = "lm") +
  labs(title = "Scatterplot of Gene Expression vs. Concentration",
       subtitle = "Coloured by Treatment (Activating Factor 42) or Placebo",
       col = "Subject Group") +
  ylab("Gene Expression") +
  xlab("Growth Factor Concentration (mg/mL)") +
  theme_bw()
p4

# Scatterplot of Gene Expression vs Concentration (Both Types)
p5 <- ggplot(
  data = data,
  aes(x = concentration, y = gene_exp)) +
  geom_point(aes(col = group), size = 3) +
  geom_smooth(aes(col = group), method = "lm") +
  facet_wrap(vars(cell_line)) +
  labs(title = "Scatterplot of Gene Expression vs. Concentration ",
       subtitle = "Coloured by Treatment Group (Placebo or Activating Factor 42)",
       col = "Treatment Group") +
  ylab("Gene Expression") +
  xlab("Growth Factor Concentration (mg/mL)") +
  theme_bw()
p5

# Scatterplot of Gene Exp. vs Conc. (treatment)
p6 <- ggplot(
  data = data |> filter(group == "AF42"),
  aes(x = concentration, y = gene_exp)) +
  geom_point(aes(col = exp_name), size = 3) +
  geom_smooth(aes(col = exp_name), method = "lm") +
  facet_wrap(vars(cell_line)) +
  labs(title = "Scatterplot of Treatment Group (AF42) Gene Expression vs. Concentration ",
       subtitle = "Coloured by Experiment Name",
       col = "Experiment Name") +
  ylab("Gene Expression") +
  xlab("Growth Factor Concentration (mg/mL)") +
  theme_bw()
p6

# Scatterplot of Gene Exp. vs Conc. (Placebo)
p7 <- ggplot(
  data = data |> filter(group == "placebo"),
  aes(x = concentration, y = gene_exp)) +
  geom_point(aes(col = exp_name), size = 3) +
  geom_smooth(aes(col = exp_name), method = "lm") +
  facet_wrap(vars(cell_line)) +
  labs(title = "Scatterplot of Placebo Group Gene Expression vs. Concentration ",
       subtitle = "Coloured by Experiment Block",
       col = "Experiment Block") +
  ylab("Gene Expression") +
  xlab("Growth Factor Concentration (mg/mL)") +
  theme_bw() +
  scale_colour_viridis_d()
p7


### Boxplots
# Boxplot of Gene Expression by Treatment Group
p8 <- ggplot(
  data = data,
  aes(x = group, y = gene_exp, fill = group)) +
  geom_boxplot() +
  facet_wrap(vars(cell_line)) +
  labs(title = "Boxplot of Gene Expression by Treatment Group (Placebo or Activating Factor 42)",
       fill = "Treatment Group") +
  ylab("Gene Expression") +
  xlab("Treatment Group") +
  theme_bw()
p8

# Boxplot of Gene Expression by Treatment Group wrapped by cell_line
p8 <- ggplot(
  data = data,
  aes(x = group, y = gene_exp, fill = group)) +
  geom_boxplot() +
  facet_wrap(vars(cell_line)) +
  labs(title = "Boxplot of Gene Expression by Treatment Group (Placebo or Activating Factor 42)",
       fill = "Treatment Group") +
  ylab("Gene Expression") +
  xlab("Treatment Group") +
  theme_bw()
p8

# Look at cell101 individual experiments
p9 <- ggplot(
  data = data |> filter(cell_line == "cell101",
                        group == "AF42"),
  aes(x = group, y = gene_exp, fill = group)) +
  geom_boxplot() +
  facet_wrap(vars(exp_name)) +
  labs(title = "Boxplots of Gene Expression for Cell101-Type with AF42 Treatment",
       subtitle = "Grouped by Experiment Name (GL-Rza or GL-xpo)",
       fill = "Treatment Group") +
  ylab("Gene Expression") +
  xlab("Experiment Name") +
  theme_bw()
p9


# TABLES ----
options(digits = 2)
# Wild Table
wild_af<- data |>
  filter(cell_line == "wild", group == "AF42") |>
  select(gene_exp) |>
  pull()

wild_plac <- data |>
  filter(cell_line == "wild", group == "placebo") |>
  select(gene_exp) |>
  pull()

wild_data <- data |>
  filter(cell_line == "wild") |>
  select(group, gene_exp) |>
  group_by(group) |>
  skimr::skim_without_charts() |>
  select(-c(n_missing,
            complete_rate,
            skim_variable,
            skim_type)) |>
  mutate(iqr = c(IQR(wild_af), IQR(wild_plac)))
wild_data

# Cell101 Table
cell101_af<- data |>
  filter(cell_line == "cell101", group == "AF42") |>
  select(gene_exp) |>
  pull()

cell101_plac <- data |>
  filter(cell_line == "cell101", group == "placebo") |>
  select(gene_exp) |>
  pull()

cell101_data <- data |>
  filter(cell_line == "cell101") |>
  select(group, gene_exp) |>
  group_by(group) |>
  skimr::skim_without_charts() |>
  select(-c(n_missing,
            complete_rate,
            skim_variable,
            skim_type)) |>
  mutate(iqr = c(IQR(cell101_af), IQR(cell101_plac)))
cell101_data

cell101_data |> gt() |>
  fmt_number(
    columns = mean
  )



groups <- c("AF42", "placebo")
stdev_wild <- c(sd(wild_af), sd(wild_plac))
mean_wild <- c(mean(wild_af), mean(wild_plac))
iqr_wild <- c(IQR(wild_af), IQR(wild_plac))

wild_stats <- tibble(groups, stdev_wild, mean_wild, iqr_wild)
wild_stats

gt(wild_stats)


  gt() |>
  tab_header(
    title = "Summary Statistics",
    subtitle = "Wild Type Cell"
  )




