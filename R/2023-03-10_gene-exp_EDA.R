## Stats Consulting - Collaboration Project
## Douglas Dally - a1767927
## 10/03/2023

# Load Libs
pacman::p_load(tidyverse, here, viridis, gt)

# Read in data
data <- read_csv(
  here::here(
    "data", "2023-06-01_gene-exp_cleaned.csv"
  )
)
data

## EDA

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

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp_hist.png"),
  plot = p1,
  width = 12,
  height = 9
)

# Histogram of Gene Expression for Wild and Type-101 Cell Lines
p2 <- ggplot(data, aes(x = gene_exp)) +
  geom_histogram(fill = "cyan4", colour = "black") +
  labs(title = "Histogram of Gene Expression for Wild and Type-101 Cell Lines") +
  xlab("Gene Expression") +
  ylab("Count") +
  theme_bw()+
  facet_wrap(vars(cell_line))
p2

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp-cell-line_hist.png"),
  plot = p2,
  width = 12,
  height = 9
)

### Scatterplots
# Scatterplot of Gene Expression vs. Concentration
p3 <- ggplot(data, aes(x = concentration,
                       y = gene_exp,
                       col = treatment)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = treatment), method = "lm") +
  labs(title = "Scatterplot of Gene Expression vs. Concentration",
       subtitle = "Coloured by Treatment (Activating Factor 42 or Placebo)",
       col = "Treatment Group",
       x = expression(paste("Growth Factor Concentration ", mu,"g/ml")),
       y = "Gene Expression") +
  theme_bw()
p3

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp_scatter.png"),
  plot = p3,
  width = 12,
  height = 9
)

# Scatterplot of Gene Expression vs Concentration (Both Types)
p4 <- ggplot(
  data = data,
  aes(x = concentration, y = gene_exp)) +
  geom_point(aes(col = treatment), size = 3) +
  geom_smooth(aes(col = treatment), method = "lm") +
  facet_wrap(vars(cell_line)) +
  labs(title = "Scatterplot of Gene Expression vs. Concentration ",
       subtitle = "Coloured by Treatment Group (Placebo or Activating Factor 42)",
       col = "Treatment Group",
       x = expression(paste("Growth Factor Concentration ", mu,"g/ml")),
       y = "Gene Expression") +
  theme_bw()
p4

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp-cell-line_scatter.png"),
  plot = p4,
  width = 12,
  height = 9
)

p5 <- ggplot(
  data = data,
  aes(x = concentration, y = gene_exp)) +
  geom_point(aes(col = cell_type, shape = treatment), size = 2) +
  geom_smooth(aes(col = cell_type), method = "lm") +
  facet_wrap(vars(cell_line)) +
  labs(title = "Scatterplot of Gene Expression vs. Concentration ",
       subtitle = "Coloured by Cell Type",
       col = "Cell Type",
       shape = "Treatment Type",
       x = expression(paste("Growth Factor Concentration (", mu,"g/ml)")),
       y = "Gene Expression") +
  theme_bw()
p5

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp-cell-type_scatter.png"),
  plot = p5,
  width = 12,
  height = 9
)

### Boxplots
# Boxplot of Gene Expression by Treatment Group
p6 <- ggplot(
  data = data,
  aes(x = treatment, y = gene_exp, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(vars(cell_line)) +
  labs(title = "Boxplot of Gene Expression by Treatment Group (Placebo or Activating Factor 42)",
       fill = "Treatment Group") +
  ylab("Gene Expression") +
  xlab("Treatment Group") +
  theme_bw()
p6

ggsave(
  filename = here::here("figs", "2023-06-01_gene-exp-cell-line_box.png"),
  plot = p6,
  width = 12,
  height = 9
)


# TABLES ----
options(digits = 3)
# Wild Table
wild_af<- data |>
  filter(cell_line == "wild", treatment == "AF42") |>
  select(gene_exp) |>
  pull()

wild_plac <- data |>
  filter(cell_line == "wild", treatment == "placebo") |>
  select(gene_exp) |>
  pull()

wild_data <- data |>
  filter(cell_line == "wild") |>
  select(treatment, gene_exp) |>
  group_by(treatment) |>
  skimr::skim_without_charts() |>
  select(-c(n_missing,
            complete_rate,
            skim_variable,
            skim_type)) |>
  mutate(iqr = c(IQR(wild_af), IQR(wild_plac)))
wild_data

# Cell101 Table
cell101_af<- data |>
  filter(cell_line == "cell101", treatment == "AF42") |>
  select(gene_exp) |>
  pull()

cell101_plac <- data |>
  filter(cell_line == "cell101", treatment == "placebo") |>
  select(gene_exp) |>
  pull()

cell101_data <- data |>
  filter(cell_line == "cell101") |>
  select(treatment, gene_exp) |>
  group_by(treatment) |>
  skimr::skim_without_charts() |>
  select(-c(n_missing,
            complete_rate,
            skim_variable,
            skim_type)) |>
  mutate(iqr = c(IQR(cell101_af), IQR(cell101_plac)))
cell101_data

# Summary Tables
wild_summary <- data.frame(
  treatment = c("placebo", "AF42"),
  min = wild_data$numeric.p0,
  max = wild_data$numeric.p100,
  median = wild_data$numeric.p50,
  mean = wild_data$numeric.mean,
  std_dev = wild_data$numeric.sd,
  iqr = wild_data$iqr)

wild_summary |> gt() |>
  tab_header(
    title = "Summary Statistics",
    subtitle = "Wild Type Cell Line"
  ) |> gtsave(
    here::here("tabs", "2023-06-01_wild-summary.docx")
  )

cell101_summary <- data.frame(
  treatment = c("placebo", "AF42"),
  min = cell101_data$numeric.p0,
  max = cell101_data$numeric.p100,
  median = cell101_data$numeric.p50,
  mean = cell101_data$numeric.mean,
  std_dev = cell101_data$numeric.sd,
  iqr = cell101_data$iqr)

cell101_summary |> gt() |>
  tab_header(
    title = "Summary Statistics",
    subtitle = "Type-101 Cell Line"
  ) |> gtsave(
    here::here("tabs", "2023-06-01_cell101-summary.docx")
  )





