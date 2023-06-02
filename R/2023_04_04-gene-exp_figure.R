## Stats Consulting - Collaboration Project (Task 2)
## Douglas Dally - a1767927
## 04/04/2023

## Load Libs
install.packages("pacman")
pacman::p_load(tidyverse, here, patchwork, ggrepel, showtext, ggpubr, extrafont)

## Add font
font_add(
  family = "times",
  regular = here::here(
    "figs","times.ttf"
  )
)

font_import()
loadfonts(device="win")

## Load Data
data <- read_csv(
  here::here(
    "data", "2023-06-01_gene-exp_cleaned.csv"
  )
)
data


## Index positions of the points we want to label
idx_lab <- c(11, 22, 33, 43, 54, 65, 76, 87)

## Create new column exp_label with empty strings for all index not in idx_lab
data <- data |>
  mutate(exp_label = cell_type)

data$exp_label[-idx_lab] <- ""
data

## Make legend order AF42 then Placebo
data <- data |>
  mutate(treatment = as.factor(treatment),
         treatment = factor(treatment, levels = rev(levels(treatment))))
data

## Remove "GL-" from treatment names
data$exp_label <- str_remove_all(string = data$exp_label, pattern = "[GL-]")
data

## Step by step plots for wild type ##

# Step 1: get basic scatterplot
wild_p1 <- data |>
  filter(cell_line == "wild") |>
  ggplot(aes(x = concentration, y = gene_exp, label = exp_label)) +
  geom_point(aes(fill = treatment), size = 3, shape = 21, colour = "black") +
  theme_bw()
wild_p1

# Step 2: Add required x-axis and grid
wild_p2 <- wild_p1 +
  scale_x_continuous(breaks = seq(0,10,1), minor_breaks = seq(0,11,0.5))
wild_p2

# Step 3: Add required colour scheme and legend names
wild_p3 <- wild_p2 +
  scale_fill_manual(values = c("#83a7d4", "#d1bf94"),
                    labels = c("Activating Factor 42", "Placebo"))
wild_p3

# Step 4: Use ggrepel to label experiment name
wild_p4 <- wild_p3 +
  geom_label_repel(aes(fill = treatment),
                   col = "black",
                   xlim = c(-Inf, NA),
                   ylim = c(-Inf, NA),
                   min.segment.length = 0,
                   max.overlaps = Inf,
                   nudge_x = 1,
                   show.legend = FALSE,
                   size = 4,
                   family = "Times New Roman")
wild_p4

# Step 5: Add titles and themes
wild_p5 <- wild_p4 +
   labs(title = "Wild Type",
       tag = "A",
       fill = "Treatment",
       x = expression(paste(mu,"g/ml")),
       y = "Gene Expression") +
  theme(text = element_text(family = "Times New Roman"),
        legend.position = "bottom",
        legend.title=element_text(size=13),
        legend.text=element_text(size=12),
        plot.title = element_text(size=16),
        plot.tag = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title =element_text(size=13),
        aspect.ratio = 1)
wild_p5

## Same plot for cell-type 101
cell101_plot <- data |>
  filter(cell_line == "cell101") |>
  ggplot(aes(x = concentration, y = gene_exp, label = exp_label)) +
  geom_point(aes(fill = treatment), size = 3, shape = 21, colour = "black") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,10,1), minor_breaks = seq(0,11,0.5)) +
  scale_fill_manual(values = c("#83a7d4", "#d1bf94"),
                    labels = c("Activating Factor 42", "Placebo")) +
  geom_label_repel(aes(fill = treatment),
                   col = "black",
                   xlim = c(-Inf, NA),
                   ylim = c(-Inf, NA),
                   min.segment.length = 0,
                   max.overlaps = Inf,
                   nudge_x = 1,
                   show.legend = FALSE,
                   size = 4,
                   family = "Times New Roman") +
  labs(title = "Cell-Type 101",
       tag = "B",
       fill = "Treatment",
       x = expression(paste(mu,"g/ml")),
       y = "Gene Expression") +
  theme(text = element_text(family = "Times New Roman"),
        legend.position = "bottom",
        legend.title=element_text(size=13),
        legend.text=element_text(size=12),
        plot.title = element_text(size=16),
        plot.tag = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title =element_text(size=13),
        aspect.ratio = 1)
 cell101_plot

 ## Now combine two plots using ggarrange
 final_plot <- ggarrange(wild_p5, cell101_plot, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
 final_plot

 tiff(here::here("figs","2023-04-05_gene-exp-figure_Karl.tif"),
      width = 9,
      height = 6,
      units = 'in',
      res = 500)
 final_plot
 dev.off()


