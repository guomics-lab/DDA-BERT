##SFigure1_1
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)

rm(list=ls())
sage <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\main_figure\\Figure2_score\\Sage_Arabidopsis_merged_all_score.csv")
sage$Method <- "Sage"


dda <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\DDA-BERT_initial_Arabidopsis_merged_all_score.csv")
dda$Method <- "DDA-BERT"


data <- bind_rows(sage, dda)


data$LengthGroup <- factor(data$LengthGroup, levels = c("7-15", "16-24", "25-33", "34-42", "43-50"))
data$Type <- factor(data$Type, levels = c("Target", "Decoy"))
data$Method <- factor(data$Method, levels = c("Sage", "DDA-BERT"))


p <- ggplot(data, aes(x = LengthGroup, y = Score, fill = Type)) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Target"),
    side = "l", scale = "width", trim = FALSE
  ) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Decoy"),
    side = "r", scale = "width", trim = FALSE
  ) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(
    x = "Peptide length (amino acids)",
    y = "Score"
  ) +
  scale_fill_manual(values = c("Target" = "#35978f", "Decoy" = "#9970ab")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    axis.line = element_line()
  )


ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\sf1_Score_Arabidopsis.pdf", p, width = 7, height = 5)



##SFigure1_2
library(ggplot2)
library(dplyr)
library(gghalves)

rm(list=ls())
sage <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\main_figure\\Figure2_score\\Sage_yeast_merged_all_score.csv")
sage$Method <- "Sage"

dda <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\DDA-BERT_initial_yeast_merged_all_score.csv")
dda$Method <- "DDA-BERT"


data <- bind_rows(sage, dda)


data$LengthGroup <- factor(data$LengthGroup, levels = c("7-15", "16-24", "25-33", "34-42", "43-50"))
data$Type <- factor(data$Type, levels = c("Target", "Decoy"))
data$Method <- factor(data$Method, levels = c("Sage", "DDA-BERT"))


p <- ggplot(data, aes(x = LengthGroup, y = Score, fill = Type)) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Target"),
    side = "l", scale = "width", trim = FALSE
  ) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Decoy"),
    side = "r", scale = "width", trim = FALSE
  ) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") + 
  labs(
    x = "Peptide length (amino acids)",
    y = "Score"
  ) +
  scale_fill_manual(values = c("Target" = "#35978f", "Decoy" = "#9970ab")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    axis.line = element_line()
  )


ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\sf1_Score_yeast.pdf", p, width = 7, height = 5)



##SFigure1_3
library(ggplot2)
library(dplyr)
library(gghalves)

rm(list=ls())
sage <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\main_figure\\Figure2_score\\Sage_fruit_fly_merged_all_score.csv")
sage$Method <- "Sage"


dda <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\DDA-BERT_initial_fruit_fly_merged_all_score.csv")
dda$Method <- "DDA-BERT"


data <- bind_rows(sage, dda)


data$LengthGroup <- factor(data$LengthGroup, levels = c("7-15", "16-24", "25-33", "34-42", "43-50"))
data$Type <- factor(data$Type, levels = c("Target", "Decoy"))
data$Method <- factor(data$Method, levels = c("Sage", "DDA-BERT"))


p <- ggplot(data, aes(x = LengthGroup, y = Score, fill = Type)) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Target"),
    side = "l", scale = "width", trim = FALSE
  ) +
  gghalves::geom_half_violin(
    data = filter(data, Type == "Decoy"),
    side = "r", scale = "width", trim = FALSE
  ) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") + 
  labs(
    x = "Peptide length (amino acids)",
    y = "Score"
  ) +
  scale_fill_manual(values = c("Target" = "#35978f", "Decoy" = "#9970ab")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    axis.line = element_line()
  )

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\suppl_figure\\sf1_Score_fruit_fly.pdf", p, width = 7, height = 5)


