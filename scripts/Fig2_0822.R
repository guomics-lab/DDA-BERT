##Figure2_1

library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)


sage <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\main_figure\\Figure2_score\\Sage_human_merged_all_score.csv")
sage$Method <- "Sage"


dda <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_initial_human_merged_all_score.csv")
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

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig2_Score_left.pdf", p, width = 7, height = 5)


##Figure2_2
library(ggplot2)
library(dplyr)

dda <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_initial_human_merged_all_score_length.csv")


dda$Type <- factor(dda$Type, levels = c("Target", "Decoy"))
dda$Length <- as.integer(dda$Length)



p <- ggplot(dda, aes(x = factor(Length), y = Score, fill = Type)) +
  geom_boxplot(outlier.size = 0.1, width = 0.6) +
  facet_wrap(~ Type, ncol = 1, scales = "free_y") +
  labs(
    x = "Peptide length (amino acids)",
    y = "Score"
  ) +
  scale_fill_manual(values = c("Target" = "#35978f", "Decoy" = "#9970ab")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    axis.line = element_line(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig2_Score_human_pep_length.pdf", p, width = 8, height = 8)

