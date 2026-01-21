##Figure 2C
library(ggplot2)
library(dplyr)

rm(list=ls())
dda <- read.csv("C:\\DDA-BERT\\DDA-BERT_0105\\score_distribution\\dda-bert\\dda-bert_human_merged.csv")

dda$Type <- factor(dda$Type, levels = c("Target", "Decoy"))

p <- ggplot(dda, aes(x = factor(sequence_len), y = Score, fill = Type)) +
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

ggsave("C:\\DDA-BERT\\DDA-BERT_0105\\main_figure\\fig2_Score_human_pep_length.pdf", p, width = 8, height = 8)

