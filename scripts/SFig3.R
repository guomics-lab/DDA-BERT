library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)

rm(list=ls())
data <-  read.csv("C:\\DDA-BERT_0105\\DBSE_all_proteins_stats_20260111.csv")

colnames(data) <- c("filename", "FragPipe", "DDA-BERT", "Sage", "AlphaPept","AlphaPeptDeep", "MS2Rescore")

df <- pivot_longer(data, cols = starts_with(c("FragPipe","Sage","MS2Rescore","AlphaPept","AlphaPeptDeep","DDA-BERT")), 
                   names_to = c("Tool"), 
                   values_to = "Proteins")

df_long <- df %>%
  mutate(
    
    dataset = factor(
      case_when(
        grepl("WT", filename) ~ "yeast",
        grepl("20210715_Exploris1", filename) ~ "human",
        grepl("Arabidopsis_Cold", filename) ~ "Arabidopsis",
        grepl("HeLa_digest", filename) ~ "trace_sample",
        grepl("Ref6496_VK", filename) ~ "fruit_fly",
        grepl("20230108_AST", filename) ~ "astral",
        grepl("Substrate_comp_Rapid", filename) ~ "single_cell",
        TRUE ~ "other"
      ),
      levels = c("human", "fruit_fly", "Arabidopsis", "yeast")  
    )
  )

df_long <- df_long %>%
  filter(! is.na(dataset))

data_summary_stats <- df_long %>%
  group_by(dataset, Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins, na.rm = TRUE),
    sd_Proteins = sd(Proteins, na.rm = TRUE),
    .groups = "drop"  
  ) %>%

  mutate(
    Tool = factor(
      Tool,
      levels = c("AlphaPept", "AlphaPeptDeep", "MS2Rescore", "Sage", "FragPipe", "DDA-BERT")
    )
  )


tool_colors <- c( 
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "MS2Rescore"             = "#4575B4",
  "Sage"                   = "#762A83",
  "FragPipe"               = "#1A9850",
  "DDA-BERT"               = "#D73027"
)


data_summary_stats <- data_summary_stats %>%
  filter(dataset %in% c("human", "fruit_fly", "Arabidopsis", "yeast"))

ggplot(data_summary_stats, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins),
    position = position_dodge(0.7), width = 0.2
  ) +
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") + 
  labs(
    y = "# protein groups",
    x = NULL,
    fill = "Tool"
  ) +
  scale_fill_manual(values = tool_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "italic", size = 12),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),#, linewidth = 0.6),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

ggsave(filename = "C:\\DDA-BERT_0105\\suppl_figure\\SFigure_3_proteins_all_20260111.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 8,
       height = 8, 
       dpi = 300,
       units = "in") 