library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)

df <- read_excel("C:\\DDA-BERT_0105\\DDA-BERT_Figures_0105_ajun.xlsx", sheet = 3)

df_long <- df %>%
  pivot_longer(
    cols = c('FragPipe (MSBooster)', Sage, MS2Rescore, AlphaPept, AlphaPeptDeep, `DDA-BERT`),
    names_to = "Tool",
    values_to = "peptide_num"
  ) %>%
  mutate(
    dataset = case_when(
      grepl("WT", file_name) ~ "yeast",
      grepl("20210715_Exploris1", file_name) ~ "human",
      grepl("Arabidopsis_Cold", file_name) ~ "Arabidopsis",
      grepl("HeLa_digest", file_name) ~ "trace_sample",
      grepl("Ref6496_VK", file_name) ~ "fruit_fly",
      grepl("20230108_AST", file_name) ~ "astral",
      grepl("Substrate_comp_Rapid", file_name) ~ "single_cell",
      TRUE ~ "other"
    )
  ) %>%
  filter(!dataset %in% c("trace_sample", "astral", "single_cell", "other"))


df_long$dataset <- recode(df_long$dataset,
                          "yeast" = "Saccharomyces cerevisiae",
                          "human" = "Homo sapiens",
                          "Arabidopsis" = "Arabidopsis thaliana",
                          "fruit_fly" = "Drosophila melanogaster"
)


df_summary <- df_long %>%
  group_by(dataset, Tool) %>%
  summarise(
    mean_val = mean(peptide_num, na.rm = TRUE),
    sd_val = sd(peptide_num, na.rm = TRUE),
    .groups = "drop"
  )

tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe (MSBooster)"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"             = "#4575B4"
)

tool_order <- c("AlphaPept", "AlphaPeptDeep", "MS2Rescore", "Sage", "FragPipe (MSBooster)", "DDA-BERT")
df_summary$Tool <- factor(df_summary$Tool, levels = tool_order)

df_summary$dataset <- factor(df_summary$dataset, 
                             levels = c("Homo sapiens", "Drosophila melanogaster", 
                                        "Arabidopsis thaliana", "Saccharomyces cerevisiae"))


p <- ggplot(df_summary, aes(x = Tool, y = mean_val, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.2,
    position = position_dodge(0.6)
  ) +

  facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = tool_colors) +
  labs(
    x = NULL,
    y = "Number of precursors at 1% FDR"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "italic", size = 12),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )

ggsave("C:\\DDA-BERT_0105\\suppl_figure\\sf2_precursor_num.pdf", p,
       width = 7, height = 6, units = "in")
