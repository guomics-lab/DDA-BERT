library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)


df <- read_excel("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_Figures.xlsx", sheet = 3) 

df_long <- df %>%
  pivot_longer(
    cols = c('FragPipe (MSBooster)', Sage, MS2Rescore, AlphaPept, AlphaPeptDeep, `DDA-BERT`), 
    names_to = "Tool",
    values_to = "peptide_num"
  )

df_long <- df_long %>%
  mutate(
    dataset = case_when(
      grepl("180min_", file_name) ~ "yeast",
      grepl("20210715_Exploris1", file_name) ~ "human",
      grepl("Arabidopsis_Cold", file_name) ~ "Arabidopsis",
      grepl("HeLa_digest", file_name) ~ "trace_sample",
      grepl("Ref6496_VK", file_name) ~ "fruit_fly",
      grepl("20230108_AST", file_name) ~ "astral",
      grepl("Substrate_comp_Rapid", file_name) ~ "single_cell",
      TRUE ~ "other"
    )
  )



df_long <- df_long %>%
  filter(!dataset %in% c("trace_sample", "astral", "single_cell"))

df_long$dataset <- recode(df_long$dataset,
                          "yeast" = "Saccharomyces cerevisiae",
                          "human" = "Homo sapiens",
                          "Arabidopsis" = "Arabidopsis thaliana",
                          "fruit_fly" = "Drosophila melanogaster"
)

tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "MS2Rescore"                = "#4575B4",
  "Sage"                   = "#762A83",
  "FragPipe (MSBooster)"   = "#1A9850",
  "DDA-BERT"               = "#D73027"
)

tool_order <- c("AlphaPept", "AlphaPeptDeep", "MS2Rescore", "Sage", "FragPipe (MSBooster)", "DDA-BERT")
df_long$Tool <- factor(df_long$Tool, levels = tool_order)

df_long$dataset <- factor(df_long$dataset, 
                          levels = c("Homo sapiens", "Drosophila melanogaster", 
                                     "Arabidopsis thaliana", "Saccharomyces cerevisiae"))


p <- ggplot(df_long, aes(x = Tool, y = peptide_num, fill = Tool)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = tool_colors) +
  labs(
    x = NULL,
    y = "Number of peptides at 1% FDR"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
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

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig4_peptide_num_20250722.pdf", p,
       width = 7, height = 6, units = "in")