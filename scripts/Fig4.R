library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)

##1. Figures 4A-D
df <- read_excel("C:\\DDA-BERT\\DDA-BERT_0105\\DDA-BERT_Figures_0105.xlsx", sheet = 4)

df_long <- df %>%
  pivot_longer(
    cols = c('FragPipe (MSBooster)', Sage, MS2Rescore, AlphaPept, AlphaPeptDeep, `DDA-BERT`), 
    names_to = "Tool",
    values_to = "peptide_num"
  )

df_long <- df_long %>%
  mutate(
    dataset = case_when(
      grepl("WT", file_name) ~ "yeast",
      grepl("20210715_Exploris1", file_name) ~ "human",
      grepl("Arabidopsis_Cold", file_name) ~ "Arabidopsis",
      grepl("HeLa_digest", file_name) ~ "trace_sample",
      grepl("Ref6496_VK", file_name) ~ "fruit_fly",
      grepl("20230108_AST", file_name) ~ "astral",
      grepl("20141208", file_name) ~ "single_cell",
      #TRUE ~ "other"
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

ggsave("C:\\DDA-BERT\\DDA-BERT_0105\\main_figure\\fig4_peptide_num_20260105.pdf", p,
       width = 7, height = 6, units = "in")



##2. Figures 4E

rm(list=ls())
library(UpSetR)

df <- read.csv("C:\\DDA-BERT\\DDA-BERT_0105\\human_overlap_20260105.csv",row.names = 1)
colnames(df) <- c("FragPipe", "Sage", "AlphaPept","AlphaPeptDeep", "MS2Rescore","DDA-BERT")


tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"                = "#4575B4"
)

pdf("C:\\DDA-BERT\\DDA-BERT_0105\\main_figure\\fig4_human_peptide_overlap_upset_0105.pdf", width = 7, height = 7)

set_sizes <- colSums(df)
ordered_sets <- names(sort(set_sizes, decreasing = TRUE))

ordered_colors <- tool_colors[ordered_sets]

upset(
  df,
  sets = ordered_sets,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = ordered_colors,
  main.bar.color = "#2c3e50",
  matrix.color = "black",
  matrix.dot.alpha = 1,
  point.size = 1.5,
  line.size = 0.7,
  number.angles = 0,
  show.numbers = "yes",
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1),
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size"
)

dev.off()


exclusive_counts <- sapply(
  colnames(df),
  function(tool) {
    sum(df[[tool]] == 1 & rowSums(df[, setdiff(colnames(df), tool)]) == 0)
  }
)

print("Exclusive peptides (identified only by each tool):")
print(exclusive_counts)
	   