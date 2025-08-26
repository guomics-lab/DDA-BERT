library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)

df <- read_excel("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_Figures_0820.xlsx", sheet = 1)

tool_labels <- c(
  rep("FragPipe (MSBooster)", 5),
  rep("Sage", 5),
  rep("MS2Rescore", 5),
  rep("AlphaPept", 5),
  rep("AlphaPeptDeep", 5),
  rep("DDA-BERT", 5)
)

colnames(df)[-1] <- paste0(tool_labels, "_", 1:5)
df_long <- df %>%
  pivot_longer(
    cols = -filename,
    names_to = "method_rep",
    values_to = "count"
  ) %>%
  mutate(
    method = sub("_.*", "", method_rep)
  )

df_long <- df_long %>%
  mutate(
    method_rep = case_when(
      grepl("_1$", method_rep) ~ "0.002",
      grepl("_2$", method_rep) ~ "0.004",
      grepl("_3$", method_rep) ~ "0.006",
      grepl("_4$", method_rep) ~ "0.008",
      grepl("_5$", method_rep) ~ "0.01",
      TRUE ~ method_rep
    )
  )

df_long <- df_long %>%
  mutate(
    dataset = case_when(
      grepl("180min_", filename) ~ "yeast",
      grepl("20210715_Exploris1", filename) ~ "human",
      grepl("Arabidopsis_Cold", filename) ~ "Arabidopsis",
      grepl("HeLa_digest", filename) ~ "trace_sample",
      grepl("Ref6496_VK", filename) ~ "fruit_fly",
      grepl("20230108_AST", filename) ~ "astral",
      grepl("Substrate_comp_Rapid", filename) ~ "single_cell",
      TRUE ~ "other"
    )
  )

colnames(df_long) <- c("file_name", "fdr_cutoff", "PSM_num", "Tools", "Group")



df_long$fdr_cutoff <- as.numeric(df_long$fdr_cutoff)

df_long$Tools <- factor(df_long$Tools, levels = c(
  "AlphaPept",
  "AlphaPeptDeep",
  "MS2Rescore",
  "Sage",
  "FragPipe (MSBooster)",
  "DDA-BERT"
))


library(ggplot2)
library(patchwork)

df_long$Group <- as.character(df_long$Group)

tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe (MSBooster)"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"                = "#4575B4"
)

tool_breaks <- levels(df_long$Tools)

plot_group <- function(group_name, title_text) {
  ggplot(df_long %>% filter(Group == group_name), 
         aes(x = fdr_cutoff, y = PSM_num, color = Tools, group = Tools)) +
    stat_summary(fun = mean, geom = "line", size = 0.6) +
    labs(x = "FDR threshold", y = "# PSMs") +
    ggtitle(title_text) +
    scale_color_manual(breaks = tool_breaks, values = tool_colors) +

    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),
      axis.title = element_text(family = "Arial", size = 10),
      axis.text = element_text(family = "Arial", size = 9),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(0.15, "cm"), 
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.background = element_rect(fill = alpha('white', 0.8), color = 'black'),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
      plot.margin = margin(3, 3, 3, 3)
    )
}

p4 <- plot_group("yeast", "Saccharomyces cerevisiae")
p1 <- plot_group("human", "Homo sapiens")
p3 <- plot_group("Arabidopsis", "Arabidopsis thaliana")
p2 <- plot_group("fruit_fly", "Drosophila melanogaster")

combined_plot <- ((p1 | p2) / (p3 | p4)) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") & 
  guides(color = guide_legend(nrow = 1))

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig3_PSM_FDR_cutoff_20250722.pdf", combined_plot,
       width = 8, height = 7.5, units = "in", device = cairo_pdf)
