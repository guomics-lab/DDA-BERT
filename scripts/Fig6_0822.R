##Figure6_1
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(gghalves)


data <- read_excel("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_Figures.xlsx", sheet = 6)
colnames(data) <- c("file_name","FragPipe (MSBooster)","Sage","MS2Rescore","AlphaPept","AlphaPeptDeep","DDA-BERT")

df_hela <- data %>%
  filter(grepl("HeLa", file_name)) %>%
  mutate(
    sample_type = case_when(
      grepl("SPME_0_1ng", file_name) ~ "0.4 cell (HeLa cell)",
      grepl("0_5ng", file_name) ~ "2 cell (HeLa cell)",
      grepl("SPME_1ng", file_name) ~ "4 cell (HeLa cell)",
      grepl("10ng", file_name) ~ "40 cell (HeLa cell)",
      TRUE ~ NA_character_
    )
  )


df <- pivot_longer(df_hela, cols = starts_with(c("AlphaPept","AlphaPeptDeep","MS2Rescore","Sage","FragPipe (MSBooster)","DDA-BERT")), 
                   names_to = c("Tool"), 
                   names_pattern = "(.*)", 
                   values_to = "PSMs")

data_filtered <- df %>%
  filter(PSMs != 0)

data_summary_stats <- data_filtered %>%
  group_by(sample_type, Tool) %>%
  summarise(
    mean_PSMs = mean(PSMs),
    sd_PSMs = sd(PSMs)
  )%>%
  mutate(Tool = factor(Tool, levels = c("AlphaPept","AlphaPeptDeep","MS2Rescore","Sage","FragPipe (MSBooster)","DDA-BERT")))



tool_colors <- c( 
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "MS2Rescore"                = "#4575B4",
  "Sage"                   = "#762A83",
  "FragPipe (MSBooster)"   = "#1A9850",
  "DDA-BERT"               = "#D73027"
)


data_summary_stats <- data_summary_stats %>%
  mutate(Tool = factor(Tool, levels = names(tool_colors)))


ggplot(data_summary_stats, aes(x = sample_type, y = mean_PSMs, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_PSMs - sd_PSMs, ymax = mean_PSMs + sd_PSMs), 
    position = position_dodge(0.7), width = 0.2
  ) +
  labs(
    y = "# PSMs at 1% FDR",
    x = "",
    fill = "Tool"
  ) +
  theme_minimal() +
  scale_fill_manual(values = tool_colors) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 12000, by = 2000), 
    labels = c("0", "2000", "4000", "6000", "8000", "10000", "12000"), 
    limits = c(0, 12000),
    expand = c(0, 0)
  )

ggsave(filename = "Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig66_trace_sample_PSMs_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 4,
       height = 4, 
       units = "in") 
	   

##Figure6_2
df <- read_excel("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_Figures.xlsx", sheet = 1)

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
df_long$Tools <- factor(df_long$Tools, levels = c("AlphaPept","AlphaPeptDeep","MS2Rescore","Sage","FragPipe (MSBooster)","DDA-BERT"))


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
p_single_cell <- ggplot(df_long %>% filter(Group == "single_cell"), 
                        aes(x = fdr_cutoff, y = PSM_num, color = Tools, group = Tools)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
  labs(x = "FDR threshold", y = "# PSMs", title = "Single-cell Proteomics") +
  scale_color_manual(breaks = tool_breaks, values = tool_colors) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.15, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha('white', 0.8), color = 'black'),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    plot.margin = margin(3, 3, 3, 3),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1))


ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig6_single_cell_PSMs_barplot.pdf", width = 4, height = 4, units = "in")


##Figure6_3
library(ggplot2)
library(readxl)
library(tidyr)

data <- read_excel("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\DDA-BERT_Figures.xlsx", sheet = 6)

df_hela <- data %>%
  filter(grepl("HeLa_digest", file_name)) %>%
  mutate(
    sample_type = case_when(
      grepl("SPME_0_1ng", file_name) ~ "0.4 cell (HeLa cell)",
      grepl("0_5ng", file_name) ~ "2 cell (HeLa cell)",
      grepl("SPME_1ng", file_name) ~ "4 cell (HeLa cell)",
      grepl("10ng", file_name) ~ "40 cell (HeLa cell)",
      TRUE ~ NA_character_
    )
  )


df <- pivot_longer(df_hela, cols = starts_with(c("FragPipe (MSBooster)","Sage","MS2Rescore","AlphaPept","AlphaPeptDeep","DDA-BERT")), 
                   names_to = c("Tool"), 
                   names_pattern = "(.*)", 
                   values_to = "PSMs")

data_filtered <- df %>%
  filter(PSMs != 0)

data_summary_stats <- data_filtered %>%
  group_by(sample_type, Tool) %>%
  summarise(
    mean_PSMs = mean(PSMs),
    sd_PSMs = sd(PSMs)
  )%>%
  mutate(Tool = factor(Tool, levels = c("AlphaPept","AlphaPeptDeep","MS2Rescore","Sage","FragPipe (MSBooster)","DDA-BERT")))


tool_colors <- c( 
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "MS2Rescore"                = "#4575B4",
  "Sage"                   = "#762A83",
  "FragPipe (MSBooster)"   = "#1A9850",
  "DDA-BERT"               = "#D73027"
)


data_summary_stats <- data_summary_stats %>%
  mutate(Tool = factor(Tool, levels = names(tool_colors)))


ggplot(data_summary_stats, aes(x = sample_type, y = mean_PSMs, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_PSMs - sd_PSMs, ymax = mean_PSMs + sd_PSMs), 
    position = position_dodge(0.7), width = 0.2
  ) +
  labs(
    y = "# peptides at 1% FDR",
    x = "",
    fill = "Tool"
  ) +
  theme_minimal() +
  scale_fill_manual(values = tool_colors) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 6000, by = 1000), 
    labels = c("0", "1000", "2000", "3000", "4000", "5000",'6000'), 
    limits = c(0, 6000),
    expand = c(0, 0)
  )

ggsave(filename = "Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig6_trace_sample_Peptides_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 4,
       height = 4, 
       units = "in") 


##Figure6_4
library(ggplot2)
library(readxl)
library(tidyr)

rm(list=ls())
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
      grepl("^Substrate_comp_Rapid", file_name) ~ "single_cell",
      TRUE ~ "other"
    )
  )


df_long <- df_long %>%
  filter(dataset == "single_cell")

tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe (MSBooster)"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"                = "#4575B4"
)

tool_order <- c("AlphaPept","AlphaPeptDeep","MS2Rescore","Sage","FragPipe (MSBooster)","DDA-BERT")
df_long$Tool <- factor(df_long$Tool, levels = tool_order)


p <- ggplot(df_long, aes(x = Tool, y = peptide_num, fill = Tool)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = tool_colors) +
  labs(
    x = NULL,
    y = "# peptides at 1% FDR"
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

ggsave("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig6_sinble-cell_peptide_num.pdf", p,
       width = 4, height = 4, units = "in")

	   