##Compare with MaxQuant
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)

##S2A: tissue_PSM
Sfigure2 <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "Max_PSMs")
sf2_data <- Sfigure2 %>% filter(Samples == "CRC FFPE samples")

df <- sf2_data  %>%
  pivot_longer(cols = c('FragPipe (v21_1)', 'Sage (v0.14.6)', 'MaxQuant (v2.6.5.0)','DDA-BERT','AlphaPpet (v0.5.2)'), 
               names_to = "Tool", 
               values_to = "PSMs")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(PSMs),
    sd_Proteins = sd(PSMs)
  )


data_summary$Tool <- factor(data_summary$Tool, levels = c('MaxQuant (v2.6.5.0)','AlphaPpet (v0.5.2)', "Sage (v0.14.6)", "FragPipe (v21_1)", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe (v21_1)" = "#1193cc", "Sage (v0.14.6)" = "#d97816",'MaxQuant (v2.6.5.0)' = '#2a9d8f', 'AlphaPpet (v0.5.2)' = '#264653')) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100000, by = 20000), 
    labels = c("0", "20000", "40000", "60000", "80000", "100000"), 
    limits = c(0, 100000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\Compare_with_MaxQuant\\SFigure_tissue_MaxQuant_PSM_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##S2B: tissue_Protein groups
Sfigure2 <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "Max_pros")
sf2_data <- Sfigure2 %>% filter(Samples == "CRC FFPE samples")

df <- sf2_data  %>%
  pivot_longer(cols = c('FragPipe (v21_1)', 'Sage (v0.14.6)', 'MaxQuant (v2.6.5.0)','DDA-BERT','AlphaPpet (v0.5.2)'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )

data_summary$Tool <- factor(data_summary$Tool, levels = c('MaxQuant (v2.6.5.0)','AlphaPpet (v0.5.2)', "Sage (v0.14.6)", "FragPipe (v21_1)", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe (v21_1)" = "#1193cc", "Sage (v0.14.6)" = "#d97816",'MaxQuant (v2.6.5.0)' = '#2a9d8f', 'AlphaPpet (v0.5.2)' = '#264653')) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 10000, by = 2000), 
    labels = c("0", "2000", "4000", "6000", "8000", "10000"), 
    limits = c(0, 10000),
    expand = c(0, 0)
  )


ggsave(filename = "E:\\github\\DDA-BERT\\Compare_with_MaxQuant\\SFigure_MaxQuant_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in")



##S2C: single cell_PSM
Sfigure2 <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "Max_PSMs")
sf2_data <- Sfigure2 %>% filter(Samples == "1 cell (Hela cell)")

df <- sf2_data  %>%
  pivot_longer(cols = c('FragPipe (v21_1)', 'Sage (v0.14.6)', 'MaxQuant (v2.6.5.0)','DDA-BERT','AlphaPpet (v0.5.2)'), 
               names_to = "Tool", 
               values_to = "PSMs")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(PSMs),
    sd_Proteins = sd(PSMs)
  )

data_summary$Tool <- factor(data_summary$Tool, levels = c('MaxQuant (v2.6.5.0)','AlphaPpet (v0.5.2)', "Sage (v0.14.6)", "FragPipe (v21_1)", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, linewidth = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe (v21_1)" = "#1193cc", "Sage (v0.14.6)" = "#d97816",'MaxQuant (v2.6.5.0)' = '#2a9d8f', 'AlphaPpet (v0.5.2)' = '#264653')) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 10000, by = 2000), 
    labels = c("0", "2000", "4000", "6000", "8000", "10000"), 
    limits = c(0, 10000),
    expand = c(0, 0)
  )


ggsave(filename = "E:\\github\\DDA-BERT\\Compare_with_MaxQuant\\SFigure_1cell_MaxQuant_PSM_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##S2D: single cell_Protein groups
Sfigure2 <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "Max_pros")
sf2_data <- Sfigure2 %>% filter(Samples == "1 cell (Hela cell)")

df <- sf2_data  %>%
  pivot_longer(cols = c('FragPipe (v21_1)', 'Sage (v0.14.6)', 'MaxQuant (v2.6.5.0)','DDA-BERT','AlphaPpet (v0.5.2)'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )

data_summary$Tool <- factor(data_summary$Tool, levels = c('MaxQuant (v2.6.5.0)','AlphaPpet (v0.5.2)', "Sage (v0.14.6)", "FragPipe (v21_1)", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe (v21_1)" = "#1193cc", "Sage (v0.14.6)" = "#d97816",'MaxQuant (v2.6.5.0)' = '#2a9d8f', 'AlphaPpet (v0.5.2)' = '#264653')) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 2400, by = 400), 
    labels = c("0", "400", "800", "1200", "1600", "2000", '2400'), 
    limits = c(0, 2400),
    expand = c(0, 0)
  )


ggsave(filename = "E:\\github\\DDA-BERT\\Compare_with_MaxQuant\\SFigure_MaxQuant_1cell_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in")	   
