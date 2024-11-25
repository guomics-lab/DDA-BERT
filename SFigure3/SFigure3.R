##SFigure3: protein group number (delete peptide_num = 1)
SFigure3 <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "pro_num_delete_peptide1")
colnames(SFigure3)[5] <- "Sample_type"

##3A: cell-line
jurkat_data <- SFigure3 %>% filter(Sample_type == "Jurkat cell")

df <- jurkat_data %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )


data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 5000, by = 1000), 
    labels = c("0", "1000", "2000", "3000", "4000", "5000"), 
    limits = c(0, 5000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_cell-line_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##3B: tissue
tissue_data <- SFigure3 %>% filter(Sample_type == "CRC FFPE samples")

df <- tissue_data %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )



data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 7000, by = 1000), 
    labels = c("0","1000", "2000", "3000", "4000", "5000", "6000", "7000"), 
    limits = c(0, 7000),
    expand = c(0, 0)
  )


ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_tissue_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##3C: single cell
single_cell <- SFigure3 %>% filter(Sample_type == "1 cell (Hela cell)")

df <- single_cell %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )


data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 1200, by = 200), 
    labels = c("0","200", "400", "600", "800", "1000", "1200"), 
    limits = c(0, 1200),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_single_cell_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##3D: timsTOF
tims_data <- SFigure3 %>% filter(Sample_type == "HCC1806 and HS578T cell lines")

df <- tims_data %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )


data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 6000, by = 1000), 
    labels = c("0", "1000", "2000", "3000", "4000", "5000", "6000"), 
    limits = c(0, 6000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_tims_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##3E: Astral
Astral_data <- SFigure3 %>% filter(Sample_type == "HeLa_200ng")

df <- Astral_data %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )


data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 6000, by = 1000), 
    labels = c("0", "1000", "2000", "3000", "4000", "5000", "6000"), 
    limits = c(0, 6000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_Astral_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##3F: mouse
mouse_data <- SFigure3 %>% filter(Sample_type == "Mouse liver")

df <- mouse_data %>%
  pivot_longer(cols = c(FragPipe, Sage, 'DDA-BERT'), 
               names_to = "Tool", 
               values_to = "Proteins")

data_summary <- df %>%
  group_by(Tool) %>%
  summarise(
    mean_Proteins = mean(Proteins),
    sd_Proteins = sd(Proteins)
  )

data_summary$Tool <- factor(data_summary$Tool, levels = c("Sage", "FragPipe", "DDA-BERT"))

ggplot(data_summary, aes(x = Tool, y = mean_Proteins, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = mean_Proteins - sd_Proteins, ymax = mean_Proteins + sd_Proteins), width = 0.1, size = 0.4) +
  labs(y = "# proteins", x = "") +
  scale_fill_manual(values = c("DDA-BERT" = "#893b8f", "FragPipe" = "#1193cc", "Sage" = "#d97816")) +
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
    breaks = seq(0, 3500, by = 500), 
    labels = c("0","500", "1000", "1500", "2000", "2500", "3000", "3500"), 
    limits = c(0, 3500),
    expand = c(0, 0)
  )


ggsave(filename = "E:\\github\\DDA-BERT\\delete_peptide1_pros\\delete_peptide1_mouse_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 
