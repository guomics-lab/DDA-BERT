library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)

##Figure5A: PSMs barplot
rm(list=ls())
df <- read_excel('E:\\github\\DDA-BERT\\DDA-BERT_figure_2024.xlsx', sheet = "mouse")
data_summary <- df %>%
  group_by(FDR, Tool) %>%
  summarise(
    mean_PSMs = mean(PSMs),
    sd_PSMs = sd(PSMs)
  )%>%
  mutate(Tool = factor(Tool, levels = c("Sage", "FragPipe", "DDA-BERT")))

ggplot(data_summary, aes(x = FDR, y = mean_PSMs / 10000, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = (mean_PSMs - sd_PSMs) / 10000, ymax = (mean_PSMs + sd_PSMs) / 10000), 
                width = 0.2, size = 0.4, position = position_dodge(0.6)) +
  labs(x = "FDR (%)", y = "# PSMs (×10^4)") +
  scale_fill_manual(values = c("#d97816", "#1193cc", "#893b8f"),
                    labels = c("Sage", "FragPipe", "DDA-BERT")) +
  theme_minimal() + 
  theme(
    legend.title = element_blank(),
    legend.position = c(0.8, 0.9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  scale_x_discrete(labels = c("0.1", "0.5", "1", "2", "3", "4", "5")) +
  scale_y_continuous(
    breaks = seq(0, 8, by = 2), 
    labels = c("0","2", "4", "6", "8"), 
    limits = c(0, 8),
    expand = c(0, 0)
  )  

ggsave(filename = "E:\\github\\DDA-BERT\\psm_barplot\\Figure5_mouse_PSMs.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 	


##Figure5B: precursor plot
mouse_data <- figure5 %>% filter(Sample_type == "Mouse liver")
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

ggplot(data_summary, aes(x = Tool, y = mean_Proteins / 10000, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = (mean_Proteins - sd_Proteins) / 10000, ymax = (mean_Proteins + sd_Proteins) / 10000), width = 0.1, size = 0.4) +
  labs(y = "# precursors (×10^4)", x = "") +
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
    breaks = seq(0, 4.2, by = 0.7), 
    labels = c("0","0.7", "1.4", "2.1", "2.8", "3.5", "4.2"), 
    limits = c(0, 4.2),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\precursor_barplot\\Figure5_mouse_precursor_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##Figure5C: peptide plot
mouse_data <- figure5 %>% filter(Sample_type == "Mouse liver")
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
  labs(y = "# peptides", x = "") +
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
    breaks = seq(0, 30000, by = 10000), 
    labels = c("0","10000", "20000", "30000"), 
    limits = c(0, 30000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\peptide_barplot\\Figure5_mouse_peptide_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##Figure5D: protein group barplot
mouse_data <- figure5 %>% filter(Sample_type == "Mouse liver")
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
    breaks = seq(0, 4000, by = 1000), 
    labels = c("0","1000", "2000", "3000", "4000"), 
    limits = c(0, 4000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\pro_barplot\\Figure5_mouse_pro_barplot.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##Figure5E: Venn diagram of protein groups
rm(list=ls())
read_and_merge_with_info <- function(directory, pattern) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE)
  data_list <- lapply(files, read.csv)
  combined_data <- bind_rows(data_list)
  unique_data <- distinct(combined_data)
  return(unique_data)
}

data_dda_bert <- read_and_merge_with_info("E:\\DDA-BERT\\peptide_per_protein\\dda_bert", "031519-MouseLiver")
data_sage <- read_and_merge_with_info("E:\\DDA-BERT\\peptide_per_protein\\sage", "031519-MouseLiver")
data_fp <- read_and_merge_with_info("E:\\DDA-BERT\\peptide_per_protein\\fp", "031519-MouseLiver")

list_dda_bert <- data_dda_bert$protein_group
list_sage <- data_sage$protein_group
list_fp <- data_fp$protein_group

only_dda_bert <- length(setdiff(list_dda_bert, union(list_sage, list_fp)))
only_sage <- length(setdiff(list_sage, union(list_dda_bert, list_fp)))
only_fp <- length(setdiff(list_fp, union(list_dda_bert, list_sage)))

dda_sage_overlap <- length(setdiff(intersect(list_dda_bert, list_sage), list_fp))
dda_fp_overlap <- length(setdiff(intersect(list_dda_bert, list_fp), list_sage))
sage_fp_overlap <- length(setdiff(intersect(list_sage, list_fp), list_dda_bert))

all_overlap <- length(intersect(list_dda_bert, intersect(list_sage, list_fp)))

fit <- euler(c(
  "DDA_BERT" = only_dda_bert, 
  "Sage" = only_sage, 
  "FragPipe" = only_fp,
  "DDA_BERT&Sage" = dda_sage_overlap,
  "DDA_BERT&FragPipe" = dda_fp_overlap,
  "Sage&FragPipe" = sage_fp_overlap,
  "DDA_BERT&Sage&FragPipe" = all_overlap
))

p <- plot(fit, 
          fills = c("#893b8f", "#d97816", "#1193cc"),
          edges = list(col = "black", lex = 1),
          labels = list(font = 1, cex = 1),
          quantities = list(font = 10, cex = 1)
)


ggsave("E:\\github\\DDA-BERT\\pro_venn\\Figure5_mouse_pro_venn.pdf", plot = p, device = "pdf", width = 1.7, height = 1.7)



##Figure5F: peptide_class_per_protein
rm(list=ls())
process_file_data <- function(file_data, cleaned_sequence_col) {
  if (!cleaned_sequence_col %in% colnames(file_data)) {
    stop(paste("Error: Specified 'cleaned_sequence' column not found in the file. Expected column:", cleaned_sequence_col))
  }
  
  categorized_data <- file_data %>%
    mutate(class = case_when(
      !!sym(cleaned_sequence_col) == 1 ~ "1",
      !!sym(cleaned_sequence_col) == 2 ~ "2",
      !!sym(cleaned_sequence_col) >= 3 & !!sym(cleaned_sequence_col) <= 5 ~ "3-5",
      !!sym(cleaned_sequence_col) >= 6 & !!sym(cleaned_sequence_col) <= 10 ~ "6-10",
      !!sym(cleaned_sequence_col) > 10 ~ ">10"
    ))
  
  return(categorized_data)
}

read_and_assign_cells <- function(directory, pattern, cleaned_sequence_col) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE)
  cell_data_list <- list()
  combined_data <- list()  
  for (i in seq_along(files)) {
    file_data <- read.csv(files[i], check.names = FALSE)
    
    if (!cleaned_sequence_col %in% colnames(file_data)) {
      print(paste("Warning: Specified 'cleaned_sequence' column not found in file:", files[i]))
      print(paste("Available columns in file:", paste(colnames(file_data), collapse = ", ")))
      next
    }
    
    new_column_name <- paste0("plasma_", letters[i], "_sequence")
    if ("sequence" %in% colnames(file_data) && cleaned_sequence_col != "sequence") {
      file_data <- file_data %>% rename(!!new_column_name := sequence)
    }
    
    assign(paste0("cell_", letters[i]), file_data, envir = .GlobalEnv)
    cell_data_list[[paste0("cell_", letters[i])]] <- file_data
    processed_data <- process_file_data(file_data, cleaned_sequence_col)
    processed_data$file <- paste0("cell_", letters[i])
    combined_data[[i]] <- processed_data
  }
  
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}


##plot
process_multiple_directories <- function() {
  combined_data_sage <- read_and_assign_cells(
    "E:\\DDA-BERT\\peptide_per_protein\\sage", 
    "031519-MouseLiver", 
    "sequence_count"
  )
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\peptide_per_protein\\dda_bert", 
    "031519-MouseLiver", 
    "sequence_count"
  )
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\peptide_per_protein\\fp", 
    "031519-MouseLiver", 
    "sequence_count"
  )
  
  final_combined_data <- bind_rows(
    combined_data_sage %>% mutate(source = "sage"),
    combined_data_dda_bert %>% mutate(source = "dda_bert"),
    combined_data_fp %>% mutate(source = "fp")
  )
  
  return(final_combined_data)
}

final_combined_data <- process_multiple_directories()

grouped_data <- final_combined_data %>%
  group_by(class, file, source) %>%
  summarise(count = n(), .groups = 'drop')

source_totals <- final_combined_data %>%
  distinct(protein_group, file, source) %>%
  group_by(file, source) %>%
  summarise(source_total_count = n(), .groups = 'drop')

grouped_data_with_counts <- grouped_data %>%
  left_join(source_totals, by = c("file", "source"))

class_source_stats <- grouped_data_with_counts %>%
  group_by(class, source) %>%
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    .groups = 'drop'
  )


class_source_stats$source <- factor(class_source_stats$source, levels = c("sage", "fp", "dda_bert"))
class_source_stats$class <- factor(class_source_stats$class, levels = c("1", "2", "3-5", "6-10", ">10"))


ggplot(class_source_stats, aes(x = class, y = mean_count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), 
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c("dda_bert" = "#893b8f", "fp" = "#1193cc", "sage" = "#d97816")) +
  labs(x = "# peptides per protein", y = "# proteins") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )  +
  scale_y_continuous(
    breaks = seq(0, 1200, by = 200), 
    labels = c("0","200", "400", '600',"800", "1000", "1200"), 
    limits = c(0, 1200),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\github\\DDA-BERT\\peptide_class_per_protein\\Figure5_mouse_peptide_class_per_protein.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##Figure5G: PSMs score distributions
rm(list=ls())
data <- read.csv("E:\\DDA-BERT\\score_distribution\\031519-MouseLiver.csv")
data$label <- ifelse(data$label == 1, "Target", "Decoy")
data$precursor_part <- sapply(strsplit(as.character(data$precursor_id), "_"), function(x) x[length(x)-2])

split_data <- split(data, data$label)
remove_duplicates <- function(df) {
  
  df <- df[order(-df$score), ]
  
  df <- df[!duplicated(df$precursor_part), ]
  return(df)
}

unique_data <- lapply(split_data, remove_duplicates)
final_data <- do.call(rbind, unique_data)

p <- ggplot(final_data, aes(x = score, y = sage_discriminant_score, color = label)) +
  geom_point(size = 0.3, shape = 16) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "orange")) +
  labs(x = "DDA-BERT score", y = "Sage score") +
  xlim(-0.2, 1.2) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha('white', 0)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

p_with_density <- ggMarginal(p, type = "density", fill = TRUE, groupFill = TRUE, size = 1.0, bw = 0.1, margins = "both", adjust = 2)

> table(final_data$label)
Decoy Target 
75682  88834 

ggsave(filename = "E:\\github\\DDA-BERT\\psm_score_distribution\\Figure5_mouse_PSMs_score.pdf", 
       plot = p_with_density,
       device = "pdf",
       width = 4,
       height = 4,
       units = "in")		   