#SFigure4: human_protein_peptide_match
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)

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
    
    assign(paste0("tissue_", letters[i]), file_data, envir = .GlobalEnv)
    cell_data_list[[paste0("tissue_", letters[i])]] <- file_data   
    processed_data <- process_file_data(file_data, cleaned_sequence_col)
    processed_data$file <- paste0("tissue_", letters[i])
    
    combined_data[[i]] <- processed_data
  }
  
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}
#MS2Rescore
process_multiple_directories <- function() {
  dirs <- list(
    sage = "C:\\DDA-BERT_0105\\protein_peptide_match\\sage",
    msgpt = "C:\\DDA-BERT_0105\\protein_peptide_match\\DDA-BERT",
    fp = "C:\\DDA-BERT_0105\\protein_peptide_match\\FP",
    alphapept = "C:\\DDA-BERT_0105\\protein_peptide_match\\AlphaPept",
    pepdeep = "C:\\DDA-BERT_0105\\protein_peptide_match\\PeptDeep",
    ms2rescore = "C:\\DDA-BERT_0105\\protein_peptide_match\\MS2Rescore"
  )
  
  combined_data_list <- lapply(names(dirs), function(source) {
    data <- read_and_assign_cells(dirs[[source]], "20210715_Exploris1", "sequence_count")
    data$source <- source
    return(data)
  })
  
  final_combined_data <- bind_rows(combined_data_list)
  return(final_combined_data)
}

final_combined_data <- process_multiple_directories()

# Group & summarize
grouped_data <- final_combined_data %>%
  group_by(class, file, source) %>%
  summarise(count = n(), .groups = 'drop')

source_totals <- final_combined_data %>%
  distinct(protein, file, source) %>%
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

class_source_stats$source <- factor(class_source_stats$source, levels = c("alphapept", "pepdeep", "ms2rescore", "sage", "fp", "msgpt"))
class_source_stats$class <- factor(class_source_stats$class, levels = c("1", "2", "3-5", "6-10", ">10"))

# Plot
ggplot(class_source_stats, aes(x = class, y = mean_count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), 
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c(
    "alphapept" = "#7A1A24",
    "ms2rescore" = "#4575B4",
    "sage" = "#762A83",
    "pepdeep" = "#FDAE61",
    "fp" = "#1A9850",
    "msgpt" = "#D73027"
  )) +
  labs(x = "# peptides per protein group", y = "# protein groups") +
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
  ) +
  scale_y_continuous(
    breaks = seq(0, 1800, by = 300), 
    labels = c("0","300", "600", "900", "1200", "1300","1800"), 
    limits = c(0, 1800),
    expand = c(0, 0)
  )

ggsave(filename = "C:\\DDA-BERT_0105\\suppl_figure\\SFigure4_human_protein_peptide_match.pdf", 
       plot = last_plot(),
       device = "pdf",
       dpi = 300,
       width = 3,
       height = 3, 
       units = "in")

##SFigure4: fruit_fly protein_peptide_match
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)

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
    
    new_column_name <- paste0("fruit_fly_", letters[i], "_sequence")
    if ("sequence" %in% colnames(file_data) && cleaned_sequence_col != "sequence") {
      file_data <- file_data %>% rename(!!new_column_name := sequence)
    }
    
    assign(paste0("fruit_fly_", letters[i]), file_data, envir = .GlobalEnv)
    cell_data_list[[paste0("fruit_fly_", letters[i])]] <- file_data   
    processed_data <- process_file_data(file_data, cleaned_sequence_col)
    processed_data$file <- paste0("fruit_fly_", letters[i])
    
    combined_data[[i]] <- processed_data
  }
  
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}

process_multiple_directories <- function() {
  dirs <- list(
    sage = "C:\\DDA-BERT_0105\\protein_peptide_match\\sage",
    msgpt = "C:\\DDA-BERT_0105\\protein_peptide_match\\DDA-BERT",
    fp = "C:\\DDA-BERT_0105\\protein_peptide_match\\FP",
    alphapept = "C:\\DDA-BERT_0105\\protein_peptide_match\\AlphaPept",
    pepdeep = "C:\\DDA-BERT_0105\\protein_peptide_match\\PeptDeep",
    ms2rescore = "C:\\DDA-BERT_0105\\protein_peptide_match\\MS2Rescore"
  )
  
  combined_data_list <- lapply(names(dirs), function(source) {
    data <- read_and_assign_cells(dirs[[source]], "Ref6496", "sequence_count")
    data$source <- source
    return(data)
  })
  
  final_combined_data <- bind_rows(combined_data_list)
  return(final_combined_data)
}

final_combined_data <- process_multiple_directories()

# Group & summarize
grouped_data <- final_combined_data %>%
  group_by(class, file, source) %>%
  summarise(count = n(), .groups = 'drop')

source_totals <- final_combined_data %>%
  distinct(protein, file, source) %>%
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

class_source_stats$source <- factor(class_source_stats$source, levels = c("alphapept", "pepdeep", "ms2rescore", "sage", "fp", "msgpt"))
class_source_stats$class <- factor(class_source_stats$class, levels = c("1", "2", "3-5", "6-10", ">10"))

# Plot
ggplot(class_source_stats, aes(x = class, y = mean_count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), 
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c(
    "alphapept" = "#7A1A24",
    "pepdeep" = "#FDAE61",
    "sage" = "#762A83",
    "fp" = "#1A9850",
    "ms2rescore" = "#4575B4",
    "msgpt" = "#D73027"
  )) +
  labs(x = "# peptides per protein group", y = "# protein groups") +
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
  ) +
  scale_y_continuous(
    breaks = seq(0, 420, by = 70), 
    labels = c("0","70", "140", "210", "280","350","420"), 
    limits = c(0, 420),
    expand = c(0, 0)
  )

ggsave(filename = "C:\\DDA-BERT_0105\\suppl_figure\\SFigure4_fruit_fly_protein_peptide_match.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 3,
       height = 3, 
       units = "in")

##SFigure4: artha_protein_peptide_match
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)

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
    
    assign(paste0("arabidopsis_", letters[i]), file_data, envir = .GlobalEnv)
    cell_data_list[[paste0("arabidopsis_", letters[i])]] <- file_data   
    processed_data <- process_file_data(file_data, cleaned_sequence_col)
    processed_data$file <- paste0("arabidopsis_", letters[i])
    
    combined_data[[i]] <- processed_data
  }
  
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}

process_multiple_directories <- function() {
  dirs <- list(
    sage = "C:\\DDA-BERT_0105\\protein_peptide_match\\sage",
    msgpt = "C:\\DDA-BERT_0105\\protein_peptide_match\\DDA-BERT",
    fp = "C:\\DDA-BERT_0105\\protein_peptide_match\\FP",
    alphapept = "C:\\DDA-BERT_0105\\protein_peptide_match\\AlphaPept",
    pepdeep = "C:\\DDA-BERT_0105\\protein_peptide_match\\PeptDeep",
    ms2rescore = "C:\\DDA-BERT_0105\\protein_peptide_match\\MS2Rescore"
  )
  
  combined_data_list <- lapply(names(dirs), function(source) {
    data <- read_and_assign_cells(dirs[[source]], "Arabidopsis", "sequence_count")
    data$source <- source
    return(data)
  })
  
  final_combined_data <- bind_rows(combined_data_list)
  return(final_combined_data)
}

final_combined_data <- process_multiple_directories()

# Group & summarize
grouped_data <- final_combined_data %>%
  group_by(class, file, source) %>%
  summarise(count = n(), .groups = 'drop')

source_totals <- final_combined_data %>%
  distinct(protein, file, source) %>%
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

class_source_stats$source <- factor(class_source_stats$source, levels = c("alphapept", "pepdeep", "ms2rescore", "sage", "fp", "msgpt"))
class_source_stats$class <- factor(class_source_stats$class, levels = c("1", "2", "3-5", "6-10", ">10"))

# Plot
ggplot(class_source_stats, aes(x = class, y = mean_count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), 
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c(
    "alphapept" = "#7A1A24",
    "pepdeep" = "#FDAE61",
    "sage" = "#762A83",
    "fp" = "#1A9850",
    "ms2rescore" = "#4575B4",
    "msgpt" = "#D73027"
  )) +
  labs(x = "# peptides per protein group", y = "# protein groups") +
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
  ) +
  scale_y_continuous(
    breaks = seq(0, 2100, by = 700), 
    labels = c("0","700", "1400", "2100"), 
    limits = c(0, 2100),
    expand = c(0, 0)
  )

ggsave(filename = "C:\\DDA-BERT_0105\\suppl_figure\\SFigure4_Arabidopsis_protein_peptide_match.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 3,
       height = 3, 
       units = "in")


##SFigure4: yeast_protein_peptide_match
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)

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
    
    assign(paste0("yeast_", letters[i]), file_data, envir = .GlobalEnv)
    cell_data_list[[paste0("yeast_", letters[i])]] <- file_data   
    processed_data <- process_file_data(file_data, cleaned_sequence_col)
    processed_data$file <- paste0("yeast_", letters[i])
    
    combined_data[[i]] <- processed_data
  }
  
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}

process_multiple_directories <- function() {
  dirs <- list(
    sage = "C:\\DDA-BERT_0105\\protein_peptide_match\\sage",
    msgpt = "C:\\DDA-BERT_0105\\protein_peptide_match\\DDA-BERT",
    fp = "C:\\DDA-BERT_0105\\protein_peptide_match\\FP",
    alphapept = "C:\\DDA-BERT_0105\\protein_peptide_match\\AlphaPept",
    pepdeep = "C:\\DDA-BERT_0105\\protein_peptide_match\\PeptDeep",
    ms2rescore = "C:\\DDA-BERT_0105\\protein_peptide_match\\MS2Rescore"
  )
  
  combined_data_list <- lapply(names(dirs), function(source) {
    data <- read_and_assign_cells(dirs[[source]], "CY8611", "sequence_count")
    data$source <- source
    return(data)
  })
  
  final_combined_data <- bind_rows(combined_data_list)
  return(final_combined_data)
}

final_combined_data <- process_multiple_directories()

# Group & summarize
grouped_data <- final_combined_data %>%
  group_by(class, file, source) %>%
  summarise(count = n(), .groups = 'drop')

source_totals <- final_combined_data %>%
  distinct(protein, file, source) %>%
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

class_source_stats$source <- factor(class_source_stats$source, levels = c("alphapept", "pepdeep", "ms2rescore", "sage", "fp", "msgpt"))
class_source_stats$class <- factor(class_source_stats$class, levels = c("1", "2", "3-5", "6-10", ">10"))

# Plot
ggplot(class_source_stats, aes(x = class, y = mean_count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), 
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c(
    "alphapept" = "#7A1A24",
    "pepdeep" = "#FDAE61",
    "sage" = "#762A83",
    "fp" = "#1A9850",
    "ms2rescore" = "#4575B4",
    "msgpt" = "#D73027"
  )) +
  labs(x = "# peptides per protein group", y = "# protein groups") +
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
  ) +
  scale_y_continuous(
    breaks = seq(0, 1000, by = 200), 
    labels = c("0","200", "400", "600", "800", "1000"), 
    limits = c(0, 1000),
    expand = c(0, 0)
  )

ggsave(filename = "C:\\DDA-BERT_0105\\suppl_figure\\SFigure4_yeast_protein_peptide_match.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 3,
       height = 3, 
       units = "in")
