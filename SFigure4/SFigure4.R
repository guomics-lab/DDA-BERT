##Sfigure_PSM_class_num
rm(list=ls())
library(tidyverse)
library(readxl)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(rlang)


##4A: cell-line
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "220627_GG", 
    "scan_num"
  )
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "220627_GG", 
    "scan_num"
  )
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "220627_GG", 
    "scan_num"
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
    breaks = seq(0, 2000, by = 400), 
    labels = c("0","400", "800", "1200", "1600","2000"), 
    limits = c(0, 2000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_cell-line_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##4B: tissue
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "20210715_Exploris1", 
    "scan_num"
  )
  
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "20210715_Exploris1", 
    "scan_num"
  )
  
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "20210715_Exploris1", 
    "scan_num"
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
    breaks = seq(0, 2500, by = 500), 
    labels = c("0","500", "1000", "1500", "2000","2500"), 
    limits = c(0, 2500),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_tissue_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 


##4C: single cell
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "Substrate_comp_Rapid", 
    "scan_num"
  )
  
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "Substrate_comp_Rapid", 
    "scan_num"
  )
  
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "Substrate_comp_Rapid", 
    "scan_num"
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
    breaks = seq(0, 1200, by = 300), 
    labels = c("0","300", "600", "900",'1200'), 
    limits = c(0, 1200),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_1cell_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##4D: timsTOF
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "N20220628wucl", 
    "scan_num"
  )
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "N20220628wucl", 
    "scan_num"
  )
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "N20220628wucl", 
    "scan_num"
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
    breaks = seq(0, 2500, by = 500), 
    labels = c("0","500", "1000", "1500", "2000","2500"), 
    limits = c(0, 2500),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_tims_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##4E: Astral
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "20230108_AST_Neo1", 
    "scan_num"
  )
  
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "20230108_AST_Neo1", 
    "scan_num"
  )
  
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "20230108_AST_Neo1", 
    "scan_num"
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
    breaks = seq(0, 2000, by = 400), 
    labels = c("0","400", "800", "1200", "1600","2000"), 
    limits = c(0, 2000),
    expand = c(0, 0)
  )

ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_Astral_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 



##4F: mouse
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
    "E:\\DDA-BERT\\psm_precursor_per_protein\\sage", 
    "031519-MouseLiver", 
    "scan_num"
  )
  
  combined_data_dda_bert <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\dda_bert", 
    "031519-MouseLiver", 
    "scan_num"
  )
  
  combined_data_fp <- read_and_assign_cells(
    "E:\\DDA-BERT\\psm_precursor_per_protein\\fp", 
    "031519-MouseLiver", 
    "scan_num"
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
    breaks = seq(0, 1200, by = 300), 
    labels = c("0","300", "600", "900",'1200'), 
    limits = c(0, 1200),
    expand = c(0, 0)
  )
ggsave(filename = "E:\\DDA-BERT\\psm_per_protein\\figures\\SFigure_mouse_PSMs_class.pdf", 
       plot = last_plot(),
       device = "pdf",
       width = 2.3,
       height = 2.3, 
       units = "in") 
