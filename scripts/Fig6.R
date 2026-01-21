##Figure 6A
library(ggplot2)
library(dplyr)
library(readr)

df_plot_csv <- "/ajun/DDA_BERT_manu/df_plot.csv"
df <- read_csv(df_plot_csv, show_col_types = FALSE)

score_col <- "score"
order <- c("train", "holdout")

if (!("label_group" %in% colnames(df)) && ("label" %in% colnames(df))) {
  df <- df %>% mutate(label_group = dplyr::recode(label, `1` = "target", `0` = "decoy"))
}

df <- df %>%
  filter(pep_group %in% order) %>%
  mutate(pep_group = factor(pep_group, levels = order, ordered = TRUE))


out_pdf <- "/ajun/DDA_BERT_manu/boxplot.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

palette <- c(train = "#e07067", holdout = "#25acae")
pep_order <- order
plot_df <- df

pdf(out_pdf, width = 6.5, height = 4.2)

for (lg in c("decoy", "target")) {
  sub <- plot_df %>%
    filter(label_group == lg) %>%
    filter(!is.na(pep_group), !is.na(score)) %>%
    mutate(pep_group = factor(pep_group, levels = pep_order, ordered = TRUE))

  if (nrow(sub) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(lg) +
      annotate("text", x = 0, y = 0, label = "No data")
  } else {
    p <- ggplot(sub, aes(x = pep_group, y = .data[["score"]], fill = pep_group)) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.5,
        color = "black",
        linewidth = 1.2
      ) +
      scale_fill_manual(values = palette, drop = FALSE) +
      labs(title = paste0("Score Boxplot \u2014 ", lg), x = NULL, y = "score") +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
      )
  }

  print(p)

  cat("\n[", lg, "] n per pep_group:\n", sep = "")
  counts <- sub %>%
    count(pep_group, name = "n") %>%
    right_join(
      data.frame(pep_group = factor(pep_order, levels = pep_order, ordered = TRUE)),
      by = "pep_group"
    ) %>%
    mutate(n = ifelse(is.na(n), 0L, as.integer(n)))

  for (i in seq_len(nrow(counts))) {
    cat(as.character(counts$pep_group[i]), " ", counts$n[i], "\n", sep = "")
  }
}

dev.off()

cat("\nSaved: ", out_pdf, "\n", sep = "")



##Figure 6B
library(ggplot2)

# Same logic: plot two KDEs (label==0 and label==1) filled, with legend, then save as PDF.
p <- ggplot() +
  geom_density(
    data = subset(df_all, label == 0),
    aes(x = score, fill = "Decoy", color = "Decoy"),
    alpha = 0.4
  ) +
  geom_density(
    data = subset(df_all, label == 1),
    aes(x = score, fill = "Target", color = "Target"),
    alpha = 0.4
  ) +
  labs(x = "score", y = "Density") +
  guides(fill = guide_legend(title = "Peptide Type"),
         color = guide_legend(title = "Peptide Type")) +
  theme_classic() +
  theme(
    figure = element_rect(size = 2),
    legend.position = "right"
  )

out_path <- "/ajun/DDA_BERT_manu/human_intensity.pdf"

ggsave(out_path, plot = p, device = "pdf", width = 2, height = 2, units = "in")



##Figure 6C
##file: 220627_GG_02_rep3_dil10
source("D:\\DDA-BERT\\FDP\\run_fdp_calc.R")
color_mapping <- c("Paired method" = "#7CAE00", "Sample method" = "#C77CFF", "Lower bound" = "#00BFC4", "Combined method" = "#F8766D")

report_file <- "D:\\DDA-BERT\\FDP\\data\\220627_GG_02_rep3_dil10_peptide_100pct_processed.tsv"
pep_file <- "D:\\DDA-BERT\\FDP\\fasta\\human_entrapment_diann_pep.txt"

pro_fdp_file1 <- run_dda_bert_fdp_analysis(report_file, level = "peptide", prefix = "220627_GG_02_rep3_dil10", pep_file = pep_file, k_fold = 1)

fdp_data <- read.csv(pro_fdp_file1)

if("combined_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "combined_fdp"] <- "FDP"
} else {
  warning("Column 'combined_fdp' not found; cannot rename")
}

if("paired_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "paired_fdp"] <- "FDP_1B"
} else {
  warning("Column 'paired_fdp' not found; cannot rename")
}

temp_fdp_file <- tempfile(fileext = ".csv")
write.csv(fdp_data, temp_fdp_file, row.names = FALSE)

gg3 <- plot_fdp_fdr_v2(temp_fdp_file, fdr_max = 0.10, fig_title = "220627_GG_02_rep3_dil10", add_numbers = TRUE, r = 1, color_mapping = color_mapping, fdr_decimal_place = 2)

ggsave("D:\\DDA-BERT\\FDP\\220627_GG_02_rep3_dil10.pdf", gg3, width = 4, height = 4, dpi = 300)



##file: 220627_GG_04_rep1_dil10
source("D:\\DDA-BERT\\FDP\\run_fdp_calc.R")
color_mapping <- c("Paired method" = "#7CAE00", "Sample method" = "#C77CFF", "Lower bound" = "#00BFC4", "Combined method" = "#F8766D")

report_file <- "D:\\DDA-BERT\\FDP\\data\\220627_GG_04_rep1_dil10_peptide_100pct_processed.tsv"
pep_file <- "D:\\DDA-BERT\\FDP\\fasta\\human_entrapment_diann_pep.txt"

pro_fdp_file1 <- run_dda_bert_fdp_analysis(report_file, level = "peptide", prefix = "220627_GG_04_rep1_dil10", pep_file = pep_file, k_fold = 1)

fdp_data <- read.csv(pro_fdp_file1)

if("combined_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "combined_fdp"] <- "FDP"
} else {
  warning("Column 'combined_fdp' not found; cannot rename")
}

if("paired_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "paired_fdp"] <- "FDP_1B"
} else {
  warning("Column 'paired_fdp' not found; cannot rename")
}

temp_fdp_file <- tempfile(fileext = ".csv")
write.csv(fdp_data, temp_fdp_file, row.names = FALSE)

gg3 <- plot_fdp_fdr_v2(temp_fdp_file, fdr_max = 0.10, fig_title = "220627_GG_04_rep1_dil10", add_numbers = TRUE, r = 1, color_mapping = color_mapping, fdr_decimal_place = 2)

ggsave("D:\\DDA-BERT\\FDP\\220627_GG_04_rep1_dil10.pdf", gg3, width = 4, height = 4, dpi = 300)



##file: 220627_GG_06_rep1_dil10
source("D:\\DDA-BERT\\FDP\\run_fdp_calc.R")
color_mapping <- c("Paired method" = "#7CAE00", "Sample method" = "#C77CFF", "Lower bound" = "#00BFC4", "Combined method" = "#F8766D")

report_file <- "D:\\DDA-BERT\\FDP\\data\\220627_GG_06_rep1_dil10_peptide_100pct_processed.tsv"
pep_file <- "D:\\DDA-BERT\\FDP\\fasta\\human_entrapment_diann_pep.txt"

pro_fdp_file1 <- run_dda_bert_fdp_analysis(report_file, level = "peptide", prefix = "220627_GG_06_rep1_dil10", pep_file = pep_file, k_fold = 1)

fdp_data <- read.csv(pro_fdp_file1)

if("combined_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "combined_fdp"] <- "FDP"
} else {
  warning("Column 'combined_fdp' not found; cannot rename")
}

if("paired_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "paired_fdp"] <- "FDP_1B"
} else {
  warning("Column 'paired_fdp' not found; cannot rename")
}

temp_fdp_file <- tempfile(fileext = ".csv")
write.csv(fdp_data, temp_fdp_file, row.names = FALSE)

gg3 <- plot_fdp_fdr_v2(temp_fdp_file, fdr_max = 0.10, fig_title = "220627_GG_06_rep1_dil10", add_numbers = TRUE, r = 1, color_mapping = color_mapping, fdr_decimal_place = 2)

ggsave("D:\\DDA-BERT\\FDP\\220627_GG_06_rep1_dil10.pdf", gg3, width = 4, height = 4, dpi = 300)



##file: 220627_GG_07_rep2_dil10
source("D:\\DDA-BERT\\FDP\\run_fdp_calc.R")
color_mapping <- c("Paired method" = "#7CAE00", "Sample method" = "#C77CFF", "Lower bound" = "#00BFC4", "Combined method" = "#F8766D")

report_file <- "D:\\DDA-BERT\\FDP\\data\\220627_GG_07_rep2_dil10_peptide_100pct_processed.tsv"
pep_file <- "D:\\DDA-BERT\\FDP\\fasta\\human_entrapment_diann_pep.txt"

pro_fdp_file1 <- run_dda_bert_fdp_analysis(report_file, level = "peptide", prefix = "220627_GG_07_rep2_dil10", pep_file = pep_file, k_fold = 1)

fdp_data <- read.csv(pro_fdp_file1)

if("combined_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "combined_fdp"] <- "FDP"
} else {
  warning("Column 'combined_fdp' not found; cannot rename")
}

if("paired_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "paired_fdp"] <- "FDP_1B"
} else {
  warning("Column 'paired_fdp' not found; cannot rename")
}

temp_fdp_file <- tempfile(fileext = ".csv")
write.csv(fdp_data, temp_fdp_file, row.names = FALSE)

gg3 <- plot_fdp_fdr_v2(temp_fdp_file, fdr_max = 0.10, fig_title = "220627_GG_07_rep2_dil10", add_numbers = TRUE, r = 1, color_mapping = color_mapping, fdr_decimal_place = 2)

ggsave("D:\\DDA-BERT\\FDP\\220627_GG_07_rep2_dil10.pdf", gg3, width = 4, height = 4, dpi = 300)