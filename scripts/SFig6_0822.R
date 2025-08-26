##SFigure6: arabidopsis peptide upset
rm(list=ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(UpSetR)
df <- read.csv("D:\\paper\\guo-summer_25\\test_7_11\\test_8_18\\5-8\\data\\arabidopsis_overlap_result0818.csv",row.names = 1)
colnames(df) <- c("FragPipe", "Sage", "AlphaPept","AlphaPeptDeep", "MS2Rescore","DDA-BERT")

df <- df %>% mutate(across(everything(), as.numeric))

tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"                = "#4575B4"
)


pdf("D:\\paper\\guo-summer_25\\test_7_11\\test_8_18\\5-8\\png\\arabidopsis_peptide_overlap_upset.pdf", width = 7, height = 7)

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
  text.scale = 1.5,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size"
)

dev.off()

