if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}


remotes::install_github("krassowski/complex-upset")


rm(list=ls())
library(UpSetR)

df <- read.csv("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\human_overlap_result0819.csv",row.names = 1)
colnames(df) <- c("FragPipe", "Sage", "AlphaPept","AlphaPeptDeep", "MS2Rescore","DDA-BERT")


tool_colors <- c(
  "AlphaPept"              = "#7A1A24",
  "AlphaPeptDeep"          = "#FDAE61",
  "DDA-BERT"               = "#D73027",
  "FragPipe"   = "#1A9850",
  "Sage"                   = "#762A83",
  "MS2Rescore"                = "#4575B4"
)

pdf("Z:\\members\\ajun\\1billion_MSGPT\\DDA-BERT_0809\\114_4_cls_epoch18\\main_figure\\fig5_human_peptide_overlap_upset_0722.pdf", width = 7, height = 7)

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
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1)
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