options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

library(readr)
library(dplyr)
library(ggplot2)
library(DBI)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)

run_dda_bert_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    if(level=="peptide"){
        b <- a %>% select(cleaned_sequence, modified_sequence, score, q_value) %>% 
            arrange(desc(score)) %>% 
            group_by(cleaned_sequence) %>%  
            filter(row_number() == 1) %>%
            rename(peptide = cleaned_sequence, mod_peptide = modified_sequence) %>%             
            ungroup()                                       
        set.seed(2024)

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")
        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b, out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-std_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar D:\\paper\\guo-summer_25\\test_7_11\\fdrbench-0.0.1\\fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:1 -miss_c 2 -enzyme 1 -maxLength 50",sep="")
            }else{
                cmd <- paste("java -jar D:\\paper\\guo-summer_25\\test_7_11\\fdrbench-0.0.1\\fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:1 -miss_c 2 -enzyme 1 -maxLength 50",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}

plot_fdp_fdr_v2=function(fdp_file="",fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE,r=1,fixed_fdr_max=FALSE,max_x=NA,max_y=NA,
                         color_mapping=NULL,
                         legend_position=NULL,
                         fdr_decimal_place=2,
                         return_data=FALSE,
                         text_position=NULL,
                         add_max_qvalue=FALSE) {
    x <- read_csv(fdp_file)
    if("FDP_1B" %in% names(x)){
        if(r>=2){
            dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
                rename(`Combined method`=FDP,`Matched method`=FDP_1B,`Lower bound`=FDP_min) %>%
                gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
        }else{
            dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
                rename(`Combined method`=FDP,`Paired method`=FDP_1B,`Lower bound`=FDP_min) %>%
                gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
        }
    }else{
        dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_min) %>% distinct() %>% 
            rename(`Combined method`=FDP,`Lower bound`=FDP_min) %>%
            gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }
    

    max_fdp <- max(c(dat$FDP,dat$q_value))
    if(!is.null(fdr_max)){
        if(fixed_fdr_max){
            max_fdp <- fdr_max
        }else{
            max_fdp <- min(c(fdr_max,max_fdp))
        }
    }
    
    gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
            geom_abline(slope = 1,intercept = 0,color="gray")+
            #rasterise(geom_line(), dpi = 300) + 
            geom_line()+
            xlab("FDR threshold")+
            ylab("Estimated FDP")+
            theme_bw()+
            #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
            theme_pubr(base_size = 12,border = TRUE)

    if(!is.null(color_mapping)){
        gg1 <- gg1 + scale_color_manual(values = color_mapping)
    }

    if(scale_xy){

        if(!is.na(max_x) || !is.na(max_y)){
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
            if(!is.na(max_x)){
                gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
            }else{
                gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
            if(!is.na(max_y)){
                gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
            }else{
                gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
        }else{
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
                xlim(0,max_fdp)+
                ylim(0,max_fdp)+
                scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
                #scale_y_pct()+
                #scale_x_pct()+
                scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
        }
    }
    #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
    ## legend on the botton right
    ## no background color for the legend
    if(is.null(legend_position)){
        gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
    }else{
        gg1 <- gg1 + theme(legend.position = legend_position,plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
    }
    
            
    if(!is.null(fig_title)){
        gg1 <- gg1 + ggtitle(fig_title)
    }

    added_numbers <- NULL
    if(add_numbers){
        # add numbers on the top left size of the figure using annotation, text align to left
        # y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
        if(fdr_decimal_place==1){
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.1f%%",FDP001*100))
        }else if(fdr_decimal_place==2){
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.2f%%",FDP001*100))
        }else{
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
        }
        #n_t <- 
        if(abs(max_fdp-0.01)<=0.02){
            # text right align
            added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
            if(add_max_qvalue){
                added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
            }
            if(is.null(text_position)){
                gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
            }else{
                gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
            }
            
        }else{
            added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
            if(add_max_qvalue){
                added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
            }
            if(is.null(text_position)){
                gg1 <- gg1 + annotate("text", x = 0.01*1.05, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
            }else{
                gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
            }
        }
    }
   
    # library(ggpubr)
    # options(repr.plot.width = 12, repr.plot.height = 6)
    # gg <- ggarrange(gg1,gg2,ncol = 2,  common.legend = TRUE)
    # options(jupyter.plot_mimetypes = "image/png")
    if(return_data){
        return(list(gg=gg1,data=dat,added_numbers=added_numbers))
    }else{
        return(gg1)
    }

    #pdf("fdp_fdr_lib_qvalue_protein.pdf",width = 8,height = 4)
    #print(gg)
    #dev.off()
}


##Preprocessing
setwd("D:/paper/guo-summer_25/test_7_11/data/DDA-BERT")

input_dir <- "D:\\paper\\guo-summer_25\\test_7_11\\test_8_19\\data"
output_dir <- "D:\\paper\\guo-summer_25\\test_7_11\\test_8_19\\tsv"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file_path in csv_files) {
  tryCatch({

    cat(paste("Processing file:", file_path, "\n"))
    data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    
    required_columns <- c("cleaned_sequence", "precursor_mz", "precursor_charge", 
                         "modified_sequence", "label", "score", "q_value", 
                         "scan_number")
    
    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      cat(paste("Warning: File", file_path, "is missing the following columns:", paste(missing_columns, collapse = ", "), "\n"))
      next 
    }

    data$modified_sequence <- gsub("n\\[42\\]", "n(UniMod:1)", data$modified_sequence)
    data$modified_sequence <- gsub("N\\[.98\\]", "N(UniMod:7)", data$modified_sequence)
    data$modified_sequence <- gsub("Q\\[.98\\]", "Q(UniMod:7)", data$modified_sequence)
    data$modified_sequence <- gsub("M\\[15.99\\]", "M(UniMod:35)", data$modified_sequence)
    data$modified_sequence <- gsub("C\\[57.02\\]", "C(UniMod:4)", data$modified_sequence)
    
    data <- data[order(-data$score), ]
    
    data_unique <- data[!duplicated(data$scan_number), ]
    
    file_name <- basename(file_path)
    file_base <- tools::file_path_sans_ext(file_name)
    output_file <- file.path(output_dir, paste0(file_base, "_processed.tsv"))
    
    output_data <- data_unique[, required_columns]
    write.table(output_data, file = output_file, sep = "\t", na = "nan", 
                row.names = FALSE, quote = FALSE)
    
    cat(paste("Saved processed file to:", output_file, "\n"))
  }, error = function(e) {
    cat(paste("Error processing file", file_path, ":", conditionMessage(e), "\n"))
  })
}

cat("All files processed!\n")



##Figure7
color_mapping <- c("Paired method" = "#7CAE00", "Sample method" = "#C77CFF", "Lower bound" = "#00BFC4", "Combined method" = "#F8766D")

report_file <- "D:\\paper\\guo-summer_25\\test_7_11\\test_8_19\\tsv\\220627_GG_06_rep1_dil10_peptide_processed.tsv"
pep_file <- "D:\\paper\\guo-summer_25\\test_7_11\\human_entrapment_diann_pep.txt"

pro_fdp_file1 <- run_dda_bert_fdp_analysis(report_file, level = "peptide", prefix = "220627_GG_06_rep1_dil10", pep_file = pep_file, k_fold = 1)

fdp_data <- read.csv(pro_fdp_file1)

if("combined_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "combined_fdp"] <- "FDP"
} else {
  warning("Column 'combined_fdp' not found, cannot rename")
}

if("paired_fdp" %in% colnames(fdp_data)) {
  colnames(fdp_data)[colnames(fdp_data) == "paired_fdp"] <- "FDP_1B"
} else {
  warning("Column 'paired_fdp' not found, cannot rename")
}

temp_fdp_file <- tempfile(fileext = ".csv")
write.csv(fdp_data, temp_fdp_file, row.names = FALSE)

gg3 <- plot_fdp_fdr_v2(temp_fdp_file, fdr_max = 0.10, fig_title = "220627_GG_06_rep1_dil10", add_numbers = TRUE, r = 1, color_mapping = color_mapping, fdr_decimal_place = 2)

print(gg3)
ggsave("D:\\paper\\guo-summer_25\\test_7_11\\test_8_19\\png\\220627_GG_06_rep1_dil10.png", gg3, width = 4, height = 4, dpi = 300)
ggsave("D:\\paper\\guo-summer_25\\test_7_11\\test_8_19\\png\\220627_GG_06_rep1_dil10.pdf", gg3, width = 4, height = 4, dpi = 300)
