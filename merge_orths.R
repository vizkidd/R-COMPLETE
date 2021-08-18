suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(require(showtext))
suppressMessages(require(ggplot2))
#library(plotly)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Org list path (2) Path to *.orths files (3) Output file name (4) E-value cut-off", call.=FALSE)
}

set.seed(123)

font_add_google("Gochi Hand", "gochi")
#font_add_google("Schoolbell", "bell")

## Automatically use showtext to render text
showtext_auto()

org_list_path <- args[1]
orths_path <- args[2] ##NAMING format of each file should be query-subject.orths
out_file <- args[3]
e_cutoff <- as.numeric(args[4])

print(e_cutoff)

orgs.list <- factor(scan(org_list_path, character()))
subject_list <- list()
for (query in orgs.list) {
  for (subject in orgs.list) {
    if(query != subject){
      file0 <- file.path(orths_path, paste(query,"-",subject,".orths", sep=""))
      #print(file0)
      if(file.exists(file0) && file.info(file0)$size > 0){
        subject_list[[query]][[subject]] <- read.table(file0)
      }else{
        print(paste("WARNING:",file0," doesn't exist"))
      }
    }
  }
}

final_df <- data.frame()
for (query in orgs.list) {
  for (subject in orgs.list) {
    if(query != subject){
      final_df <- rbind(final_df, rbind(subject_list[[query]][[subject]][na.omit(match(subject_list[[subject]][[query]][,c(2)],subject_list[[query]][[subject]][,c(1)])),], subject_list[[subject]][[query]][na.omit(match(subject_list[[subject]][[query]][,c(1)],subject_list[[query]][[subject]][,c(2)])),]))
    }
  }
}
final_df <- final_df[unique(match(final_df$V1,final_df$V2)),] ##ONLY KEEP HITS WHICH have reciprocal hits
ggplot(final_df, aes(x=V3)) +   geom_histogram(color="black", fill="white")  ## % identity histogram
ggplot(final_df, aes(x=V11)) +   geom_histogram(color="black", fill="white")  + scale_x_continuous(limits = c(NA,max(final_df$V11))) ## E-value pre-filter histogram
nrow(final_df)
final_df <- final_df[which(final_df$V11 < e_cutoff),] ## E-value filtering  (1e-05)
nrow(final_df)
ggplot(final_df, aes(x=V11)) +   geom_histogram(color="black", fill="white") # + scale_x_continuous(limits = c(NA,max(final_df$V11)/1000)) ## POST FILtering e-value hist ,  binwidth=min(final_df$V11)
ggplot(final_df, aes(x=V3)) +   geom_histogram(color="black", fill="white")  ##  POST FILtering % identity histogram

write.table(final_df, file=out_file,quote = F,row.names = F, col.names = F)