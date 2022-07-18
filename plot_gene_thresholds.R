require(ggplot2)
require(scales)
require(hrbrthemes)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

param_file <- "parameters.txt"
param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^") #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")
plots_out_path <- as.character(param_table[which(param_table=="plot_path"),c(2)])
gene_drop_thresh <- as.numeric(param_table[which(param_table=="gene_drop_thresh"),c(2)])

dir.create(plots_out_path)
gene.thresholds <- read.csv("files/gene_thresholds.txt", header = F)
rownames(gene.thresholds) <- gene.thresholds$V1
#gene.thresholds <- gene.thresholds[,-1]

gene.thresholds$V5 <- rescale(gene.thresholds[, 4])
gene.thresholds$V6 <- gene.thresholds$V4
gene.thresholds$V4 <- range01(gene.thresholds$V4)

padding <- 5
gene.thresholds$V5 <- gene.thresholds$V5 * max(gene.thresholds$V6)/min(gene.thresholds$V6) + padding
gene.thresholds$right <- cumsum(gene.thresholds$V5) + 30*c(0:(nrow(gene.thresholds)-1))
gene.thresholds$left <- gene.thresholds$right - gene.thresholds$V5 

# Plot
p <- ggplot(gene.thresholds, aes(ymin = 0)) + 
  geom_rect(aes(xmin = left - 0.5 , xmax = right +0.5 , ymax = V2, fill = V3)) + # , color = V3)) + 
  xlab("") + 
  ylab("Threshold") +
  theme_ipsum() +
  theme(axis.text.x = element_blank(), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  geom_label(data=gene.thresholds, aes(x=right - V5 / 2, y=V2+2, label=V2), color="black", fontface="bold",alpha=0.6, size=2.7, inherit.aes = FALSE, angle = 60 ) +
  geom_text(data=gene.thresholds, aes(x=right - V5 / 2, y= -10, label= paste(V1,"(",V6,")")), color="black", fontface="bold",alpha=0.6, size=2.5, inherit.aes = FALSE, angle = 90 ) +
  geom_hline(yintercept=gene_drop_thresh)

p$labels$fill <- "Orthologous between reference species?"

png(filename =paste(plots_out_path,"/gene_thresholds.png",sep="/"),width = 15, height = 15 , units = "in", res = 100)
p
dev.off()

