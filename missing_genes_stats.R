.libPaths("D:/Work/r_libs")
.libPaths("/mnt/d/Work/r_libs_ubuntu")
.libPaths("/data/meyer/viz/r_libs")
library(tidyverse)
library(reshape2)
library(ggplot2)

file_paths <- read.csv("MISSING_GENES_FINAL", header = F)
file_paths <- split(file_paths$V2, file_paths$V1)
org_list <- c()

org_list <- lapply(file_paths, FUN = function(x){
    factor(as.vector(t(read.table(x,header = F))))
})

org_list <- melt(org_list)

stats <- apply(org_list, 2, table)

number_of_bar <- nrow(stats$value)

stats$value <- as.data.frame(cbind(stats$value, seq(1,number_of_bar)))

label_data <- stats$value

angle <-  90 - 360 * (label_data$V2-0.5) /number_of_bar

label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

p <- ggplot(stats$value, aes(x=as.factor(V2), y=V1)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-80,50) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=V2, y=V1+10, label=rownames(label_data), hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 
