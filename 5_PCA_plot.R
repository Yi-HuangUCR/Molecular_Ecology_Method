library(factoextra)
sample <- read.table(“PCAinput.txt”)
factor <- read.csv(“population.csv”)
m <- procomp(sample[,2:35])
PCAdata <- data.frame(factor$pop,factor$subspecies,m$x)

#color code individual sample based on populations
library(ggplot2)
g <- ggplot(data = PCAdata, aes(x=PC1, y=PC2, color = pop))
+ geom_point()
+ stat_ellipse()
+ theme_classic()

#color code individual sample based on subspecies classification
library(ggplot2)
g <- ggplot(data = PCAdata, aes(x=PC1, y=PC2, color = subspecies))
+ geom_point()
+ stat_ellipse()
+ theme_classic()
