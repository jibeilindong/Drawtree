
library(ggpubr)
library(ggplot2)
library(reshape2)
library(extrafont)
loadfonts(device = "pdf")
library(ggshapes)
library(ape)
library(ggdendro)
library(dendextend)
library(RColorBrewer)
library(scales)

testpath <- outpath <- 'C://XXX' #change your work dir
setwd(testpath) 

df <- read.csv("XXX.txt",header = TRUE,skip=0,stringsAsFactors=F,row.names=NULL,sep = "\t") #change your txt file
colnames(df)[2] <- "V2"
df <- df[1:23,]
df_gid <- df$V2

rawdata <- df[, 4:dim(df)[2]]
rawdata <- rawdata[,-45]
rownames(rawdata) <- df[,2]

D <- dist(rawdata, method = "euclidean") 
a <- as.matrix(rawdata)
D <- a%*%t(a)
D_1 <- apply(D,c(1,2),function(point) {1/(point+1)})
D_1 <- as.dist(D_1)
hc <- hclust(D_1, method="average")

colors <- hue_pal()(4)
rlabels <- hc$labels[hc$order]
ro <- apply(df, 1, function(row){which(row[1]==rlabels)})
rMIC <- as.factor(df[ro,c("V2")])

df$V2 <- factor(x=df$V2,levels=df$V2[hc$order],ordered=T)
df2 <- df[hc$order,c(-3,-48)]

data <- melt(df2,id.vars=c("V2","X"))
data$variable <- gsub("X","syf",data$variable)
my_cols <- c("white", "red")

colMIC <- colors[as.factor(data$No)]

xx <- hc %>% as.dendrogram %>% set("labels_col", value=colors[1], k=23)
gg0 <- ggdendrogram(xx, rotate=T, labels=FALSE, leaf_labels = TRUE) + 
  scale_y_reverse() +
  theme(text=element_text(size=9, family="Times"),
        axis.text.x = element_blank(),
        plot.margin = margin(0, -2.5, 0.8, 0, "cm")) #plot.margin = margin(-0.05, -2.5, 1.25, 0, "cm")
gg0

gg1 <- ggballoonplot(data,x="variable",y="V2",size=3,fill="value",shape = 22)+
  scale_fill_gradientn(colors = my_cols) + theme_bw() + 
  scale_x_discrete(position = "bottom") + scale_y_discrete(position = "left") +
  theme(text=element_text(size=9, family="Times"), legend.position="none", panel.spacing=unit(0.7,"lines"),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.placement="outside", strip.text = element_text(face="bold", size=6, lineheight=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black",angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("")
gg1
plot_tree <- ggarrange(gg0, gg1, widths = c(1,1.5), nrow=1)
plot_tree
ggsave(filename=paste(outpath,"//","plot_tree.pdf", sep=""),plot=plot_tree,limitsize=T,width = 10,height = 4) # change your program name
