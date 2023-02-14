
library(ggplot2)
library(ggseqlogo)
library(Biostrings)
filepath<-"C:\\Users\\Rongze chen\\Documents\\WeChat Files\\wxid_1avvzl2v7e6922\\FileStorage\\File\\2021-07\\UP000002717_1140.fasta" #根据实际情况进行修改
setwd(testpath)
seqdata = readAAStringSet("C:\\Users\\Rongze chen\\Documents\\WeChat Files\\wxid_1avvzl2v7e6922\\FileStorage\\File\\2021-07\\UP000002717_1140.fasta\\UP000002717_1140.fasta", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
names_seq <- names(seqdata)
length_seq <- width(seqdata)

AA_seq <- c("Y","W","F")
AA_seq_Frequency <- letterFrequency(seqdata,AA_seq)
AA_seq_Frequency <- as.data.frame(AA_seq_Frequency)
AA_seq_Frequency <- cbind(names_seq = names_seq,length_seq = length_seq,AA_seq_Frequency)
Name_group <- c("A","seq","note")
meta_data <-  as.data.frame(tstrsplit(AA_seq_Frequency$names_seq,"[|]"),col.names = Name_group)
AA_seq_Frequency <- cbind(meta_data,AA_seq_Frequency[,-1])
Name_group <- c("note","describe")
meta_data_note <-  as.data.frame(tstrsplit(AA_seq_Frequency$note," OS"),col.names = Name_group)
AA_seq_Frequency$note <- meta_data_note$note
AA_seq_Frequency$describe <- meta_data_note$describe

AA_seq_Frequency <- mutate(AA_seq_Frequency,
                           Y_freq = AA_seq_Frequency$Y/AA_seq_Frequency$length_seq,
                           W_freq = AA_seq_Frequency$W/AA_seq_Frequency$length_seq,
                           F_freq = AA_seq_Frequency$F/AA_seq_Frequency$length_seq,
                           total_freq = Y_freq+W_freq+F_freq)
fwrite(AA_seq_Frequency,paste(filepath,"\\","AA_seq_Frequency",".txt",sep=""),
       row.names=F,col.names=T,quote=F,sep = "\t")


AA_seq_group <- c("Photosystem II","Photosystem I","quinone oxidoreductase","ATP synthase")
AA_seq_Frequency_all <- row_number_all <- c()
for ( i in AA_seq_group)
{
  row_number <- grep(i,AA_seq_Frequency$note)
  row_number_all <- c(row_number_all,row_number)
  AA_seq_Frequency_i <- AA_seq_Frequency[grep(i,AA_seq_Frequency$note),]
  AA_seq_Frequency_i$group <- i
  AA_seq_Frequency_all <- rbind(AA_seq_Frequency_all,AA_seq_Frequency_i)
}
AA_seq_Frequency_other <- AA_seq_Frequency[-row_number_all,]
AA_seq_Frequency_other$group <- "Other group"
AA_seq_Frequency_all_1 <- AA_seq_Frequency
AA_seq_Frequency_all_1$group <- "All"
AA_seq_Frequency_all <- rbind(AA_seq_Frequency_all,AA_seq_Frequency_other)
AA_seq_Frequency_all <- rbind(AA_seq_Frequency_all,AA_seq_Frequency_all_1)

AA_seq_Frequency_all$group <- gsub("quinone oxidoreductase","NDH complex",AA_seq_Frequency_all$group)
AA_seq_Frequency_all$group <- factor(AA_seq_Frequency_all$group,levels=c("Photosystem I","Photosystem II","NDH complex","ATP synthase","Other group","All"))
AA_seq_Frequency_all <- filter(AA_seq_Frequency_all,group != "Other group")
plot<-ggplot(data=AA_seq_Frequency_all,aes(x=group,y=total_freq,color= group))+
  geom_point(alpha = 1,position = position_jitter(0.15),size=1)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_bw()+ labs(y = "Total frequency")+#
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title =  element_blank(),
    legend.position = "nane",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90,family = "Times"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15,family = "Times"),
    axis.text.x = element_text(size = 15,angle = 45,family = "Times",hjust = 1,vjust = 1),
    axis.ticks = element_line(size = 1),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"))+
  theme(axis.title= element_text(family = "Times" ,size=20))
plot 
ggsave(filename=paste(filepath,"//","AA_seq_Frequency_all_bloxplot.png", sep=""),plot=plot,
       limitsize=T,width = 5,height = 5)
