#==========================统计肿瘤-免疫互作数据===========================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search1_new/search1_final/")
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search1_new/search1_final/"
data <- read.table("search_detail.txt",sep='\t',header=T,as.is = T,fill = T,quote='')
dim(data)
#739 23
####species####
table(data$Species)
#house mouse       human  Norway rat 
#84         654           1 
# species_t <- table(data$Species)
# label = paste(names(species_t),"(",unname(species_t),")",sep='')
# colour <- c("#fafad2","#ffb6c1","#90ee90")
# names(colour) <- names(species_t)
# 
# pdf("species.pdf")
# pie(species_t,labels = label,col=colour,main="species")
# dev.off()
####不分物种####
data_t <- table(data$Interaction_Type1)
# Direct interaction        Indirect interaction Interaction through exosome 
#186                 440                          113 
label = paste(names(data_t),"(",unname(data_t),")",sep='')
colour <- c("#dda0dd","#f5deb3","#ffa07a")
names(colour) <- names(data_t)

pdf(paste0(d,"/all_species.pdf"))
pie(data_t,labels = label,col = colour)
dev.off()
####human####
human <- data[data$Species=="human",]
human_t <- table(human$Interaction_Type1)
# Direct interaction        Indirect interaction Interaction through exosome 
#169                 389                          99 
label = paste(names(human_t),"(",unname(human_t),")",sep='')
colour <- c("#db7093","#ffa500","#48d1cc")
names(colour) <- names(human_t)

pdf("human.pdf")
pie(human_t,labels = label,col = colour,main="human")
dev.off()
####mouse####
mouse <- data[data$Species=="house mouse",]
mouse_t <- table(mouse$Interaction_Type1)
# Direct interaction        Indirect interaction Interaction through exosome 
#19                 52                          14 
label = paste(names(mouse_t),"(",unname(mouse_t),")",sep='')
colour <- c("#db7093","#ffa500","#48d1cc")
names(colour) <- names(mouse_t)

pdf("mouse.pdf")
pie(mouse_t,labels = label,col = colour,main="mouse")
legend("right",legend=names(mouse_t),pch=15,col=c("#db7093","#ffa500","#48d1cc"))
dev.off()
####immune cells####
cell <- unique(names(table(data$Cell)),names(table(data$Cell.1)))
immune_cell <- setdiff(cell,c("cancer cell","other cell"))

c_matrix <- matrix(0,nrow = length(immune_cell),ncol=3)
rownames(c_matrix) <- immune_cell
colnames(c_matrix) <- c("Direct interaction","Indirect interaction","Interaction through exosome")

for(i in immune_cell){
  imm_data <- data[data$Cell %in% i|data$Cell_1 %in% i,,drop=F]
  interact_type <- imm_data$Interaction_Type1
  interact_type_count <- table(interact_type)
  c_matrix[i,names(interact_type_count)] <- unname(interact_type_count)
}
c_plot <- data.frame(cell=rep(rownames(c_matrix),times=3),
                     type=rep(colnames(c_matrix),each=nrow(c_matrix)),
                     num=as.numeric(c_matrix))
write.table(c_plot,paste0(d,"/immune_count.txt"),sep='\t',quote = F,row.names = F)
library(ggplot2)
colour <- c("#dda0dd","#f5deb3","#ffa07a")
names(colour) <- colnames(c_matrix)

pdf(paste0(d,"/immune_bar.pdf"))
ggplot(c_plot,aes(x=cell,y=num,fill=type))+geom_bar(stat = "identity")+
  scale_fill_manual(values=colour)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,size=14,hjust = 1,colour = "black"),
        axis.text.y = element_text(size=12,colour = "black"),
        axis.title.y = element_text(size=14,colour = "black"),
        panel.grid = element_blank())+
  labs(x="",y="number of interactions")
dev.off()
####癌症类型####
table(data$Cancer)
library(plyr)
library(ggplot2)

cancer_count <- ddply(data,.(Cancer,Interaction_Type1),nrow)
write.table(cancer_count,"cancer_count.txt",sep="\t",quote = F,row.names = F)
cols <- c("#dda0dd","#f5deb3","#ffa07a")
names(cols) <- c("Direct interaction","Indirect interaction","Interaction through exosome")
pdf("cancer_count.pdf",width = 8,height = 12)

cc <- ddply(data,.(Cancer),nrow)
cancer_sort <- cc[order(cc$V1),]$Cancer
ggplot(cancer_count,aes(x=factor(Cancer,levels = cancer_sort),y=V1,fill=Interaction_Type1))+geom_bar(stat = "identity",width =0.8)+coord_flip()+
  theme(axis.text.x = element_text(colour = "black",size = 10,hjust = 1),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = cols)+
  labs(y="Number of interactions")
dev.off()

cancer_count <- ddply(data,.(Cancer,Cancer_Type,Interaction_Type1),nrow)
cc_unin <- unique(cancer_count[,c(1,2)])

res <- matrix(0,nrow=nrow(cc_unin),ncol = 3)
colnames(res) <- c("Direct interaction","Indirect interaction","Interaction through exosome")
for(i in 1:nrow(cc_unin)){
  cc <- cancer_count[cancer_count$Cancer == cc_unin[i,1]&cancer_count$Cancer_Type == cc_unin[i,2],,drop=F]
  res[i,cc[,3]] <- cc[,4]
}
res2 <- cbind(cc_unin,res)
res3 <- ddply(res2,.(Cancer),function(x) colSums(x[,c(3,4,5)]))
res4 <- c()
for(j in 1:nrow(res3)){
  ct <- unique(cancer_count[cancer_count$Cancer == res3[j,]$Cancer,]$Cancer_Type)
  res4 <- rbind(res4,cbind(res3[j,],ct))
}

write.table(res4,"Experiment_TIC.txt",sep='\t',quote = F,row.names = F)
