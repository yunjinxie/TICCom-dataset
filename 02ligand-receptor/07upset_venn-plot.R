############根据LR database，绘制upset图
library(UpSetR)

setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search2_new/")
d <- getwd()

data <- read.table("union_source_curated.txt",sep='\t',header=T,as.is = T)
data <- as.matrix(data)
head(data)
data2_h <- data[,c("cellchatDB","celltalkDB","icellnetDB",
                 "iTALKDB","nichenetDB","ramilowskiDB",
                 "singlecellsignalRDB")]

data2_h[data2_h[,]=="Yes"] <- 1
data2_h[data2_h[,]=="-"] <- 0

data3_h <- matrix(as.numeric(data2_h),nrow=nrow(data2_h))
colnames(data3_h) <- c("cellchatDB","celltalkDB","icellnetDB",
                       "iTALKDB","nichenetDB","ramilowskiDB",
                       "singlecellsignalRDB")

data3_h <- as.data.frame(data3_h)

###########统计有多少特异的
count <- c()
for(i in 1:ncol(data3_h)){
  index <- c()
  for(j in 1:nrow(data3_h)){
    index <- c(index,which(data3_h[j,i]==1&sum(data3_h[j,-i]==1)==0))
  }
  ind <- length(index)
  count <- c(count,ind)
}

#############
rs1 <- c()

for(k in 2:ncol(data3_h)){
  group <- combn(ncol(data3_h),k)
  
  for(n in 1:ncol(group)){
    
    index <- group[,n]
    col_method <- colnames(data3_h)[index]

    rs <- apply(data3_h,1,function(x) length(which(length(which(x[index]==1))==length(index)&sum(x[-c(index)]==1)==0)))
    
    len <- length(which(rs==1))
    
    if(len!=0){
      rs1 <- rbind(rs1,cbind(paste(col_method,collapse = ","),length(index),len))
    }
  }
} 
count2 <- rbind(cbind(colnames(data3_h),1,count),rs1)
colnames(count2) <- c("method_nam","method_num","count")
write.table(count2,"count.txt",sep='\t',quote = F,row.names = F)
################画饼图
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search2")
count2 <- read.table("count.txt",sep='\t',header=T,as.is = T)
count2 <- as.matrix(count2)
count3 <-count2[count2[,"method_num"]!="1",] 
num<- as.numeric(count3[,"count"])
c_plot <- aggregate(num,list(count3[,"method_num"]),function(x) sum(x))

num1 <- count2[count2[,"method_num"]==1,c(1,3)]
num1[,2] <- as.numeric(num1[,2])
num2 <- num1[num1[,2]!=0,]
colnames(num2) <- NULL
c_plot <- as.matrix(c_plot)
colnames(c_plot) <- NULL
c_plot2 <-  rbind(num2,c_plot)
num_plot <- as.numeric(c_plot2[,2])
per <- round(num_plot/sum(num_plot),4)

pdf("pie.pdf",width = 20,height = 7)
par(mai=c(.5,.5,.5,4))
pie(num_plot,labels=paste(c_plot2[,1],": ",c_plot2[,2],"(",per,")",sep=''),
    col=colorRampPalette(c("#fdefec","#ffab88", "#ff7f50"))(15),xpd=T,cex=1.5)
legend(1.9,1,legend =c_plot2[,1],col = colorRampPalette(c("#fdefec","#ffab88", "#ff7f50"))(15),ncol=1,pch = 15,cex=1.2)
dev.off()
############
pdf("human_upset.pdf",width=20,height=8)
upset(data3_h,nsets=7,order.by="freq",
      main.bar.color = "green4",sets.bar.color = "skyblue3",
      point.size=2,text.scale=c(2,1.5,2,1.5,2,2),sets.x.label="",
      queries = list(list(query  = intersects,params = list("icellnetDB","cellchatDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB"),
                          color = "orange2",active=TRUE)))
dev.off()
######显示大于50的交集
pdf("human_upset2.pdf",width=13,height=7)
upset(data3_h,nsets=7,order.by="freq",
      intersections=list(list("icellnetDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB"),
                         list("celltalkDB","nichenetDB"),
                         list("ramilowskiDB","iTALKDB","singlecellsignalRDB"),
                         list("cellchatDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB"),
                         list("cellchatDB","nichenetDB"),
                         list("singlecellsignalRDB","celltalkDB"),
                         list("icellnetDB","cellchatDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB"),
                         list("singlecellsignalRDB","celltalkDB","nichenetDB"),
                         list("cellchatDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB"),
                         list("ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB"),
                         list("ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB")
                         ),
      main.bar.color = "green4",sets.bar.color = "skyblue3",point.size=2,text.scale=c(2,1.5,2,1.5,2,2),sets.x.label="",
      queries = list(list(query  = intersects,params = list("icellnetDB","cellchatDB","ramilowskiDB","iTALKDB","singlecellsignalRDB","celltalkDB","nichenetDB"),
                          color = "orange2",active=TRUE)))
dev.off()
##===========================house mouse============================
data2_m <- data[,c("Gene1Id","Gene2Id","celltalkDB_mouse","RNAMagnetDB_mouse")]
pair_ID <- paste(data2_m[,"Gene1Id"],data2_m[,"Gene2Id"],sep='_')

#######画一个Venn图
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    celltalkDB = pair_ID[data2_m[,"celltalkDB_mouse"]=="Yes"],
    RNAMagnetDB = pair_ID[data2_m[,"RNAMagnetDB_mouse"]=="Yes"]
  ),filename = NULL,lwd=3,fill = c("cornflowerblue", "darkorchid1"),
  label.col = "white",
  cex=3,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20),
  cat.cex=4,
  main = "House Mouse",
  main.cex = 2,
  inverted = F)

pdf("house_mouse_venn.pdf",height = 7,width=9)
grid.draw(venn.plot)
dev.off()

