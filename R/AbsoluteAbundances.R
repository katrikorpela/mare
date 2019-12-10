AbsoluteAbundances <- function(qPCR.file, qPCR.var, folder.name = "", copynumberDB = NULL) {
  
qpcr <- read.delim(qPCR.file)
if(length(copynumberDB)>0) rrn <- read.delim(copynumberDB)

for (j in grep(pattern="read",x=list.files(path=paste(folder.name,"TaxonomicTables",sep=""),pattern = "table.txt"),value=T,invert=T)){
table <- read.delim(paste(folder.name,"TaxonomicTables/",j,sep=""),check.names = F)
rownames(table)<-table[,1]
table <- t(table[,-1])
rtable <- table/rowSums(table)
rtable[rowSums(table)==0,]<-0
rtableq <- rtable*qpcr[,qPCR.var]

if(length(copynumberDB)>0){
rtableqc <- rtableq
tableqname <- sapply(colnames(rtableq),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
for(i in setdiff(tableqname,rrn$name)) suppressWarnings(tryCatch(tableqname[tableqname==i]<- unlist(sapply(colnames(rtableq)[tableqname==i],function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])-1])),error=function(e) NULL))
for(i in setdiff(tableqname,rrn$name)) suppressWarnings(tryCatch(tableqname[tableqname==i]<- unlist(sapply(colnames(rtableq)[tableqname==i],function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])-2])),error=function(e) NULL))
for(i in setdiff(tableqname,rrn$name)) suppressWarnings(tryCatch(tableqname[tableqname==i]<- unlist(sapply(colnames(rtableq)[tableqname==i],function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])-3])),error=function(e) NULL))
for(i in setdiff(tableqname,rrn$name)) suppressWarnings(tryCatch(tableqname[tableqname==i]<- unlist(sapply(colnames(rtableq)[tableqname==i],function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])-4])),error=function(e) NULL))
for(i in setdiff(tableqname,rrn$name)) suppressWarnings(tryCatch(tableqname[tableqname==i]<- unlist(sapply(colnames(rtableq)[tableqname==i],function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])-5])),error=function(e) NULL))
for(i in 1:ncol(rtableq))  rtableqc[,i] <- rtableq[,i]/rrn$mean[rrn$name==tableqname[i]][1]
rtableqc[,colnames(rtableqc)[is.na(colSums(rtableqc))]]<-rtableq[,colnames(rtableqc)[is.na(colSums(rtableqc))]]
rtableq <- rtableqc
}

tablef <- as.data.frame(t(rtableq))
tablef$taxon <- rownames(tablef)
tablef <- tablef[,c(ncol(tablef),2:(ncol(tablef)-1))]
write.table(tablef,paste(folder.name,"TaxonomicTables/absolute_",j,sep=""),sep="\t",quote=F,row.names = F)
}
  
  }