TimeSeries <- function(timepoints, species.table=NULL, genus.table=NULL, family.table=NULL, order.table=NULL, class.table=NULL, phylum.table=NULL, 
                      meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
                      subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, 
                      pdf = F,  min.prevalence = 0.1, min.abundance = 0, 
                      label.direction = 1,  nonzero = T, relative = T, 
                      time, log.time = F){


cm <- read.delim(meta)
cm$group <- cm[,group]
result <- list()

for(i in timepoints){
  result[[i]] <-   GroupTest(genus.table = genus.table,
                             family.table = family.table,
                             order.table = order.table,
                             class.table = class.table,
                             phylum.table = phylum.table,
                             meta=meta,
                             group=group,
                             compare.to = compare.to,
                             readcount.cutoff = readcount.cutoff,
                             confounders = confounders,
                             select.by = time,select = i,
                             keep.result = T, nonzero = nonzero,
                             subject.ID = subject.ID, outlier.cutoff = outlier.cutoff, p.cutoff = p.cutoff, 
                             min.prevalence = min.prevalence, min.abundance = min.abundance, pdf = pdf,
                             label.direction = label.direction, relative = relative)  
}


library(scales)
library(ggplot2)
library(reshape2)

colorlist <- c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", "#A65628", "#F781BF", "#999999")[1:length(levels(as.factor(cm$group)))]
grouplist <- levels(as.factor(cm$group))

sig <- list()
for(i in timepoints){
sig[[i]]<- result[[i]][apply(na.omit(result[[i]][,grep(names(result[[i]]),pattern="_p")]),MARGIN = 1, FUN = min)<p.cutoff,"taxon"]
}
sigall <- unique(unlist(sig))
siglevel <- sapply(sigall,function(x) length(strsplit(x,split="_")[[1]]))
 
if(relative) yl <- "Relative abundance" else yl <- "Abundance"

if(length(genus.table)>0){
rgen <- sigall[siglevel==5]

lgen <- length(rgen)

if(lgen>0){
  
if(lgen < 5) {
  ncols <- lgen
  nrows <- 1
} else {
ncols <-  round(sqrt(lgen)) + 1
nrows <-  floor(sqrt(lgen))
}

gen <- read.delim(genus.table)

colnames(gen)[grepl(pattern="incertae_sedis",x=colnames(gen))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(gen)[grepl(pattern="incertae_sedis",x=colnames(gen))])
colnames(gen)[grepl(pattern="Incertae_Sedis_",x=colnames(gen))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(gen)[grepl(pattern="Incertae_Sedis_",x=colnames(gen))])
#colnames(gen)[grepl(pattern="Erysipelotrichi_",x=colnames(gen))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(gen)[grepl(pattern="Erysipelotrichi_",x=colnames(gen))])

for(i in grep(pattern="_IncertaeSedis",x=colnames(gen))){
  colnames(gen)[i] <- gsub(x=colnames(gen)[i] ,pattern="_IncertaeSedis",fixed=T,
replacement = paste("_",strsplit(x=colnames(gen)[i],split = "_")[[1]][4],strsplit(x=colnames(gen)[i],split = "_")[[1]][5],sep=""))
}
    
for(i in grep(pattern="_incertaesedis",x=colnames(gen))){
  colnames(gen)[i] <- gsub(x=colnames(gen)[i] ,pattern="_incertaesedis",fixed=T,
replacement = paste("_",strsplit(x=colnames(gen)[i],split = "_")[[1]][4],strsplit(x=colnames(gen)[i],split = "_")[[1]][5],sep=""))
}
    
for(i in grep(pattern="_uncultured",x=colnames(gen))){
  colnames(gen)[i] <- gsub(x=colnames(gen)[i] ,pattern="_uncultured",fixed=T,
replacement = paste("_",strsplit(x=colnames(gen)[i],split = "_")[[1]][4],strsplit(x=colnames(gen)[i],split = "_")[[1]][5],sep=""))
}

if(relative) gen <- (gen+1)/rowSums(gen)
res <- as.data.frame(gen[,rgen])
res$group <-cm$group
if(log.time) res$time <- log(cm[,time]+1) else res$time <- cm[,time]

if(lgen>1){
resagr <- aggregate((res[,rgen]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,rgen]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
} else {
resagr <- aggregate((res[,1]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,1]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
names(resagr)[3]<- rgen
names(ressd)[3]<- rgen
}



df <- reshape2::melt(resagr,id=c('time','group'))
dfsd <- reshape2::melt(ressd,id=c('time','group'))
df$min <- df$value-dfsd$value
df$max <- df$value+dfsd$value
df$varname <- sapply(as.character(df$variable),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
df <- na.omit(df)
df <- df[df$time %in% timepoints,]

for(k in unique(df$group)[unique(df$group)!=compare.to]){
  df[,paste(k,"p",sep="_")] <- NA
  df[,paste("star",k,sep="_")] <- NA
  for(i in timepoints){
  for(j in intersect(df$variable,result[[i]]$taxon)){
  df[df$time==i&df$variable==j&df$group==k,paste(k,"p",sep="_")] <- try(result[[i]][result[[i]]$taxon==j,paste(k,"_p",sep="")])
  }}
  df[df[,paste(k,"p",sep="_")]<0.05&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "*"
  df[df[,paste(k,"p",sep="_")]<0.01&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "**"
  df[df[,paste(k,"p",sep="_")]<0.001&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "***"
  }

p <- ggplot(df, aes(y=value,x=time,group=as.factor(group),color=as.factor(group)))+
  facet_wrap(~varname,ncol=ncols,scales='free')+
  geom_line(position = position_dodge(width = 0.2))+
 geom_point(position = position_dodge(width = 0.2))+
  scale_color_manual(values=colorlist,name="Group")+
#scale_y_continuous(labels = trans_format("exp", format = scientific_format(digits=3)))+
geom_linerange(position = position_dodge(width = 0.2),
                 aes(x=time,group=as.factor(group),color=as.factor(group),ymin=min,ymax=max),lwd=0.2)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.background = element_rect(colour="white", fill="white"))+
 xlab(time)+
  ylab(yl)

for(k in unique(df$group)[unique(df$group)!=compare.to]){
p <- p + geom_text(aes_(label=(df[,paste("star",k,sep="_")])), 
                   colour=colorlist[grouplist==k], vjust=0, 
                   position = position_dodge(width = 0.2),size=5)
}

if(pdf) {pdf(width=ncols*4,height=nrows*3,paste(group,compare.to,"TimeSeries_genus.pdf",sep="_")); plot(p); dev.off()} 
quartz(width=ncols*4,height=nrows*3); plot(p)
}
}

#plot family level
if(length(family.table)>0){
rfam <- sigall[siglevel==4]

lfam <- length(rfam)

if(lfam>0){

if(lfam < 5) {
  ncols <- lfam
  nrows <- 1
} else {
ncols <-  round(sqrt(lfam)) + 1
nrows <-  floor(sqrt(lfam))
}

fam <- read.delim(family.table)

#colnames(fam)[grepl(pattern="Erysipelotrichi_",x=colnames(fam))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(fam)[grepl(pattern="Erysipelotrichi_",x=colnames(fam))])

if(relative) fam <- (fam+1)/rowSums(fam)
res <- as.data.frame(fam[,rfam])
res$group <-cm$group
if(log.time) res$time <- log(cm[,time]+1) else res$time <- cm[,time]

if(lfam>1){
resagr <- aggregate((res[,rfam]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,rfam]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
} else {
resagr <- aggregate((res[,1]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,1]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
names(resagr)[3]<- rfam
names(ressd)[3]<- rfam
}


df <- reshape2::melt(resagr,id=c('time','group'))
dfsd <- reshape2::melt(ressd,id=c('time','group'))
df$min <- df$value-dfsd$value
df$max <- df$value+dfsd$value
df$varname <- sapply(as.character(df$variable),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
df <- na.omit(df)
df <- df[df$time %in% timepoints,]

for(k in unique(df$group)[unique(df$group)!=compare.to]){
  df[,paste(k,"p",sep="_")] <- NA
  df[,paste("star",k,sep="_")] <- NA
  for(i in timepoints){
  for(j in intersect(df$variable,result[[i]]$taxon)){
  df[df$time==i&df$variable==j&df$group==k,paste(k,"p",sep="_")] <- try(result[[i]][result[[i]]$taxon==j,paste(k,"_p",sep="")])
  }}
  df[df[,paste(k,"p",sep="_")]<0.05&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "*"
  df[df[,paste(k,"p",sep="_")]<0.01&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "**"
  df[df[,paste(k,"p",sep="_")]<0.001&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "***"
  }
  
p <- ggplot(df, aes(y=value,x=time,group=as.factor(group),color=as.factor(group)))+
  facet_wrap(~varname,ncol=ncols,scales='free')+
  geom_line(position = position_dodge(width = 0.2))+
 geom_point(position = position_dodge(width = 0.2))+
  scale_color_manual(values=colorlist,name="Group")+
#scale_y_continuous(labels = trans_format("exp", format = scientific_format(digits=1)))+
geom_linerange(position = position_dodge(width = 0.2),
                 aes(x=time,group=as.factor(group),color=as.factor(group),ymin=min,ymax=max),lwd=0.2)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.background = element_rect(colour="white", fill="white"))+
 xlab(time)+
  ylab(yl)

for(k in unique(df$group)[unique(df$group)!=compare.to]){
p <- p + geom_text(aes_(label=(df[,paste("star",k,sep="_")])), 
                   colour=colorlist[grouplist==k], vjust=0, 
                   position = position_dodge(width = 0.2),size=5)
}

if(pdf) {pdf(width=ncols*4,height=nrows*3,paste(group,compare.to,"TimeSeries_family.pdf",sep="_")); plot(p); dev.off()} 
quartz(width=ncols*4,height=nrows*3); plot(p)
}
}

#plot order level
if(length(order.table)>0){
rord <- sigall[siglevel==3]

lord <- length(rord)

if(lord>0){
  
if(lord < 5) {
  ncols <- lord
  nrows <- 1
} else {
ncols <-  round(sqrt(lord)) + 1
nrows <-  floor(sqrt(lord))
}


ord <- read.delim(order.table)

#colnames(ord)[grepl(pattern="Erysipelotrichi_",x=colnames(ord))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(ord)[grepl(pattern="Erysipelotrichi_",x=colnames(ord))])

if(relative) ord <- (ord+1)/rowSums(ord)
res <- as.data.frame(ord[,rord])
res$group <-cm$group
if(log.time) res$time <- log(cm[,time]+1) else res$time <- cm[,time]

if(lord>1){
resagr <- aggregate((res[,rord]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,rord]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
} else {
resagr <- aggregate((res[,1]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,1]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
names(resagr)[3]<- rord
names(ressd)[3]<- rord
}


df <- reshape2::melt(resagr,id=c('time','group'))
dfsd <- reshape2::melt(ressd,id=c('time','group'))
df$min <- df$value-dfsd$value
df$max <- df$value+dfsd$value
df$varname <- sapply(as.character(df$variable),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
df <- na.omit(df)
df <- df[df$time %in% timepoints,]

for(k in unique(df$group)[unique(df$group)!=compare.to]){
  df[,paste(k,"p",sep="_")] <- NA
  df[,paste("star",k,sep="_")] <- NA
  for(i in timepoints){
  for(j in intersect(df$variable,result[[i]]$taxon)){
  df[df$time==i&df$variable==j&df$group==k,paste(k,"p",sep="_")] <- try(result[[i]][result[[i]]$taxon==j,paste(k,"_p",sep="")])
  }}
  df[df[,paste(k,"p",sep="_")]<0.05&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "*"
  df[df[,paste(k,"p",sep="_")]<0.01&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "**"
  df[df[,paste(k,"p",sep="_")]<0.001&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "***"
  }
p <- ggplot(df, aes(y=value,x=time,group=as.factor(group),color=as.factor(group)))+
  facet_wrap(~varname,ncol=ncols,scales='free')+
  geom_line(position = position_dodge(width = 0.2))+
 geom_point(position = position_dodge(width = 0.2))+
  scale_color_manual(values=colorlist,name="Group")+
#scale_y_continuous(labels = trans_format("exp", format = scientific_format(digits=1)))+
geom_linerange(position = position_dodge(width = 0.2),
                 aes(x=time,group=as.factor(group),color=as.factor(group),ymin=min,ymax=max),lwd=0.2)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.background = element_rect(colour="white", fill="white"))+
 xlab(time)+
  ylab(yl)

for(k in unique(df$group)[unique(df$group)!=compare.to]){
p <- p + geom_text(aes_(label=(df[,paste("star",k,sep="_")])), 
                   colour=colorlist[grouplist==k], vjust=0, 
                   position = position_dodge(width = 0.2),size=5)
}

if(pdf){ pdf(width=ncols*4,height=nrows*3,paste(group,compare.to,"TimeSeries_order.pdf",sep="_")); plot(p); dev.off()} 
quartz(width=ncols*4,height=nrows*3); plot(p)
}
}

#plot class level
if(length(class.table)>0){
rcla <- sigall[siglevel==2]

lcla <- length(rcla)

if(lcla > 0) {
  
if(lcla < 5) {
  ncols <- lcla
  nrows <- 1
} else {
ncols <-  round(sqrt(lcla)) + 1
nrows <-  floor(sqrt(lcla))
}

cla <- read.delim(class.table)

#colnames(cla)[grepl(pattern="Erysipelotrichi_",x=colnames(cla))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(cla)[grepl(pattern="Erysipelotrichi_",x=colnames(cla))])

if(relative) cla <- (cla+1)/rowSums(cla)
res <- as.data.frame(cla[,rcla])
res$group <-cm$group
if(log.time) res$time <- log(cm[,time]+1) else res$time <- cm[,time]

if(lcla>1){
resagr <- aggregate((res[,rcla]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,rcla]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
} else {
resagr <- aggregate((res[,1]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,1]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
names(resagr)[3]<- rcla
names(ressd)[3]<- rcla
}


df <- reshape2::melt(resagr,id=c('time','group'))
dfsd <- reshape2::melt(ressd,id=c('time','group'))
df$min <- df$value-dfsd$value
df$max <- df$value+dfsd$value
df$varname <- sapply(as.character(df$variable),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
df <- na.omit(df)
df <- df[df$time %in% timepoints,]

for(k in unique(df$group)[unique(df$group)!=compare.to]){
  df[,paste(k,"p",sep="_")] <- NA
  df[,paste("star",k,sep="_")] <- NA
  for(i in timepoints){
  for(j in intersect(df$variable,result[[i]]$taxon)){
  df[df$time==i&df$variable==j&df$group==k,paste(k,"p",sep="_")] <- try(result[[i]][result[[i]]$taxon==j,paste(k,"_p",sep="")])
  }}
  df[df[,paste(k,"p",sep="_")]<0.05&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "*"
  df[df[,paste(k,"p",sep="_")]<0.01&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "**"
  df[df[,paste(k,"p",sep="_")]<0.001&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "***"
  }
  
p <- ggplot(df, aes(y=value,x=time,group=as.factor(group),color=as.factor(group)))+
  facet_wrap(~varname,ncol=ncols,scales='free')+
  geom_line(position = position_dodge(width = 0.2))+
 geom_point(position = position_dodge(width = 0.2))+
  scale_color_manual(values=colorlist,name="Group")+
#scale_y_continuous(labels = trans_format("exp", format = scientific_format(digits=1)))+
geom_linerange(position = position_dodge(width = 0.2),
                 aes(x=time,group=as.factor(group),color=as.factor(group),ymin=min,ymax=max),lwd=0.2)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.background = element_rect(colour="white", fill="white"))+
 xlab(time)+
  ylab(yl)

for(k in unique(df$group)[unique(df$group)!=compare.to]){
p <- p + geom_text(aes_(label=(df[,paste("star",k,sep="_")])), 
                   colour=colorlist[grouplist==k], vjust=0, 
                   position = position_dodge(width = 0.2),size=5)
}

if(pdf) {pdf(width=ncols*4,height=nrows*3,paste(group,compare.to,"TimeSeries_class.pdf",sep="_")); plot(p); dev.off()} 
quartz(width=ncols*4,height=nrows*3); plot(p)
}
}


#plot phylum level
if(length(phylum.table)>0){
rphy <- sigall[siglevel==1]

lphy <- length(rphy)

if(lphy > 0) {
  
if(lphy < 5) {
  ncols <- lphy
  nrows <- 1
} else {
ncols <-  round(sqrt(lphy)) + 1
nrows <-  floor(sqrt(lphy))
}


phy <- read.delim(phylum.table)

if(relative) phy <- (phy+1)/rowSums(phy)
res <- as.data.frame(phy[,rphy])
res$group <-cm$group
if(log.time) res$time <- log(cm[,time]+1) else res$time <- cm[,time]

if(lphy>1){
resagr <- aggregate((res[,rphy]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,rphy]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
} else {
resagr <- aggregate((res[,1]),by=list(time=res$time,group=res$group),mean,na.rm=T)
ressd <- aggregate((res[,1]),by=list(time=res$time,group=res$group),FUN=function(x) sd(x)/sqrt(length(x)))
names(resagr)[3]<- rphy
names(ressd)[3]<- rphy
}

df <- reshape2::melt(resagr,id=c('time','group'))
dfsd <- reshape2::melt(ressd,id=c('time','group'))
df$min <- df$value-dfsd$value
df$max <- df$value+dfsd$value
df$varname <- sapply(as.character(df$variable),function(x) strsplit(x,split="_")[[1]][length(strsplit(x,split="_")[[1]])])
df <- na.omit(df)
df <- df[df$time %in% timepoints,]

for(k in unique(df$group)[unique(df$group)!=compare.to]){
  df[,paste(k,"p",sep="_")] <- NA
  df[,paste("star",k,sep="_")] <- NA
  for(i in timepoints){
  for(j in intersect(df$variable,result[[i]]$taxon)){
  df[df$time==i&df$variable==j&df$group==k,paste(k,"p",sep="_")] <- try(result[[i]][result[[i]]$taxon==j,paste(k,"_p",sep="")])
  }}
  df[df[,paste(k,"p",sep="_")]<0.05&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "*"
  df[df[,paste(k,"p",sep="_")]<0.01&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "**"
  df[df[,paste(k,"p",sep="_")]<0.001&!is.na(df[,paste(k,"p",sep="_")]),paste("star",k,sep="_")] <- "***"
  }

  
p <- ggplot(df, aes(y=value,x=time,group=as.factor(group),color=as.factor(group)))+
  facet_wrap(~varname,ncol=ncols,scales='free')+
  geom_line(position = position_dodge(width = 0.2))+
 geom_point(position = position_dodge(width = 0.2))+
  scale_color_manual(values=colorlist,name="Group")+
#scale_y_continuous(labels = trans_format("exp", format = scientific_format(digits=3)))+
geom_linerange(position = position_dodge(width = 0.2),
                 aes(x=time,group=as.factor(group),color=as.factor(group),ymin=min,ymax=max),lwd=0.2)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.background = element_rect(colour="white", fill="white"))+
 xlab(time)+
  ylab(yl)

for(k in unique(df$group)[unique(df$group)!=compare.to]){
p <- p + geom_text(aes_(label=(df[,paste("star",k,sep="_")])), 
                   colour=colorlist[grouplist==k], vjust=0, 
                   position = position_dodge(width = 0.2),size=5)
}


if(pdf) {
pdf(width=ncols*4,height=nrows*3,paste(group,compare.to,"TimeSeries_phylum.pdf",sep="_"));  plot(p); dev.off()
}
    quartz(width=ncols*4,height=nrows*3); plot(p)
}
}
}

