GroupPlot <- function(taxa, group = NULL, taxonomic.table, meta, readcount.cutoff = 0, 
         stacked = T, bar = T, box = T, bean = T, covariate = NULL, 
         smooth.method ='loess',select.by = NULL, select = NULL, pdf = F, quartz = T, label.direction = 1){

  taxatable <- read.delim(taxonomic.table)
  taxatable <- taxatable[,taxa]
  names(taxatable)<-sapply(names(taxatable),function(x) gsub('_NA', '.',x))
  names(taxatable) <- sapply(names(taxatable),function(x) strsplit(x,split='_')[[1]][length(strsplit(x,split='_')[[1]])])
  names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])] <- paste(names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])],c(1:length(names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])])),sep="")
  taxa <- names(taxatable)

  meta <- read.delim(meta)
  dataset <- data.frame(meta,(taxatable/meta$ReadCount)*100)
  dataset <- dataset[dataset$ReadCount>readcount.cutoff,]
  if (length(select.by)!=0){
    dataset$selection <- dataset[,select.by]
    dataset <- dataset[dataset$selection==select,]
  } 
 dataset[,group]<-dataset[,group][drop=T]

if(length(covariate)!=0){
if(group=='') {
dataset[,'group']<-1
group <- 'group'
legpos <- 'none'}
else(legpos <- 'right')

df = na.omit(reshape2::melt(dataset[,c(covariate,group,taxa)], id=c(covariate,group)))
names(df) <- c("x","gr","variable","value") 
p<- ggplot2::ggplot(df, ggplot2::aes(y=value, x=x,color=factor(gr)),environment = environment()) +
  ggplot2::stat_smooth(method = smooth.method, formula = y ~ x, ggplot2::aes(fill=factor(gr)),se=T) +
  ggplot2::facet_wrap(~variable,ncol=floor(sqrt(length(taxa))),scales='free')+
  ggplot2::theme_bw()+
  ggplot2::xlab(covariate)+
  ggplot2::ylab('Relative abundance (%)')+
  ggplot2::theme(legend.position=legpos)+
  ggplot2::scale_color_manual(name=group,values=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])

if(pdf){
pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_",group, "_", select.by,select, "_", "Covariateplot.pdf", sep = ""));
plot(p)
dev.off() 
}

if (quartz) quartz()
plot(p)
}

if (box){
if(pdf){
pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_", select.by,select, "Boxplot.pdf", sep = ""))
par(mfrow=c(floor(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.2,0),mar=c(5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
for(i in taxa) {
boxplot(dataset[,i]~dataset[,group],ylab='',xlab="",outpch=21,axes=F, las=label.direction,
col=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray')[1:length(unique(dataset[,group]))],
outbg=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray')[1:length(unique(dataset[,group]))],
outcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
boxcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
medcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
whiskcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
staplecol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
mtext(side=1,line=0,at=c(1:length(levels(as.factor(dataset[,group])))),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
axis(side=2)
mtext(side=2,text=i,line=1.5,cex=1,font=3)
 }
dev.off() }
if (quartz) quartz()
  else x11()
par(mfrow=c(floor(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.2,0),mar=c(5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
for(i in taxa) {
boxplot(dataset[,i]~dataset[,group],ylab="",xlab="",outpch=21,axes=F, las=label.direction,
col=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray')[1:length(unique(dataset[,group]))],
outbg=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray')[1:length(unique(dataset[,group]))],
outcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
boxcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
medcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
whiskcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
staplecol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
mtext(side=1,line=0,at=c(1:length(levels(as.factor(dataset[,group])))),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
axis(side=2)
mtext(side=2,text=i,line=1.5,cex=1,font=3)}
}

if(bar){
if(pdf){
pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_", select.by,select, "Barplot.pdf", sep = ""));par(mfrow=c(round(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.2,0),mar=c(3.5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
  palette(c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray','royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
trmeans<-data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1]))
trse<- data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),FUN=sd,na.rm=T)[,-1])/c(sqrt(table(as.factor(dataset[,group])))))
names(trmeans)<-levels(as.factor(dataset[,group]))
names(trse)<-levels(as.factor(dataset[,group]))    
for(i in rownames(trmeans)){ 
plot(unlist(trmeans[i,]),type="h",lwd=15,lend=1,axes=F,ylab="",xlab="",
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]),na.rm=T)),
xlim=c(0.5,ncol(trmeans)+0.5),
col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
lines(unlist(trmeans[i,]-(0.001*max(unlist(1.05*trmeans[i,])))),type="h",lwd=14,lend=1,
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]))),
xlim=c(0.5,ncol(trmeans)+0.5),
col=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray'))
arrows(x0=c(c(1:ncol(trmeans))),x1=c(c(1:ncol(trmeans))),
       y0=unlist(trmeans[i,])-unlist(trse[i,]),
       y1=unlist(trmeans[i,])+unlist(trse[i,]),
       length=0.05,angle=90,code=3,
       col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
mtext(side=1,line=0.5,at=c(1:ncol(trmeans)),text=levels(as.factor(dataset[,group])),cex=0.7,las=label.direction)
mtext(side=2,text=i,line=1.5,cex=1,font=3)
axis(side=2,at=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10)),
     labels=pretty(seq(0,max(unlist(trmeans[i,])),
                      max(unlist(trmeans[i,]))/10),digits=2))}
dev.off()}  
if (quartz) quartz()
  else x11()
par(mfrow=c(floor(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.2,0),mar=c(5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
  palette(c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray','royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
trmeans<-data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1]))
trse<- data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),FUN=sd,na.rm=T)[,-1])/c(sqrt(table(as.factor(dataset[,group])))))
names(trmeans)<-levels(as.factor(dataset[,group]))
names(trse)<-levels(as.factor(dataset[,group]))    
for(i in rownames(trmeans)){ 
plot(unlist(trmeans[i,]),type="h",lwd=15,lend=1,axes=F,ylab="",xlab="",
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]),na.rm=T)),
xlim=c(0.5,ncol(trmeans)+0.5),#col=c(1:ncol(trmeans)),
col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
lines(unlist(trmeans[i,]-(0.001*max(unlist(1.05*trmeans[i,])))),type="h",lwd=14,lend=1,
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]))),
xlim=c(0.5,ncol(trmeans)+0.5),
col=c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray'))
arrows(x0=c(c(1:ncol(trmeans))),x1=c(c(1:ncol(trmeans))),
       y0=unlist(trmeans[i,])-unlist(trse[i,]),
       y1=unlist(trmeans[i,])+unlist(trse[i,]),
       length=0.05,angle=90,code=3,
       col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
mtext(side=1,line=0.5,at=c(1:ncol(trmeans)),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
mtext(side=2,text=i,line=1.5,cex=1,font=3)
axis(side=2,at=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10)),
     labels=pretty(seq(0,max(unlist(trmeans[i,])),
                      max(unlist(trmeans[i,]))/10),digits=2))}
}

if(stacked){
if(pdf){
pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_", select.by,select, "Stackedplot.pdf", sep = ""))
op<-par(xpd=T,mar=c(10,20,2,2),cex.axis=1.5,cex.lab=1.5)
barplot(as.matrix(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1])),
col=(c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')),
legend=T,args.legend=list(x=-4,text.font=3,bty="n"),
names.arg=levels(as.factor(dataset[,group])),las=label.direction,xlab="",
ylab="% of total microbiota", ylim=c(0,100))
dev.off()
par(op)}
if (quartz) quartz()
  else x11()
op<-par(xpd=T,mar=c(10,20,2,2),cex.axis=1.5,cex.lab=1.5)
barplot(as.matrix(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1])),
col=(c('skyblue','yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')),
legend=T,args.legend=list(x=-4,text.font=3,bty="n"),
names.arg=levels(as.factor(dataset[,group])),las=label.direction,xlab="",
ylab="% of total microbiota", ylim=c(0,100))
par(op)
}

if(bean){
if (quartz) quartz()
  else x11()
  par(mfrow=c(floor(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.5,0),mar=c(5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
for(i in taxa) {
beanplot::beanplot(dataset[,i]~dataset[,group],ll=0.1,ylab=i,las=label.direction,xlab="",
col=list(c('skyblue','royalblue','royalblue','royalblue'),
         c('yellowgreen','olivedrab4','olivedrab4','olivedrab4'),
          c('pink','red','red','red'),
          c('turquoise2','turquoise4','turquoise4','turquoise4'),
          c('plum','purple','purple','purple'),
          c('darkorange','darkorange3','darkorange3','darkorange3'),
          c('lightyellow','lightyellow4','lightyellow4','lightyellow4'),
          c('gray','black','black','black')),
border=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))}

if(pdf){
pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_", select.by,select, "Beanplot.pdf", sep = ""))
  par(mfrow=c(floor(sqrt(length(taxa))),round(sqrt(length(taxa)))+1),mgp=c(2,0.2,0),mar=c(5,3.5,1,1),tck=-0.01,cex.axis=1.5,cex.lab=1.5)
for(i in taxa) {
beanplot::beanplot(dataset[,i]~dataset[,group],ll=0.1,ylab=i,las=label.direction,xlab="",
col=list(c('skyblue','royalblue','royalblue','royalblue'),
         c('yellowgreen','olivedrab4','olivedrab4','olivedrab4'),
          c('pink','red','red','red'),
          c('turquoise2','turquoise4','turquoise4','turquoise4'),
          c('plum','purple','purple','purple'),
          c('darkorange','darkorange3','darkorange3','darkorange3'),
          c('lightyellow','lightyellow4','lightyellow4','lightyellow4'),
          c('gray','black','black','black')),
border=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))}
 dev.off()
par(op)} 
}
 palette("default")
}
