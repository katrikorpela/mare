GroupPlot <- function(taxa = NULL, group = NULL, taxonomic.table, meta, readcount.cutoff = 0, 
         stacked = T, bar = T, box = T, bean = T, covariate = NULL, 
         smooth.method ='loess',select.by = NULL, select = NULL, pdf = F, label.direction = 1, logtrans = F){

  if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}

  taxatable <- read.delim(taxonomic.table)
  if(length(taxa)>0) taxatable <- taxatable[,taxa]
  
  colnames(taxatable)[grepl(pattern="incertae_sedis",x=colnames(taxatable))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxatable)[grepl(pattern="incertae_sedis",x=colnames(taxatable))])
  colnames(taxatable)[grepl(pattern="Incertae_Sedis_",x=colnames(taxatable))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(taxatable)[grepl(pattern="Incertae_Sedis_",x=colnames(taxatable))])
  colnames(taxatable)[grepl(pattern="Erysipelotrichi_",x=colnames(taxatable))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(taxatable)[grepl(pattern="Erysipelotrichi_",x=colnames(taxatable))])
  
  names(taxatable)<- sapply(names(taxatable),function(x) gsub('_NA', '.',x))
  names(taxatable) <- sapply(names(taxatable),function(x) strsplit(x,split='_')[[1]][length(strsplit(x,split='_')[[1]])])
  #names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])] <- paste(names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])],c(1:length(names(taxatable)[names(taxatable)==names(table(names(taxatable))[table(names(taxatable))>1])])),sep="")
  for(i in names(table(names(taxatable))[table(names(taxatable))>1])) names(taxatable)[grep(names(taxatable),pattern=i)] <- paste(i,c(1:length(names(taxatable)[grep(names(taxatable),pattern=i)])),sep="")
    
   taxa <- names(taxatable)

  meta <- read.delim(meta)
  if(logtrans) {
    dataset <- data.frame(meta,((taxatable+1)/meta$ReadCount)*100)
    } else  dataset <- data.frame(meta,((taxatable)/meta$ReadCount)*100)
  dataset <- dataset[dataset$ReadCount>readcount.cutoff,]
  if (length(select.by)!=0){
    dataset$selection <- dataset[,select.by]
    dataset <- dataset[dataset$selection==select,]
  } 
 dataset <-dataset[!is.na(dataset[,group]),] 
 dataset[,group]<-as.factor(dataset[,group][drop=T])

 if(length(taxa)>10) {legsize = 1.5  } else legsize <- 2 
 
 
# picwidth <- log(length(levels(as.factor(dataset[,group]))))*10
#if(picwidth>10) picwidth2 = 10 else picwidth2 = picwidth
  
 palette(c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))

 
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
  ggplot2::facet_wrap(~variable,ncol=round(sqrt(length(taxa))),scales='free')+
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

quartz()
plot(p)
}
 
if(stacked){

if(pdf){
pdf(width=picwidth,paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Stackedplot.pdf", sep = ""))
op<-par(xpd=T,mar=c(10,5,2,20),cex.axis=1.5,cex.lab=1.5)
barplot(as.matrix(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1])),
col=(c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')),
#legend=T,args.legend=list(x=-(picwidth/4),text.font=3,bty="n"),
legend=T,args.legend=list(x="left",inset=1,text.font=3,bty="n"),
names.arg=levels(as.factor(dataset[,group])),las=label.direction,xlab="",
#yaxt="n",
ylab="% of total microbiota")#, ylim=c(0,100))
#axis(side=2, line=-5)
#mtext(side=2,line=-2,text="% of total microbiota",cex=1.5)
dev.off()
par(op)}

#quartz(width=picwidth2)
  quartz()
op<-par(xpd=T,mar=c(10,5,2,20),cex.axis=1.5,cex.lab=1.5)
barplot(as.matrix(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1])),
col=(c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')),
#legend=T,args.legend=list(x=-(picwidth/4),text.font=3,bty="n"),
legend=T,args.legend=list(x="left",inset=1,text.font=3,bty="n"),
names.arg=levels(as.factor(dataset[,group])),las=label.direction,xlab="",
#yaxt="n",
ylab="% of total microbiota")#, ylim=c(0,100))
#axis(side=2, line=-5)
#mtext(side=2,line=-2,text="% of total microbiota",cex=1.5)
par(op)

} 
 
 
if(bean){
if(logtrans) { lg="y" } else lg=""
  
#quartz(width=picwidth2)
   quartz()
   par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.5,0),mar=c(1,3.5,1,1),tck=-0.01,cex.axis=1,cex.lab=1)
  for(i in taxa) {
try(beanplot::beanplot(dataset[,i]~dataset[,group],ll=0.1,ylab=i,
                       xlab="",xaxt="n",#las=label.direction,
col=list(c('#E41A1C','black','black','black'),
         c('orange','black','black','black'),
          c('#377EB8','black','black','black'),
          c('skyblue','black','black','black'),
          c("#4DAF4A",'black','black','black'),
         c("darkolivegreen2",'black','black','black'),
          c('#984EA3','black','black','black'),
          c('#FFFF33','black','black','black'),
         c('#A65628','black','black','black'),
         c('#F781BF','black','black','black'),
         c('#999999','black','black','black'),
          c('blue','black','black','black'),
          c('firebrick4','black','black','black'),
          c('yellowgreen','black','black','black'),
          c('pink','black','black','black'),
          c('turquoise2','black','black','black'),
          c('plum','black','black','black'),
          c('darkorange','black','black','black'),
          c('lightyellow','black','black','black'),
          c('gray','black','black','black')),
log=lg,
border=c('black')))
}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)

df = na.omit(reshape2::melt(dataset[,c(group,taxa)], id=c(group)))
names(df)[1] <- c("group")

#quartz(width=picwidth2)
  quartz()
  par(mar=c(20,4,3,3))
tryCatch(beanplot::beanplot((value)~paste(variable,group),data=df,las=2,
        yaxt="n",xaxt="n",
        ylab="Relative abundance (%)",
        col=list(c('#E41A1C','black','black','black'),
         c('orange','black','black','black'),
          c('#377EB8','black','black','black'),
          c('skyblue','black','black','black'),
          c("#4DAF4A",'black','black','black'),
         c("darkolivegreen2",'black','black','black'),
          c('#984EA3','black','black','black'),
          c('#FFFF33','black','black','black'),
         c('#A65628','black','black','black'),
         c('#F781BF','black','black','black'),
         c('#999999','black','black','black'),
          c('blue','black','black','black'),
          c('firebrick4','black','black','black'),
          c('yellowgreen','black','black','black'),
          c('pink','black','black','black'),
          c('turquoise2','black','black','black'),
          c('plum','black','black','black'),
          c('darkorange','black','black','black'),
          c('lightyellow','black','black','black'),
          c('gray','black','black','black'))[1:length(levels(dataset[,group]))],
border=c('black'),
log=lg,
at=rep(rank(colMeans(dataset[,taxa][order(taxa)]),ties.method = "first")*length(levels(dataset[,group])),
       each=length(levels(dataset[,group])))-c(sequence(length(levels(dataset[,group])))/2)-length(levels(dataset[,group]))/4),
error = function(e) NULL)
if(logtrans){ tryCatch(axis(side=2,at=c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100),
     labels=c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100),las=2),error = function(e) NULL)
} else { tryCatch(axis(side=2,at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100),las=2),error = function(e) NULL)
     }
tryCatch(axis(side=1,
     at=seq(length(levels(dataset[,group]))/2,length(unique(paste(df$variable,df$group)))-length(levels(dataset[,group]))/2,length(levels(dataset[,group]))),
     labels=taxa[order(colMeans(dataset[,taxa]))],las=label.direction),error = function(e) NULL)
tryCatch(legend("topleft",title=group,legend=levels(dataset[,group]),fil=c(1:length(levels(dataset[,group]))),bty="n"),error = function(e) NULL)


if(pdf){
pdf(width=picwidth, paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Beanplots.pdf", sep = ""))
  par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.5,0),mar=c(1,3.5,1,1),tck=-0.01,cex.axis=1,cex.lab=1)
  for(i in taxa) {
  try(beanplot::beanplot(dataset[,i]~dataset[,group],ll=0.1,ylab=i,xlab="",xaxt="n",#las=label.direction,
col=list(c('#E41A1C','black','black','black'),
         c('orange','black','black','black'),
          c('#377EB8','black','black','black'),
          c('skyblue','black','black','black'),
          c("#4DAF4A",'black','black','black'),
         c("darkolivegreen2",'black','black','black'),
          c('#984EA3','black','black','black'),
          c('#FFFF33','black','black','black'),
         c('#A65628','black','black','black'),
         c('#F781BF','black','black','black'),
         c('#999999','black','black','black'),
          c('blue','black','black','black'),
          c('firebrick4','black','black','black'),
          c('yellowgreen','black','black','black'),
          c('pink','black','black','black'),
          c('turquoise2','black','black','black'),
          c('plum','black','black','black'),
          c('darkorange','black','black','black'),
          c('lightyellow','black','black','black'),
          c('gray','black','black','black')),
log=lg,
border=c('black')))
}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)

 dev.off()

 #-------
df = na.omit(reshape2::melt(dataset[,c(group,taxa)], id=c(group)))
names(df)[1] <- c("group")
pdf(width=picwidth,paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Beanplot.pdf", sep = ""))
par(mar=c(20,4,3,3))
tryCatch(beanplot::beanplot((value)~paste(variable,group),data=df,las=2,
        yaxt="n",xaxt="n",
        ylab="Relative abundance (%)",
        col=list(c('#E41A1C','black','black','black'),
         c('orange','black','black','black'),
          c('#377EB8','black','black','black'),
          c('skyblue','black','black','black'),
          c("#4DAF4A",'black','black','black'),
         c("darkolivegreen2",'black','black','black'),
          c('#984EA3','black','black','black'),
          c('#FFFF33','black','black','black'),
         c('#A65628','black','black','black'),
         c('#F781BF','black','black','black'),
         c('#999999','black','black','black'),
          c('blue','black','black','black'),
          c('firebrick4','black','black','black'),
          c('yellowgreen','black','black','black'),
          c('pink','black','black','black'),
          c('turquoise2','black','black','black'),
          c('plum','black','black','black'),
          c('darkorange','black','black','black'),
          c('lightyellow','black','black','black'),
          c('gray','black','black','black'))[1:length(levels(dataset[,group]))],
border=c('black'),
log=lg,
at=rep(rank(colMeans(dataset[,taxa][order(taxa)]),ties.method = "first")*length(levels(dataset[,group])),
       each=length(levels(dataset[,group])))-c(sequence(length(levels(dataset[,group])))/2)-length(levels(dataset[,group]))/4),
error = function(e) NULL)
if(logtrans){ tryCatch(axis(side=2,at=c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100),
     labels=c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100),las=2),error = function(e) NULL)
} else { tryCatch(axis(side=2,at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100),las=2),error = function(e) NULL)
     }
tryCatch(axis(side=1,
     at=seq(length(levels(dataset[,group]))/2,length(unique(paste(df$variable,df$group)))-length(levels(dataset[,group]))/2,length(levels(dataset[,group]))),
     labels=taxa[order(colMeans(dataset[,taxa]))],las=label.direction),error = function(e) NULL)
tryCatch(legend("topleft",title=group,legend=levels(dataset[,group]),fil=c(1:length(levels(dataset[,group]))),bty="n"),error = function(e) NULL)
dev.off()
#-------

 
 par(op)

} 
} 
 
if(bar){
if(pdf){
pdf(width=picwidth*0.75, paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Barplot.pdf", sep = ""))
par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.2,0),mar=c(1,3.5,1,1),tck=-0.01,cex.axis=1.5)
trmeans<-data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1]))
trse<- data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),FUN=sd,na.rm=T)[,-1])/c(sqrt(table(as.factor(dataset[,group])))))
names(trmeans)<-levels(as.factor(dataset[,group]))
names(trse)<-levels(as.factor(dataset[,group]))    
for(i in rownames(trmeans)){ 
plot(unlist(trmeans[i,]),type="h",lwd=15,lend=1,axes=F,ylab="",xlab="",
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]),na.rm=T)),
xlim=c(0.5,ncol(trmeans)+0.5),#col=c(1:ncol(trmeans)),
#col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
col="black")
lines(unlist(trmeans[i,]-(0.001*max(unlist(1.05*trmeans[i,])))),type="h",lwd=14,lend=1,
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]))),
xlim=c(0.5,ncol(trmeans)+0.5),
col=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
arrows(x0=c(c(1:ncol(trmeans))),x1=c(c(1:ncol(trmeans))),
       y0=unlist(trmeans[i,])-unlist(trse[i,]),
       y1=unlist(trmeans[i,])+unlist(trse[i,]),
       length=0.05,angle=90,code=3,
       col="black")#c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
#mtext(side=1,line=0.5,at=c(1:ncol(trmeans)),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
mtext(side=2,text=i,line=1.5,cex=0.75,font=3)
axis(side=2,at=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10)),
  labels=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10),digits=2))
}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)

dev.off()
}  
#quartz(width=picwidth2*0.75) 
  quartz()
  par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.2,0),mar=c(1,3.5,1,1),tck=-0.01,cex.axis=1.5)
trmeans<-data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),mean,na.rm=T)[,-1]))
trse<- data.frame(t(aggregate(dataset[,taxa],by=list(group=dataset[,group]),FUN=sd,na.rm=T)[,-1])/c(sqrt(table(as.factor(dataset[,group])))))
names(trmeans)<-levels(as.factor(dataset[,group]))
names(trse)<-levels(as.factor(dataset[,group]))    
for(i in rownames(trmeans)){ 
plot(unlist(trmeans[i,]),type="h",lwd=15,lend=1,axes=F,ylab="",xlab="",
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]),na.rm=T)),
xlim=c(0.5,ncol(trmeans)+0.5),#col=c(1:ncol(trmeans)),
#col=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
col="black")
lines(unlist(trmeans[i,]-(0.001*max(unlist(1.05*trmeans[i,])))),type="h",lwd=14,lend=1,
ylim=c(0,max(unlist(1.05*trmeans[i,])+unlist(trse[i,]))),
xlim=c(0.5,ncol(trmeans)+0.5),
col=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
arrows(x0=c(c(1:ncol(trmeans))),x1=c(c(1:ncol(trmeans))),
       y0=unlist(trmeans[i,])-unlist(trse[i,]),
       y1=unlist(trmeans[i,])+unlist(trse[i,]),
       length=0.05,angle=90,code=3,
       col="black")#c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black'))
#mtext(side=1,line=0.5,at=c(1:ncol(trmeans)),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
mtext(side=2,text=i,line=1.5,cex=0.75,font=3)
axis(side=2,at=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10)),
  labels=pretty(seq(0,max(unlist(trmeans[i,])),max(unlist(trmeans[i,]))/10),digits=2))
}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)

} 
 
 
if(logtrans) dataset[,taxa] <- log(dataset[,taxa]+0.00001)
 
if (box){
if(pdf){
pdf(width=picwidth, paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Boxplots.pdf", sep = ""))
par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.2,0),mar=c(1,3.5,1,1),tck=-0.01,
    cex.axis=1.5,cex.lab=0.5)
for(i in taxa) {
boxplot(dataset[,i]~dataset[,group],ylab="",xlab="",outpch=21,axes=F, #las=label.direction,
col=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
outbg=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
#outcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
##boxcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
#medcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
#whiskcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
#staplecol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
#mtext(side=1,line=0,at=c(1:length(levels(as.factor(dataset[,group])))),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
mtext(side=2,text=i,line=1.5,cex=0.75,font=3)
if(logtrans) { axis(side=2, at=log(c(0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01, 0.05,0.1,0.5,1,5,10,50,100)), 
                    labels=c(0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01, 0.05,0.1,0.5,1,5,10,50,100))
} else axis(side=2)
}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)

dev.off() 

#-------
df = na.omit(reshape2::melt(dataset[,c(group,taxa)], id=c(group)))
names(df)[1] <- c("group")
pdf(width=picwidth, paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",group,"_", select.by,select, "Boxplot.pdf", sep = ""))
par(mar=c(20,4,3,3))
boxplot((value)~paste(variable,group),
        data=df,las=2,col=c(1:length(levels(dataset[,group]))),#log="y",
        outpch=21,outbg=c(1:length(levels(dataset[,group]))),cex=0.5,
        yaxt="n",xaxt="n",
        ylab="Relative abundance (%)",
       # boxwex=1/(length(levels(dataset[,group]))),
at=rep(rank(colMeans(dataset[,taxa][order(taxa)]),ties.method = "first")*length(levels(dataset[,group])),
       each=length(levels(dataset[,group])))-c(sequence(length(levels(dataset[,group])))/2)-length(levels(dataset[,group]))/4)
if(logtrans){ axis(side=2,at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
     labels=c(0.0001,0.001,0.01,0.1,1,10,100),las=2)
} else { axis(side=2,at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100),las=2)
     }
axis(side=1,at=seq(1.5,length(unique(paste(df$variable,df$group)))-0.5,length(levels(dataset[,group]))),
         labels=taxa[order(colMeans(dataset[,taxa]))],las=label.direction)
legend("topleft",title=group,legend=levels(dataset[,group]),fil=c(1:length(levels(dataset[,group]))),bty="n")
dev.off()
#-------
}

  
#quartz(width=picwidth2)
  quartz()
  par(mfcol=c(round(sqrt(length(taxa)+1)),round(sqrt(length(taxa)+1))+1),mgp=c(2,0.2,0),mar=c(1,3.5,1,1),tck=-0.01,
    cex.axis=1)
for(i in taxa) {
boxplot(dataset[,i]~dataset[,group],ylab="",xlab="",outpch=21,axes=F, #las=label.direction,
col=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
outbg=c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
#outcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
##boxcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
#medcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
#whiskcol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], 
#staplecol=c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))])
#mtext(side=1,line=0,at=c(1:length(levels(as.factor(dataset[,group])))),text=levels(as.factor(dataset[,group])),cex=1,las=label.direction)
if(logtrans) { axis(side=2, at=log(c(0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01, 0.05,0.1,0.5,1,5,10,50,100)), 
                    labels=c(0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01, 0.05,0.1,0.5,1,5,10,50,100))
} else axis(side=2)
mtext(side=2,text=i,line=1.5,cex=0.75,font=3)}
par(mar=c(0,0,0,0))
plot.new()
legend("topleft",bty="n", legend = levels(as.factor(dataset[,group])), 
       #border= c('royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))],
       fil = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,'darkolivegreen2',"#984EA3","#FFFF33", 
       "#A65628", "#F781BF", "#999999","dodgerblue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset[,group]))], cex=legsize)


#-------
df = na.omit(reshape2::melt(dataset[,c(group,taxa)], id=c(group)))
names(df)[1] <- c("group")
#quartz(width=picwidth2)
  quartz()
  par(mar=c(20,4,3,3))
boxplot((value)~paste(variable,group),
        data=df,las=2,col=c(1:length(levels(dataset[,group]))),#log="y",
        outpch=21,outbg=c(1:length(levels(dataset[,group]))),cex=0.5,
        yaxt="n",xaxt="n",
        ylab="Relative abundance (%)",
at=rep(rank(colMeans(dataset[,taxa][order(taxa)]),ties.method = "first")*length(levels(dataset[,group])),
       each=length(levels(dataset[,group])))-c(sequence(length(levels(dataset[,group])))/2)-length(levels(dataset[,group]))/4)
if(logtrans){ axis(side=2,at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
     labels=c(0.0001,0.001,0.01,0.1,1,10,100),las=2)
} else { axis(side=2,at=c(0,20,40,60,80,100),
     labels=c(0,20,40,60,80,100),las=2)
     }
axis(side=1,
     at=seq(length(levels(dataset[,group]))/2,length(unique(paste(df$variable,df$group)))-length(levels(dataset[,group]))/2,length(levels(dataset[,group]))),
     labels=taxa[order(colMeans(dataset[,taxa]))],las=label.direction)
legend("topleft",title=group,legend=levels(dataset[,group]),fil=c(1:length(levels(dataset[,group]))),bty="n")
#-------


}

 palette("default")
}
