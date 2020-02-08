Clusters <- function(taxonomic.table, meta, N.taxa = NULL,  readcount.cutoff = 0,
                     minimum.correlation = 0.5, minimum.network = 1,
                      select.by = NULL, select = NULL, keep.result = F, pdf = F, relative = T){
if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}

cluster.similarity = 1-minimum.correlation 

taxatable <- read.delim(taxonomic.table)
metadata <- read.delim(meta)
if(relative) taxatable <- taxatable/metadata$ReadCount
taxatable <- taxatable[metadata$ReadCount > readcount.cutoff, ]

if (length(select.by) != 0) {
        metadata$selection <- metadata[, select.by]
        taxatable <- taxatable[metadata$selection == select, ]
        metadata <- metadata[metadata$selection == select, ]
    }

if (length(N.taxa) == 0) N.taxa = ncol(taxatable)
vars <- c(rev(names(colSums(taxatable,na.rm=T)[order(colSums(taxatable,na.rm=T))])))[1:N.taxa]
gs<-taxatable[,vars]
tgs <-data.frame(t(scale(gs)))
g2.cor<-cor(t(tgs),method="spearman",use="pairwise.complete.obs")
g2.cor[is.na(g2.cor)] <- 0

g2.cor2 <- g2.cor
g2.cor2[abs(g2.cor2)<minimum.correlation] <- 0 
g2.cor2 <- g2.cor2[rowSums(abs(g2.cor2)>0)>minimum.network,rowSums(abs(g2.cor2)>0)>minimum.network]

spnames1 <- rownames(g2.cor)
spnames1 <- sapply(spnames1, function(x) gsub("_NA", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_1", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_2", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_3", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_4", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_5", ".", x))
spnames1 <- sapply(spnames1, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

spnames <- rownames(g2.cor2)
classnames <- sapply(spnames, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
spnames <- sapply(spnames, function(x) gsub("_NA", ".", x))
spnames <- sapply(spnames, function(x) gsub("_1", ".", x))
spnames <- sapply(spnames, function(x) gsub("_2", ".", x))
spnames <- sapply(spnames, function(x) gsub("_3", ".", x))
spnames <- sapply(spnames, function(x) gsub("_4", ".", x))
spnames <- sapply(spnames, function(x) gsub("_5", ".", x))
spnames <- sapply(spnames, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

clusters <- hclust(as.dist(1-g2.cor),"average")
clus<-cutree(clusters,h=cluster.similarity)


if (pdf){
pdf(paste("CorrelatingTaxa_",select.by,select,".pdf",sep=""))
plot(clusters, ylab="",labels=spnames1,xlab="",cex=0.5) 
abline(h=cluster.similarity, lty=2, col="gray")

qgraph::qgraph(g2.cor2,vsize=5,rescale=T,repulsion=0.8,
               labels=substr(spnames,start=1,stop=4),
               layout="spring",diag=F,
       legend.cex=0.5,
       groups=classnames,
       color=c("#E41A1C","#FFA500","#377EB8","#87CEFA","#4DAF4A" ,'#9ACD32',"#984EA3",'#DA70D6', "#999999","gainsboro",
      "#008080","#00CED1","#F781BF","thistle1","#8DA0CB","lightsteelblue1","#FFD92F","#FFFFB3",
"#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
"#CCEBC5","#FFED6F","#C71585","#EE82EE","#66C2A5","#FC8D62","#A65628")[1:length(unique(classnames))],
       label.prop=0.99)
mtext(side=3,text="Correlations",line=2)

dev.off()  
} 

quartz()
plot(clusters, ylab="",labels=spnames1,xlab="",cex=0.5) 
abline(h=cluster.similarity, lty=2, col="gray")

quartz()
qgraph::qgraph(g2.cor2,vsize=5,rescale=T,repulsion=0.8,
               labels=substr(spnames,start=1,stop=4),
               layout="spring",diag=F,
       legend.cex=0.5,
       groups=classnames,
       color=c("#E41A1C","#FFA500","#377EB8","#87CEFA","#4DAF4A" ,'#9ACD32',"#984EA3",'#DA70D6', "#999999","gainsboro",
      "#008080","#00CED1","#F781BF","thistle1","#8DA0CB","lightsteelblue1","#FFD92F","#FFFFB3",
"#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
"#CCEBC5","#FFED6F","#C71585","#EE82EE","#66C2A5","#FC8D62","#A65628")[1:length(unique(classnames))],
       label.prop=0.99)

networks <- data.frame(metadata,taxatable)
for(i in names(table(clus)[table(clus)>1])) networks[,paste('cluster',i,sep="")] <- rowSums(networks[, names(clus)[clus==i]],na.rm=T)
for(i in names(table(clus)[table(clus)==1])) networks[,paste('cluster',i,sep="")] <- networks[, names(clus)[clus==i]]
networks <- list(networks, clus)
write.table(networks[[1]], file = "Clusters.txt", quote=F, sep="\t")
if (keep.result) return(networks)
}

