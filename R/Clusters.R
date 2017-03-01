Clusters <- function(taxonomic.table, meta, N.taxa = NULL,  readcount.cutoff = 0,
                     minimum.correlation = 0.5, minimum.network = 1, cluster.similarity = 1,
                      select.by = NULL, select = NULL,
                     quartz = T, pdf = F){

taxatable <- read.delim(taxonomic.table)
metadata <- read.delim(meta)
taxatable <- taxatable/metadata$ReadCount
taxatable <- taxatable[metadata$ReadCount > readcount.cutoff, ]

if (length(select.by) != 0) {
        metadata$selection <- metadata[, select.by]
        taxatable <- taxatable[metadata$selection == select, ]
        metadata <- metadata[metadata$selection == select, ]
    }

#pt2 <- function(q,df,log.p=F) 2*pt(-abs(q),df,log.p=log.p)
if (length(N.taxa) == 0) N.taxa = ncol(taxatable)
vars <- c(rev(names(colSums(taxatable,na.rm=T)[order(colSums(taxatable,na.rm=T))])))[1:N.taxa]
gs<-taxatable[,vars]
#n<-nrow(gs)
#df<-n-2
tgs <-data.frame(t(scale(gs)))
g2.cor<-cor(t(tgs),method="spearman",use="pairwise.complete.obs")
g2.cor[is.na(g2.cor)] <- 0
#g2.tstat<-g2.cor*sqrt((n-2)/(1-g2.cor^2))#t = r*sqrt((n-2)/(1-r^2))
#g2.p<-pt2(g2.tstat,df)
#g2.adjp <- p.adjust(g2.p,method="fdr")
#dim(g2.adjp) <- c(length(vars),length(vars))
g2.cor2 <- g2.cor#;g2.cor2[g2.p>0.05] <- 0
g2.cor2[abs(g2.cor2)<minimum.correlation] <- 0 
g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]
#g2.cor2[g2.cor2>minimum.correlation]<-1
#g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]
#g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]

spnames1 <- rownames(g2.cor)
spnames1 <- sapply(spnames1, function(x) gsub("_NA", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_1", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_2", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_3", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_4", ".", x))
spnames1 <- sapply(spnames1, function(x) gsub("_5", ".", x))
classnames <- sapply(spnames1, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
spnames1 <- sapply(spnames1, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

spnames <- rownames(g2.cor2)
spnames <- sapply(spnames, function(x) gsub("_NA", ".", x))
spnames <- sapply(spnames, function(x) gsub("_1", ".", x))
spnames <- sapply(spnames, function(x) gsub("_2", ".", x))
spnames <- sapply(spnames, function(x) gsub("_3", ".", x))
spnames <- sapply(spnames, function(x) gsub("_4", ".", x))
spnames <- sapply(spnames, function(x) gsub("_5", ".", x))
classnames <- sapply(spnames, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
spnames <- sapply(spnames, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

clusters <- hclust(as.dist(1-g2.cor),"ward.D2")
clus<-cutree(clusters,h=cluster.similarity)

g2.cor3 <- g2.cor
for(i in names(table(clus)[table(clus)>1])) g2.cor3[names(clus)[clus==i],names(clus)[clus==i]]<-1
g2.cor3[g2.cor3<1]<-0
g2.cor3<-g2.cor3[vegan::specnumber(g2.cor3)>minimum.network,vegan::specnumber(g2.cor3)>minimum.network]

spnames3 <- rownames(g2.cor3)
spnames3 <- sapply(spnames3, function(x) gsub("_NA", ".", x))
spnames3 <- sapply(spnames3, function(x) gsub("_1", ".", x))
spnames3 <- sapply(spnames3, function(x) gsub("_2", ".", x))
spnames3 <- sapply(spnames3, function(x) gsub("_3", ".", x))
spnames3 <- sapply(spnames3, function(x) gsub("_4", ".", x))
spnames3 <- sapply(spnames3, function(x) gsub("_5", ".", x))
classnames <- sapply(spnames3, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
spnames3 <- sapply(spnames3, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

if (pdf){
pdf("CorrelatingTaxa.pdf") 
plot(clusters, ylab="",labels=spnames1,xlab="",cex=0.5) 
abline(h=cluster.similarity, lty=2, col="gray")

qgraph::qgraph(g2.cor2,vsize=5,rescale=T,repulsion=0.8,
               labels=substr(spnames,start=1,stop=4),
               layout="spring",diag=F,
       legend.cex=0.5,
       groups=classnames,
       color=c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray","royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black"),
       label.prop=0.99)
mtext(side=3,text="Correlations",line=2)

qgraph::qgraph(g2.cor3,vsize=5,rescale=T,repulsion=0.7,
          labels=substr(spnames3,start=1,stop=4),layout="spring",diag=F,
       legend.cex=0.5,
       groups=classnames,
       color=c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray","royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black"),
       label.prop=0.99)
mtext(side=3,text="Clusters",line=2)

dev.off()  
} 

if (quartz) quartz() else x11()
plot(clusters, ylab="",labels=spnames1,xlab="",cex=0.5) 
abline(h=cluster.similarity, lty=2, col="gray")

if (quartz) quartz() else x11()
qgraph::qgraph(g2.cor2,vsize=5,rescale=T,repulsion=0.8,
          labels=substr(spnames,start=1,stop=4),layout="spring",diag=F,
       legend.cex=0.5,label.prop=0.99,
       groups=classnames,
       color=c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray","royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black"))
mtext(side=3,text="Correlations",line=2)

if (quartz) quartz() else x11()
qgraph::qgraph(g2.cor3,vsize=5,rescale=T,repulsion=0.7,
          labels=substr(spnames3,start=1,stop=4),layout="spring",diag=F,
       legend.cex=0.5,label.prop=0.99,
      groups=classnames,
       color=c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray","royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black"))
mtext(side=3,text="Clusters",line=2)

networks <- data.frame(metadata,taxatable)
for(i in names(table(clus)[table(clus)>1])) networks[,paste('cluster',i,sep="")] <- rowSums(networks[, names(clus)[clus==i]],na.rm=T)
for(i in names(table(clus)[table(clus)==1])) networks[,paste('cluster',i,sep="")] <- networks[, names(clus)[clus==i]]
networks <- list(networks, clus)
write.table(networks[[1]], file = "Clusters.txt", quote=F, sep="\t")
return(networks)
}

