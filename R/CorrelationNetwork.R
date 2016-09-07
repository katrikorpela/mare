CorrelationNetwork <- function(taxonomic.table, meta, N.taxa = NULL, minimum.correlation = 0.5, 
                               minimum.network = 1,  cluster.similarity = 1, quartz = T){

taxatable <- read.delim(taxonomic.table)
meta <- read.delim(meta)
taxatable <- taxatable/meta$ReadCount
pt2 <- function(q,df,log.p=F) 2*pt(-abs(q),df,log.p=log.p)
if (length(N.taxa) == 0) N.taxa = ncol(taxatable)
vars <- c(rev(names(colSums(taxatable,na.rm=T)[order(colSums(taxatable,na.rm=T))])))[1:N.taxa]
gs<-taxatable[,vars]
n<-nrow(gs)
df<-n-2
tgs <-data.frame(t(scale(gs)))
g2.cor<-cor(t(tgs[,1:n]),method="spearman",use="pairwise.complete.obs")
g2.tstat<-g2.cor*sqrt((n-2)/(1-g2.cor^2))#t = r*sqrt((n-2)/(1-r^2))
g2.p<-pt2(g2.tstat,df)
g2.adjp <- p.adjust(g2.p,method="fdr")
dim(g2.adjp) <- c(length(vars),length(vars))
g2.cor2 <- g2.cor;g2.cor2[g2.p>0.05] <- 0
g2.cor2[abs(g2.cor2)<minimum.correlation] <- 0 
g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]
g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]
g2.cor2<-g2.cor2[vegan::specnumber(g2.cor2)>minimum.network,vegan::specnumber(g2.cor2)>minimum.network]

spnames <- rownames(g2.cor)
spnames <- sapply(spnames, function(x) gsub("_NA", ".", x))
spnames <- sapply(spnames, function(x) gsub("_1", ".", x))
spnames <- sapply(spnames, function(x) gsub("_2", ".", x))
spnames <- sapply(spnames, function(x) gsub("_3", ".", x))
spnames <- sapply(spnames, function(x) gsub("_4", ".", x))
spnames <- sapply(spnames, function(x) gsub("_5", ".", x))
classnames <- sapply(spnames, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
spnames <- sapply(spnames, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])

if (quartz) quartz()
qgraph::qgraph(g2.cor2,vsize=5,rescale=T,labels=substr(spnames,start=1,stop=4),layout="spring",diag=F,
       legend.cex=0.5,label.prop=0.99,groups=classnames,min=0.2,
       color=c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray","royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black"))

clusters <- hclust(as.dist(1-g2.cor),"ward.D2")
if (quartz) quartz(); plot(clusters, xlab="", ylab="",labels=spnames); abline(h=cluster.similarity, lty=2, col="gray",xlab="")
clus<-cutree(clusters,h=cluster.similarity)
networks <- data.frame(meta,taxatable)
for(i in names(table(clus)[table(clus)>1])) networks[,paste('network',i,sep="")] <- rowMeans(networks[, names(clus)[clus==i]],na.rm=T)
networks <- list(networks, clus)
write.table(networks[[1]], file = "Networks.txt", quote=F, sep="\t")
return(networks)
}

