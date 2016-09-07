CorrelationMap <- function(taxonomic.table, meta, variables, select.by = NULL, 
                           selection = NULL,  outlier.cutoff = 3, readcount.cutoff = 0, 
                           min.abundance = 0, min.prevalence = 0, quartz = T, pdf = F){
    
taxa <- read.delim(taxonomic.table)
taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
metadata <- read.delim(meta)
taxa <- taxa[metadata$ReadCount > readcount.cutoff, ]
metadata <- metadata[metadata$ReadCount > readcount.cutoff, ]
    
if (length(select.by) != 0) {
        metadata$selection <- metadata[, select.by]
        taxa <- taxa[metadata$selection == select, ]
        metadata <- metadata[metadata$selection == select, ]
    }
    
reltaxa <- (1 + taxa)/metadata$ReadCount
for (i in names(reltaxa)) {
    for (j in 1:nrow(taxa)) {
         reltaxa[j, i][reltaxa[j, i] > (mean(reltaxa[, i]) + outlier.cutoff * 
                sd(reltaxa[, i]))] <- mean(reltaxa[, i]) + outlier.cutoff * 
                sd(reltaxa[, i])
        }
}

spnames <- names(reltaxa)
spnames <- sapply(spnames, function(x) gsub("_NA", ".", x))
spnames <- sapply(spnames, function(x) gsub("_1", ".", x))
spnames <- sapply(spnames, function(x) gsub("_2", ".", x))
spnames <- sapply(spnames, function(x) gsub("_3", ".", x))
spnames <- sapply(spnames, function(x) gsub("_4", ".", x))
spnames <- sapply(spnames, function(x) gsub("_5", ".", x))


classnamesN <- rep(1,length(spnames))
if (length(strsplit(names(reltaxa)[1], split = "_", fixed = T)[[1]])>1){
   classnames <- sapply(spnames, function(x) strsplit(x, split = "_", fixed = T)[[1]][2])
   classes <- levels(as.factor(classnames))
  classesN <- order(levels(as.factor(classnames)))   
  classnamesN <-classnames
  for(i in 1:length(classnames)) classnamesN[i] <-  classesN[classes==classnames[i]]
}

spnames <- sapply(spnames, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
names(reltaxa) <- spnames
metadata <- metadata[,variables]
 
palette(c("black","skyblue","yellowgreen", "turquoise2", "plum", "darkorange", "gray","royalblue", "olivedrab4", "red", 
                      "turquoise4", "purple", "darkorange3", "lightyellow4"))
if (quartz) 
  quartz()
gplots::heatmap.2(cor(log(reltaxa+0.0001),metadata,use="pairwise.complete.obs"),
                  col=rainbow(256, start=0,end=0.34),density.info = "none",trace="none",
          keysize=1,key.xlab = "Correlation",margins=c(10,10),colRow=as.numeric(classnamesN))

if (pdf)
pdf("CorrelatioMap.pdf")
 gplots::heatmap.2(cor(log(reltaxa+0.0001),metadata,use="pairwise.complete.obs"),
                  col=rainbow(256, start=0,end=0.34),density.info = "none",trace="none",
          keysize=1,key.xlab = "Correlation",margins=c(10,10),colRow=as.numeric(classnamesN))
dev.off()
}

    
