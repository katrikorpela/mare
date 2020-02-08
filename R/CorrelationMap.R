CorrelationMap <- function(taxonomic.table, meta, variables, select.by = NULL, 
                           selection = NULL,  outlier.cutoff = 3, readcount.cutoff = 0, 
                           min.abundance = 0, min.prevalence = 0, pdf = F, relative = T, logtrans = T){
if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}
 

taxa <- read.delim(taxonomic.table)
taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
  
 if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
  if(ncol(taxa)>0) {

metadata <- read.delim(meta)
taxa <- taxa[metadata$ReadCount > readcount.cutoff, ]
metadata <- metadata[metadata$ReadCount > readcount.cutoff, ]
    
if (length(select.by) != 0) {
        taxa <- taxa[metadata[, select.by] == selection, ]
        metadata <- metadata[metadata[, select.by] == selection, ]
    }
    
if(relative) reltaxa <- (1 + taxa)/metadata$ReadCount else reltaxa <- taxa + 1

for (i in names(reltaxa)) {
    for (j in 1:nrow(taxa)) {
         reltaxa[j, i][reltaxa[j, i] > (mean(reltaxa[, i]) + outlier.cutoff * sd(reltaxa[, i]))] <- mean(reltaxa[, i]) + outlier.cutoff * sd(reltaxa[, i])
        }
}

spnames <- names(reltaxa)

classnamesN <- rep(1,length(spnames))
if (taxonomic.table!="CAZy_table.txt"){classnames <- sapply(spnames, function(x) strsplit(x, split = "_", fixed = T)[[1]][1])
} else { classnames <- sapply(spnames, function(x) substr(x, start = 1, stop=2)[[1]][1])}
   classes <- levels(as.factor(classnames))
  classesN <- order(levels(as.factor(classnames)))   
  classnamesN <-classnames
  for(i in 1:length(classnames)) classnamesN[i] <-  classesN[classes==classnames[i]]


spnames <- sapply(spnames, function(x) gsub("_NA", ".", x))
spnames <- sapply(spnames, function(x) gsub("_1", ".", x))
spnames <- sapply(spnames, function(x) gsub("_2", ".", x))
spnames <- sapply(spnames, function(x) gsub("_3", ".", x))
spnames <- sapply(spnames, function(x) gsub("_4", ".", x))
spnames <- sapply(spnames, function(x) gsub("_5", ".", x))

spnames <- sapply(spnames, function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
names(reltaxa) <- spnames
metadata <- metadata[,variables]
 
if(logtrans) reltaxa <- log(reltaxa)

n<-nrow(reltaxa)
df<-n-2
correl <- cor(reltaxa,metadata,use="pairwise.complete.obs")[,colSums(abs(cor(reltaxa,metadata,use="pairwise.complete.obs")),na.rm=T)>0]
correl[is.na(correl)]<-0
pt2 <- function(q,df,log.p=F) 2*pt(-abs(q),df,log.p=log.p)
tstat<-correl*sqrt((n-2)/(1-correl^2))
correl.p<-pt2(tstat,df)
correl.sym <- correl.p
correl.sym[correl.p<0.05&correl.p>0.00999999]<-"*"
correl.sym[correl.p<0.01&correl.p>0.000999999]<-"**"
correl.sym[correl.p<0.001]<-"***"
correl.sym[correl.p>0.05]<-""

palette(c("#E41A1C","#FFA500","#377EB8","#87CEFA","#4DAF4A" ,'#9ACD32',"#984EA3",'#DA70D6', "#999999","gainsboro",
      "#008080","#00CED1","#F781BF","thistle1","#8DA0CB","lightsteelblue1","#FFD92F","#FFFFB3",
"#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
"#CCEBC5","#FFED6F","#C71585","#EE82EE","#66C2A5","#FC8D62","#A65628"))


if (pdf){
pdf(paste("CorrelationMap_",select.by,select,".pdf",sep=""))
gplots::heatmap.2(correl, col=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",
                  RowSideColors=classnamesN,
         cellnote=correl.sym,notecol = "black", keysize=1,key.xlab = "Correlation",margins=c(10,10),
         colRow=as.numeric(classnamesN))
dev.off()
}
quartz();par(mar=c(2,2,2,5)) 
gplots::heatmap.2(correl, col=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",
                  RowSideColors=classnamesN,
         cellnote=correl.sym,notecol = "black", keysize=1,key.xlab = "Correlation",margins=c(10,10),
         colRow=as.numeric(classnamesN))

}
}
    
