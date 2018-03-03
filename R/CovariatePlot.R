CovariatePlot <- function(meta, taxonomic.table, covariate, taxa=NULL, smooth.method = "loess", 
    readcount.cutoff = 0, select.by = NULL, select = NULL, pdf = F, logtrans = F, nozeros = F) {

   if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}
   
    meta <- read.delim(meta)
    taxatable <- read.delim(taxonomic.table)
   
    names(taxatable) <- sapply(names(taxatable), function(x) gsub("_NA", ".", x))
    names(taxatable) <- sapply(names(taxatable), function(x) strsplit(x, split = "_")[[1]][length(strsplit(x, split = "_")[[1]])])
    if(length(taxa)==0) taxa <- colnames(taxatable)
    taxa <- sapply(taxa, function(x) gsub("_NA", ".", x))
    taxa <- sapply(taxa, function(x) strsplit(x, split = "_")[[1]][length(strsplit(x, split = "_")[[1]])])
    
    if(nozeros) taxatable[taxatable==0]<-NA
    
    if(logtrans) { dataset <- data.frame(meta, log(((taxatable+1)/meta$ReadCount)))
    } else dataset <- data.frame(meta, (taxatable/meta$ReadCount) * 100)

    if (length(select.by) != 0) {
        dataset$selection <- dataset[, select.by]
        dataset <- dataset[dataset$selection == select, ]
    }
    dataset <- dataset[dataset$ReadCount > readcount.cutoff, ]
    
    
    
    df = na.omit(reshape2::melt(dataset[, c(covariate, taxa)], id = c(covariate)))
    names(df) <- c("x", "variable", "value")
    p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
          ggplot2::stat_smooth(method = smooth.method, formula = y ~ x, se = T,color="gray30",  fill = "cornflowerblue") + 
          ggplot2::geom_point(pch=20,color="gray30")+
          ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(taxa))), scales = "free") + 
          ggplot2::theme_bw() + 
          ggplot2::xlab(covariate) + 
          ggplot2::ylab("Relative abundance (%)") + 
          ggplot2::theme(legend.position = "none",strip.background =  ggplot2::element_rect(color = "white",fill="white"))
          
 quartz()
    plot(p)
    
    if (pdf) {
        pdf(paste(covariate,"Plot.pdf",sep=""))
       plot(p)
        dev.off()
    }
} 
