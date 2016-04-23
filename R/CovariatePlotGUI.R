CovariatePlotGUI <- function(taxonomic.table, meta, mean.abundance = NULL, prevalence = NULL,
                             covariate,  smooth.method = "loess", 
                    readcount.cutoff = 0, select.by = NULL, select = NULL, pdf = F, quartz = T) {
    
  h <- help("CovariatePlot"); print(h)
    
    meta <- read.delim(meta)
    taxatable <- read.delim(taxonomic.table)
    names(taxatable) <- sapply(names(taxatable), function(x) gsub("_NA", ".", 
        x))
    names(taxatable) <- sapply(names(taxatable), function(x) strsplit(x, split = "_")[[1]][length(strsplit(x, 
        split = "_")[[1]])])
    
    if(length(mean.abundance)!=0){
    abu <- colMeans(taxatable/rowSums(taxatable))
    Abundant <- names(abu)[abu > mean.abundance] 
    taxa <- Abundant 
    }  
   
  if(length(prevalence)!=0){
    prev <- colSums(taxatable > 0)/nrow(taxatable)
    Prevalent <- names(prev)[prev > prevalence]
    taxa <- Prevalent 
  }
    
    taxa <- sapply(taxa, function(x) gsub("_NA", ".", x))
    taxa <- sapply(taxa, function(x) strsplit(x, split = "_")[[1]][length(strsplit(x, 
        split = "_")[[1]])])
    
    dataset <- data.frame(meta, (taxatable/meta$ReadCount) * 100)
    if (length(select.by) != 0) {
        dataset$selection <- dataset[, select.by]
        dataset <- dataset[dataset$selection == select, ]
    }
    dataset <- dataset[dataset$ReadCount > readcount.cutoff, ]
    
    df = na.omit(reshape2::melt(dataset[, c(covariate, taxa)], id = c(covariate)))
    names(df) <- c("x", "variable", "value")
    p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
        ggplot2::stat_smooth(method = smooth.method, formula = y ~ x, se = T, 
            fill = "skyblue") + ggplot2::geom_point(pch=20,color="royalblue")+
      ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(taxa))), 
        scales = "free") + ggplot2::theme_bw() + ggplot2::xlab(covariate) + 
        ggplot2::ylab("Relative abundance (%)") + ggplot2::theme(legend.position = "none")  
    
    if (quartz) quartz()
    
    plot(p)
    
    if (pdf) {
        pdf("LinePlot.pdf")
        
        df = na.omit(reshape2::melt(dataset[, c(covariate, taxa)], id = c(covariate)))
        names(df) <- c("x", "variable", "value")
        p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
            ggplot2::stat_smooth(method = smooth.method, formula = y ~ x, se = T, 
                fill = "skyblue") + ggplot2::geom_point(pch=20,color="royalblue")+
          ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(taxa))), 
            scales = "free") + ggplot2::theme_bw() + ggplot2::xlab(covariate) + 
            ggplot2::ylab("Relative abundance (%)") + ggplot2::theme(legend.position = "none")
        plot(p)
        dev.off()
    }
} 
