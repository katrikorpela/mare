CommonTaxa <- function(taxonomic.table, mean.abundance, prevalence) {
    taxatable <- read.delim(taxonomic.table)
    abu <- colMeans(taxatable/rowSums(taxatable),na.rm=T)
    Abundant <- names(abu)[abu > mean.abundance]
    prev <- colSums(taxatable > 0)/nrow(taxatable)
    Prevalent <- names(prev)[prev > prevalence]
    ABU <- list(Abundant, Prevalent)
    names(ABU) <- c("Abundant", "Prevalent")
    return(ABU)
} 
