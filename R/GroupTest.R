GroupTest <- function(taxonomic.table, meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
    subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, zinf.cutoff = 0, 
    select.by = NULL, select = NULL, pdf = F,  min.prevalence = 0, min.abundance = 0, label.direction = 1) {
    
    taxa <- read.delim(taxonomic.table)
    taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
    if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {
    meta <- read.delim(meta)
    taxa <- taxa[meta$ReadCount > readcount.cutoff, ]
    meta <- meta[meta$ReadCount > readcount.cutoff, ]
    
    if (length(select.by) != 0) {
        meta$selection <- meta[, select.by]
        taxa <- taxa[meta$selection == select, ]
        meta <- meta[meta$selection == select, ]
    }
    
    reltaxa <- (1 + taxa)/meta$ReadCount
    for (i in names(reltaxa)) {
        for (j in 1:nrow(taxa)) {
            reltaxa[j, i][reltaxa[j, i] > (mean(reltaxa[, i]) + outlier.cutoff * 
                sd(reltaxa[, i]))] <- mean(reltaxa[, i]) + outlier.cutoff * 
                sd(reltaxa[, i])
        }
    }
    taxa <- round(reltaxa * meta$ReadCount - 1)
    taxa[taxa<0]<-0
    dataset <- data.frame(meta, taxa)
    dataset[, group] <- dataset[, group][drop = T]
    if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
    dataset$G <- as.factor(dataset[,group])
    dataset[, group] <- as.character(dataset[, group])
    dataset[, group][dataset[, group] == compare.to] <- "00"
    dataset[, group] <- as.factor(dataset[, group])
    
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    group_test <- data.frame(array(dim = c(length(names(taxa)), (length(levels(as.factor(dataset[,group])))))))
    rownames(group_test) <- names(taxa)
    names(group_test) <- levels(as.factor(dataset[, group]))
    if (length(subject.ID) != 0) {
        dataset$ID <- as.factor(dataset[, subject.ID])
        for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * 
            nrow(taxa)) - 1)]) {
            tryCatch(group_test[i, -1] <- summary(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ as.factor(", group, 
                ")+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = na.omit(dataset[, c(group, confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], i, "ID", 
                  "ReadCount")[c(group, confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], i, "ID", "ReadCount") != ""]]), 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)))$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4],
                #[-c(1:(1 + length(confounders[confounders != ""]))), 4], 
                error = function(e) NULL)
        }
        for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff * 
            nrow(taxa))]) {
            tryCatch(group_test[i, -1] <- summary(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ as.factor(", group, 
                ")+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = na.omit(dataset[, c(group, confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], i, "ID", 
                  group, "ReadCount")[c(group, confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], i, "ID", 
                  group, "ReadCount") != ""]]), zeroInflation = TRUE, admb.opts = glmmADMB::admbControl(shess = F, 
                  noinit = FALSE)))$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4],
                #[-c(1:(1 + length(confounders[confounders != ""]))), 4], 
                error = function(e) NULL)
        }
    } else {
        for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * nrow(taxa)) - 1)]) {
            tryCatch(group_test[i, -1] <- summary(MASS::glm.nb(as.formula(paste("round(", 
                i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ as.factor(", group, 
                ")+", "offset(log(ReadCount))")), data = dataset), 
                control = glm.control(maxit = 5000))$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4],
                #[-c(1:(1 + length(confounders[confounders != ""]))), 4], 
                error = function(e) NULL)
        }
        for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff *  nrow(taxa))]) {
            tryCatch(group_test[i, -1] <- summary(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ as.factor(", group, 
                ")+", "offset(log(ReadCount))")), family = "nbinom", data = na.omit(dataset[, 
                c(group, confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], i, group, "ReadCount")[c(group, confounders[1], 
                  confounders[2], confounders[3], confounders[4], confounders[5], 
                  i, group, "ReadCount") != ""]]), zeroInflation = TRUE, admb.opts = glmmADMB::admbControl(shess = F, 
                noinit = FALSE)))$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4],
                #[-c(1:(1 + length(confounders[confounders != ""]))), 4], 
                error = function(e) NULL)
        }
        
    }
    group_test[, 1] <- 1
    group_test <- na.omit(group_test)
    sig <- rownames(group_test)[sapply(data.frame(t(group_test[, -1])), min) < p.cutoff]
    sig <- sapply(sig, function(x) gsub("_NA", ".", x))
    sig <- sapply(sig, function(x) gsub("_1", ".", x))
    sig <- sapply(sig, function(x) gsub("_2", ".", x))
    sig <- sapply(sig, function(x) gsub("_3", ".", x))
    sig <- sapply(sig, function(x) gsub("_4", ".", x))
    sig <- sapply(sig, function(x) gsub("_5", ".", x))
    sig <- sapply(sig, function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, 
        split = "_", fixed = T)[[1]])])
    sig <- as.character(sig)
    names(group_test)[-1] <- paste("p",names(group_test)[-1],sep="_")
    for (k in names(group_test)[-1]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[, 
        k], "fdr")
    names(group_test)[1] <- "taxon"
    group_test$taxon <- rownames(group_test)
    
    for(i in levels(dataset$G)){
      for(j in rownames(group_test)){
        group_test[j,paste("Mean",i,sep="_")] <-  mean(dataset[dataset$G==i,j]/dataset$ReadCount[dataset$G==i])
        }}
    for(i in levels(dataset$G)[levels(dataset$G)!= compare.to]) {
      group_test[,paste("FoldChange",i,sep="_")] <- group_test[,paste("Mean",i,sep="_")] / group_test[,paste("Mean",compare.to,sep="_")]
    }
    
    write.table(group_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
        "_GroupTest_", group, compare.to, "_", select.by, select, ".txt", sep = ""), 
        quote = F, row.names = F, sep = "\t")
    
    taxa2 <- taxa
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_NA", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_1", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_2", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_3", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_4", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) gsub("_5", ".", x))
    names(taxa2) <- sapply(names(taxa2), function(x) strsplit(x, split = "_", 
        fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
    dataset2 <- data.frame(meta, (taxa2/meta$ReadCount) * 100)
    dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, ]
    if (length(select.by) != 0) {
        dataset2$selection <- dataset2[, select.by]
        dataset2 <- dataset2[dataset2$selection == select, ]
        dataset2[, group] <- dataset2[, group][drop = T]
    }
    if (length(sig) > 0) {
        if (pdf) {
            pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_", 
                group, compare.to, "_", select.by, select, "_Barplot.pdf", 
                sep = ""))
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, 
                cex.axis = 1.3, cex.lab = 1.5)
            for (i in sig) {
                boxplot(dataset2[, i] ~ dataset2[, group], ylab = i, xlab = "", las = label.direction, 
                  col = c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                    group]))], outpch = 21, outbg = c("skyblue", "yellowgreen", 
                    "pink", "turquoise2", "plum", "darkorange", "lightyellow", 
                    "gray")[1:length(unique(dataset2[, group]))], outcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, group]))], 
                  boxcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    group]))], medcol = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, group]))], whiskcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, group]))], 
                  staplecol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    group]))])
            }
            dev.off()
            pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_", 
                group, compare.to, "_", select.by, select, "_Beanplot.pdf", 
                sep = ""))
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, 
                cex.axis = 1.3, cex.lab = 1.5)
            for (i in sig) {
                tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, group],xlab="", las = label.direction,
                  ll = 0.1, ylab = i, col = list(c("skyblue", 
                    "royalblue", "royalblue", "royalblue"), c("yellowgreen", 
                    "olivedrab4", "olivedrab4", "olivedrab4"), c("pink", "red", 
                    "red", "red"), c("turquoise2", "turquoise4", "turquoise4", 
                    "turquoise4"), c("plum", "purple", "purple", "purple"), 
                    c("darkorange", "darkorange3", "darkorange3", "darkorange3"), 
                    c("lightyellow", "lightyellow4", "lightyellow4", "lightyellow4"), 
                    c("gray", "black", "black", "black")), border = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")), error = function(e) NULL)
            }
            dev.off()
        }
        
        quartz()
        par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
            1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, cex.axis = 1.3, 
            cex.lab = 1.5)
        for (i in sig) {
            boxplot(dataset2[, i] ~ dataset2[, group], ylab = i,xlab="", las = label.direction,
                col = c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                  "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                  group]))], outpch = 21, outbg = c("skyblue", "yellowgreen", 
                  "pink", "turquoise2", "plum", "darkorange", "lightyellow", 
                  "gray")[1:length(unique(dataset2[, group]))], outcol = c("royalblue", 
                  "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                  "lightyellow4", "black")[1:length(unique(dataset2[, group]))], 
                boxcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                  "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  group]))], medcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                  "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  group]))], whiskcol = c("royalblue", "olivedrab4", "red", 
                  "turquoise4", "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  group]))], staplecol = c("royalblue", "olivedrab4", "red", 
                  "turquoise4", "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  group]))])
        }
        quartz()
        par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
            1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, cex.axis = 1.3, 
            cex.lab = 1.5)
        for (i in sig) {
            tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, group],xlab="", las = label.direction,
                ll = 0.1, ylab = i, col = list(c("skyblue", "royalblue", 
                  "royalblue", "royalblue"), c("yellowgreen", "olivedrab4", 
                  "olivedrab4", "olivedrab4"), c("pink", "red", "red", "red"), 
                  c("turquoise2", "turquoise4", "turquoise4", "turquoise4"), 
                  c("plum", "purple", "purple", "purple"), c("darkorange", 
                    "darkorange3", "darkorange3", "darkorange3"), c("lightyellow", 
                    "lightyellow4", "lightyellow4", "lightyellow4"), c("gray", 
                    "black", "black", "black")), border = c("royalblue", "olivedrab4", 
                  "red", "turquoise4", "purple", "darkorange3", "lightyellow4", 
                  "black")), error = function(e) NULL)
        }
    }
    return(group_test)
    }
} 
