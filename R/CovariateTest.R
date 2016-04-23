CovariateTest <- function(taxonomic.table, meta, covariate, readcount.cutoff = 0, 
    confounders = NULL, subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, 
    group = NULL, zinf.cutoff = 0, select.by = NULL, select = NULL, pdf = F, 
    min.prevalence = 0) {
    
    metadata <- read.delim(meta)
    taxa <- read.delim(taxonomic.table)
    taxa <- taxa[, colSums(taxa > 0, na.rm = T) > min.prevalence * nrow(taxa)]
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
    taxa <- round(reltaxa * metadata$ReadCount - 1)
    taxa[taxa<0]<-0
    dataset <- data.frame(metadata, taxa)
    
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    covariate_test <- data.frame(taxon = names(taxa), estimate = rep(NA, length(names(taxa))), p = rep(NA, length(names(taxa))))
    names(covariate_test)[2:3] <- c(paste(covariate,"estimate",sep="_"),paste(covariate,"p",sep="_"))
    rownames(covariate_test) <- names(taxa)
    
    if (length(group) != 0) {
        dataset$group <- as.factor(dataset[, group])
        dataset$group <- dataset$group[drop = T]
        covariate_test <- data.frame(array(dim = c(length(names(taxa)), (1 + 2*length(levels(dataset$group))))))
        names(covariate_test) <- c("taxon", c(paste(covariate, levels(dataset$group),"estimate", sep = "_"),
                                              paste(covariate, levels(dataset$group),"p", sep = "_")))
        covariate_test$taxon <- names(taxa)
        rownames(covariate_test) <- names(taxa)
        
        if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])
            if (max(table(dataset[, group], dataset[, subject.ID])) > 1) {
                for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * 
                  nrow(taxa)) - 1)]) {
                  for (j in levels(dataset$group)) {
                    tryCatch(covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                      i, "+1)~", confounders[1], "+", confounders[2], "+", 
                      confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1 | 
                      ID, family = "nbinom", data = na.omit(dataset[dataset$group == 
                      j, c(confounders[1], confounders[2], confounders[3], 
                      confounders[4], confounders[5], i, "ID", covariate, "ReadCount")[c(confounders[1], 
                      confounders[2], confounders[3], confounders[4], confounders[5], 
                      i, "ID", covariate, "ReadCount") != ""]]), admb.opts = glmmADMB::admbControl(shess = F, 
                      noinit = FALSE)))$coef[-c(1:(1 + length(confounders[confounders != 
                      ""]))), c(1,4)], error = function(e) NULL)
                  }
                }
                for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff * 
                  nrow(taxa))]) {
                  for (j in levels(dataset$group)) {
                    tryCatch(covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                      i, "+1)~", confounders[1], "+", confounders[2], "+", 
                      confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1 | 
                      ID, family = "nbinom", data = na.omit(dataset[dataset$group == 
                      j, c(confounders[1], confounders[2], confounders[3], 
                      confounders[4], confounders[5], i, "ID", covariate, "ReadCount")[c(confounders[1], 
                      confounders[2], confounders[3], confounders[4], confounders[5], 
                      i, "ID", covariate, "ReadCount") != ""]]), zeroInflation = TRUE, 
                      admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)))$coef[-c(1:(1 + 
                      length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
                  }
                }
            }
        } else {
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * 
                nrow(taxa)) - 1)]) {
                for (j in levels(dataset$group)) {
                  tryCatch(covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(MASS::glm.nb(as.formula(paste("(", 
                    i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate, 
                    "+", "offset(log(ReadCount))")), data = dataset[dataset$group == 
                    j, ], control = glm.control(maxit = 1000)))$coef[-c(1:(1 + 
                    length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
                }
            }
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff * 
                nrow(taxa))]) {
                for (j in levels(dataset$group)) {
                  tryCatch(covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                    i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate, 
                    "+", "offset(log(ReadCount))")), family = "nbinom", data = na.omit(dataset[dataset$group == 
                    j, c(confounders[1], confounders[2], confounders[3], confounders[4], 
                    confounders[5], i, covariate, "ReadCount")[c(confounders[1], 
                    confounders[2], confounders[3], confounders[4], confounders[5], 
                    i, covariate, "ReadCount") != ""]]), zeroInflation = TRUE, 
                    admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)))$coef[-c(1:(1 + 
                    length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
                }
            }
        }
        for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, i])
        covariate_test[, 1] <- 1
        sig <- names(sapply(rownames(covariate_test), function(x) min(covariate_test[x, 
            -c(1,2,4,6,8,10,12,14,16,18,20)], na.rm = T))[sapply(rownames(covariate_test), function(x) min(covariate_test[x, 
            -c(1,2,4,6,8,10,12,14,16,18,20)], na.rm = T)) < p.cutoff])
        sig <- sapply(sig, function(x) gsub("_NA", ".", x))
        sig <- sapply(sig, function(x) gsub("_1", ".", x))
        sig <- sapply(sig, function(x) gsub("_2", ".", x))
        sig <- sapply(sig, function(x) gsub("_3", ".", x))
        sig <- sapply(sig, function(x) gsub("_4", ".", x))
        sig <- sapply(sig, function(x) gsub("_5", ".", x))
        sig <- sapply(sig, function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, 
            split = "_", fixed = T)[[1]])])
        sig <- as.character(sig)
        for (j in levels(dataset$group)) {
            covariate_test[, paste(covariate, j, "FDR", sep = "_")] <- p.adjust(covariate_test[, 
                paste(covariate, j, "p", sep = "_")], "fdr")
        }
        covariate_test$taxon <- rownames(covariate_test)
        write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
            "_CovariateTest_", covariate, "_", group, "_", select.by, select, ".txt", sep = ""), 
            quote = F, row.names = F, sep = "\t")
        # taxa2 <- taxa names(taxa2) <- sapply(names(taxa2), function(x)
        # gsub('_NA', '.', x)) names(taxa2) <- sapply(names(taxa2), function(x)
        # strsplit(x, split = '_', fixed = T)[[1]][length(strsplit(x, split = '_',
        # fixed = T)[[1]])]) dataset2 <- data.frame(meta, (taxa2/rowSums(taxa2)) *
        # 100) dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, ]
        # dataset2$selection <- dataset2[,select.by] dataset2 <-
        # dataset2[dataset2$selection==select,]
        # dataset2[,group]<-dataset2[,group][drop=T]
        if (length(sig) > 0) {
            GroupPlot(taxa = sig, group = as.character(group), taxonomic.table = as.character(taxonomic.table), 
                meta = as.character(meta), bar = F, box = F, stacked = F, bean = F, 
                covariate = covariate, smooth.method = "lm", pdf = pdf, select.by = select.by, 
                select = select)
        }
    } else {
        if (length(subject.ID) != 0) {
            dataset$ID <- as.factor(dataset[, subject.ID])
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * 
                nrow(taxa)) - 1)]) {
                tryCatch(covariate_test[i, c(2,3)] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                  i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                  data = na.omit(dataset[, c(confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], i, "ID", 
                    covariate, "ReadCount")[c(confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], i, "ID", 
                    covariate, "ReadCount") != ""]]), admb.opts = glmmADMB::admbControl(shess = F, 
                    noinit = FALSE)))$coef[-c(1:(1 + length(confounders[confounders != 
                  ""]))), c(1,4)], error = function(e) NULL)
            }
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff * 
                nrow(taxa))]) {
                tryCatch(covariate_test[i, c(2,3)] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                  i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                  data = na.omit(dataset[, c(confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], i, "ID", 
                    covariate, "ReadCount")[c(confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], i, "ID", 
                    covariate, "ReadCount") != ""]]), zeroInflation = TRUE, 
                  admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)))$coef[-c(1:(1 + 
                  length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
            }
        } else {
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) > ((zinf.cutoff * 
                nrow(taxa)) - 1)]) {
                tryCatch(covariate_test[i, c(2,3)] <- summary(MASS::glm.nb(as.formula(paste("(", 
                  i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), data = dataset, control = glm.control(maxit = 1000)))$coef[-c(1:(1 + 
                  length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
            }
            for (i in names(taxa)[colSums(taxa > 0, na.rm = T) < (zinf.cutoff * 
                nrow(taxa))]) {
                tryCatch(covariate_test[i, c(2,3)] <- summary(glmmADMB::glmmadmb(as.formula(paste("(", 
                  i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), family = "nbinom", data = na.omit(dataset[, 
                  c(confounders[1], confounders[2], confounders[3], confounders[4], 
                    confounders[5], i, covariate, "ReadCount")[c(confounders[1], 
                    confounders[2], confounders[3], confounders[4], confounders[5], 
                    i, covariate, "ReadCount") != ""]]), zeroInflation = TRUE, 
                  admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)))$coef[-c(1:(1 + 
                  length(confounders[confounders != ""]))), c(1,4)], error = function(e) NULL)
            }
        }
        for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[,i])
        covariate_test[, 1] <- 1
        # covariate_test <- na.omit(covariate_test)
        sig <- rownames(covariate_test)[covariate_test[, 3] < p.cutoff & !is.na(covariate_test[, 3])]
        sig <- sapply(sig, function(x) gsub("_NA", ".", x))
        sig <- sapply(sig, function(x) gsub("_1", ".", x))
        sig <- sapply(sig, function(x) gsub("_2", ".", x))
        sig <- sapply(sig, function(x) gsub("_3", ".", x))
        sig <- sapply(sig, function(x) gsub("_4", ".", x))
        sig <- sapply(sig, function(x) gsub("_5", ".", x))
        sig <- sapply(sig, function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, 
            split = "_", fixed = T)[[1]])])
        sig <- as.character(sig)
        covariate_test[, paste(covariate, "FDR", sep = "_")] <- p.adjust(covariate_test[,3], "fdr")
        covariate_test$taxon <- rownames(covariate_test)
        write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
            "_CovariateTest_", covariate, "_", select.by, select, ".txt", sep = ""), quote = F, 
            row.names = F, sep = "\t")
        taxa2 <- taxa
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_NA", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_1", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_2", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_3", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_4", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_5", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
        dataset2 <- data.frame(metadata, (taxa2/metadata$ReadCount) * 100)
        dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, ]
        if (length(select.by) != 0) {
            dataset2$selection <- dataset2[, select.by]
            dataset2 <- dataset2[dataset2$selection == select, ]
        }
        if (length(sig) > 0) {
            df = na.omit(reshape2::melt(dataset2[, c(covariate, sig)], id = covariate))
            names(df) <- c("x", "variable", "value")
            p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
                ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = T, 
                  fill = "skyblue") +
              ggplot2::geom_point(pch=20,color = "royalblue") + ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(sig))), 
                scales = "free") + ggplot2::theme_bw() + ggplot2::xlab(covariate) + 
                ggplot2::ylab("Relative abundance (%)") + ggplot2::theme(legend.position = "right")
            quartz()
            plot(p)
            
            if (pdf) {
                pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_", 
                  covariate, "_", select.by, select, "_", "CovariatePlot.pdf", 
                  sep = ""))
                plot(p)
                dev.off()
            }
        }
    }
} 
