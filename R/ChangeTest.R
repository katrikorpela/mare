ChangeTest <- function(taxonomic.table, meta, group = NULL, compare.to = NULL, 
    covariate = NULL, readcount.cutoff = 0, confounders = NULL, subject.ID,  time,
    outlier.cutoff = 3, p.cutoff = 0.05, select.by = NULL, select = NULL, pdf = F, 
    consecutive = T, min.prevalence = 0, label.direction = 1) {
    
    taxa <- read.delim(taxonomic.table)
    taxa <- taxa[, names(colSums(taxa > 0)[colSums(taxa > 0) > (min.prevalence * 
        nrow(taxa))])]
    
    metadata <- read.delim(meta)
    
    if (length(select.by) != 0) {
        metadata$selection <- metadata[, select.by]
        taxa <- taxa[metadata$selection == select, ]
        metadata <- metadata[metadata$selection == select, ]
    }
    
    metadata$ID <- metadata[, subject.ID]
    metadata$time <- as.numeric(metadata[, time])
    
    reltaxa <- (1 + taxa)/metadata$ReadCount
    
    deltataxa <- matrix(nrow = nrow(reltaxa), ncol = ncol(reltaxa))
    rownames(deltataxa) <- rownames(reltaxa)
    colnames(deltataxa) <- colnames(reltaxa)
    deltataxa <- as.data.frame(deltataxa)
    
    if (consecutive) {
        for (k in colnames(deltataxa)) deltataxa[, paste("baseline", k, sep = "")] <- NA
        for (i in names(table(metadata$ID)[table(metadata$ID) > 1])) {
            metadata$time[metadata$ID == i] <- order(metadata$time[metadata$ID == i])
            for (j in unique(metadata$time[metadata$ID == i][order(metadata$time[metadata$ID == 
                i])])[-1]) {
                for (k in colnames(reltaxa)) {
                  deltataxa[metadata$ID == i & metadata$time == j, k] <- reltaxa[metadata$ID == 
                    i & metadata$time == j, k] - reltaxa[metadata$ID == i & 
                    metadata$time == (j - 1), k]
                  tryCatch(deltataxa[metadata$ID == i & metadata$time == j, 
                    paste("baseline", k, sep = "")] <- reltaxa[metadata$ID == 
                    i & metadata$time == (j - 1), k], error = function(e) NULL)
                }
            }
        }
    } else {
        for (k in colnames(deltataxa)) deltataxa[, paste("baseline", k, sep = "")] <- NA
        for (i in names(table(metadata$ID)[table(metadata$ID) > 1])) {
            #metadata$time[metadata$ID == i] <- order(metadata$time[metadata$ID == i])
            for (j in unique(metadata$time[metadata$ID == i][order(metadata$time[metadata$ID == i])])[-1]) {
                for (k in colnames(reltaxa)) {
                  deltataxa[metadata$ID == i & metadata$time == j, k] <- reltaxa[metadata$ID == 
                    i & metadata$time == j, k] - reltaxa[metadata$ID == i & 
                    metadata$time == 1, k]
                  tryCatch(deltataxa[metadata$ID == i & metadata$time == j, 
                    paste("baseline", k, sep = "")] <- reltaxa[metadata$ID == 
                    i & metadata$time == 1, k], error = function(e) NULL)
                }
            }
        }
    }
    
    deltataxa <- deltataxa[metadata$time != levels(as.factor(metadata$time))[1], ]
    
    for (i in names(reltaxa)) {
        for (j in rownames(deltataxa)) {
            deltataxa[j, i][deltataxa[j, i] > (mean(deltataxa[, i]) + outlier.cutoff * 
                sd(deltataxa[, i]))] <- mean(deltataxa[, i]) + outlier.cutoff * 
                sd(deltataxa[, i])
            deltataxa[j, i][deltataxa[j, i] < (mean(deltataxa[, i]) - outlier.cutoff * 
                sd(deltataxa[, i]))] <- mean(deltataxa[, i]) - outlier.cutoff * 
                sd(deltataxa[, i])
        }
    }
    
    dataset <- data.frame(metadata[metadata$time != levels(as.factor(metadata$time))[1], ], deltataxa)
    dataset <- dataset[dataset$ReadCount > readcount.cutoff, ]
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    
    if (length(group) != 0 & length(covariate) == 0) {
        dataset[, group] <- dataset[, group][drop = T]
        dataset$G <- as.factor(dataset[,group])
        dataset[, group] <- as.character(dataset[, group])
        if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
        dataset[, group][dataset[, group] == compare.to] <- "00"
        dataset[, group] <- as.factor(dataset[, group])
        
        grouptime <- paste(dataset[dataset[, time] != 1 & dataset[, group] != 
            "00", time], dataset[dataset[, time] != 1 & dataset[, group] != 
            "00", group], sep = "/")
        
        group_test <- data.frame(array(dim = c(length(names(taxa)), 1 + (length(levels(as.factor(grouptime)))))))
        rownames(group_test) <- names(taxa)
        names(group_test) <- c("taxon", levels(as.factor(grouptime)))
        for (i in names(taxa)) {
            for (k in unique(dataset[, time])) {
                tryCatch(group_test[i, paste(k, levels(as.factor(dataset[, 
                  group]))[-1], sep = "/")] <- summary(lm(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+ as.factor(", group, ")", sep = "")), data = na.omit(dataset[dataset[, 
                  time] == k, c(group, confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], i, paste("baseline", i, sep = ""))[c(group, 
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], i, paste("baseline", i, sep = "")) != ""]])))$coef[-c(1:(2 + 
                  length(confounders[confounders != ""]))), 4], error = function(e) NULL)
            }
        }
        
        sig <- rownames(group_test)[sapply(data.frame(t(group_test)), min, 
            na.rm = T) < p.cutoff]
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
        group_test$taxon <- rownames(group_test)
      for(i in levels(dataset$G)){
      for(j in rownames(group_test)){
        group_test[j,paste("Mean",i,sep="_")] <-  mean(dataset[dataset$G==i,j]/dataset$ReadCount[dataset$G==i])
        }}
    for(i in levels(dataset$G)[levels(dataset$G)!= compare.to]) {
      group_test[,paste("Difference",i,sep="_")] <- abs(group_test[,paste("Mean",i,sep="_")] - group_test[,paste("Mean",compare.to,sep="_")])
    }
        write.table(group_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
            "_", "ChangeTest_", group, compare.to, "_", select.by, select, ".txt", 
            sep = ""), quote = F, row.names = F, sep = "\t")
   
        if (length(sig) > 0) {     
        taxa2 <- deltataxa
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_NA", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_1", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_2", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_3", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_4", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) gsub("_5", ".", x))
        names(taxa2) <- sapply(names(taxa2), function(x) strsplit(x, split = "_", 
            fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
        dataset2 <- data.frame(metadata[metadata$time != 1, ], taxa2)
        dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, ]
        if (length(select.by) != 0) {
            dataset2$selection <- dataset2[, select.by]
            dataset2 <- dataset2[dataset2$selection == select & !is.na(dataset2$selection), 
                ]
            dataset2[, group] <- dataset2[, group][drop = T]
        }
        
            if (pdf) {
                pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_change_", 
                  group, compare.to, "_", select.by, select, "_Boxplot.pdf", 
                  sep = ""))
                par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                  1), mgp = c(3, 0.5, 0), mar = c(5, 5, 1, 1), tck = -0.01, 
                  cex.axis = 1.5, cex.lab = 1.5)
                for (i in sig) {
                  boxplot(dataset2[, i] ~ paste(dataset2[dataset2[, time] != 
                    1, time], dataset2[, group], sep="/"), ylab = i, xlab = group, col = c("skyblue", 
                    "yellowgreen", "pink", "turquoise2", "plum", "darkorange", 
                    "lightyellow", "gray")[1:length(unique(dataset2[, group]))], 
                    outpch = 21, outbg = c("skyblue", "yellowgreen", "pink", 
                      "turquoise2", "plum", "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                      group]))], outcol = c("royalblue", "olivedrab4", "red", 
                      "turquoise4", "purple", "darkorange3", "lightyellow4", 
                      "black")[1:length(unique(dataset2[, group]))], boxcol = c("royalblue", 
                      "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                      "lightyellow4", "black")[1:length(unique(dataset2[, group]))], 
                    medcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                      "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                      group]))], whiskcol = c("royalblue", "olivedrab4", "red", 
                      "turquoise4", "purple", "darkorange3", "lightyellow4", 
                      "black")[1:length(unique(dataset2[, group]))], staplecol = c("royalblue", 
                      "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                      "lightyellow4", "black")[1:length(unique(dataset2[, group]))], las=label.direction)
                }
                dev.off()
                pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_change_", 
                  group, compare.to, "_", select.by, select, "_Beanplot.pdf", 
                  sep = ""))
                par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                  1), mgp = c(3, 0.5, 0), mar = c(3.5, 5, 1, 1), tck = -0.01, 
                  cex.axis = 1.5, cex.lab = 1.5)
                for (i in sig) {
                  tryCatch(beanplot::beanplot(dataset2[, i] ~ paste(dataset2[, 
                    time], dataset2[, group], sep="/"), ll = 0.1, ylab = i, xlab = group, 
                    col = list(c("skyblue", "royalblue", "royalblue", "royalblue"), 
                      c("yellowgreen", "olivedrab4", "olivedrab4", "olivedrab4"), 
                      c("pink", "red", "red", "red"), c("turquoise2", "turquoise4", 
                        "turquoise4", "turquoise4"), c("plum", "purple", "purple", 
                        "purple"), c("darkorange", "darkorange3", "darkorange3", 
                        "darkorange3"), c("lightyellow", "lightyellow4", "lightyellow4", 
                        "lightyellow4"), c("gray", "black", "black", "black"))[1:length(unique(dataset2[, 
                      group]))], border = c("royalblue", "olivedrab4", "red", 
                      "turquoise4", "purple", "darkorange3", "lightyellow4", 
                      "black")[1:length(unique(dataset2[, group]))]), error = function(e) NULL)
                }
                dev.off()
            }
            
            quartz()
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(3, 0.5, 0), mar = c(5, 5, 1, 1), tck = -0.01, 
                cex.axis = 1.5, cex.lab = 1.5)
            for (i in sig) {
                boxplot(dataset2[, i] ~ paste(dataset2[dataset2[, time] != 
                  1, time], dataset2[, group], sep="/"), ylab = i, xlab = group, col = c("skyblue", 
                  "yellowgreen", "pink", "turquoise2", "plum", "darkorange", 
                  "lightyellow", "gray")[1:length(unique(dataset2[, group]))], 
                  outpch = 21, outbg = c("skyblue", "yellowgreen", "pink", 
                    "turquoise2", "plum", "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                    group]))], outcol = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, group]))], boxcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, group]))], 
                  medcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    group]))], whiskcol = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, group]))], staplecol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, group]))],las=label.direction)
            }
            quartz()
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(3, 0.5, 0), mar = c(3.5, 5, 1, 1), tck = -0.01, 
                cex.axis = 1.5, cex.lab = 1.5)
            for (i in sig) {
                tryCatch(beanplot::beanplot(dataset2[, i] ~ paste(dataset2[, 
                  time], dataset2[, group], sep="/"), ll = 0.1, ylab = i, xlab = group, 
                  col = list(c("skyblue", "royalblue", "royalblue", "royalblue"), 
                    c("yellowgreen", "olivedrab4", "olivedrab4", "olivedrab4"), 
                    c("pink", "red", "red", "red"), c("turquoise2", "turquoise4", 
                      "turquoise4", "turquoise4"), c("plum", "purple", "purple", 
                      "purple"), c("darkorange", "darkorange3", "darkorange3", 
                      "darkorange3"), c("lightyellow", "lightyellow4", "lightyellow4", 
                      "lightyellow4"), c("gray", "black", "black", "black"))[1:length(unique(dataset2[, 
                    group]))], border = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, group]))]), error = function(e) NULL)
            }
        }
    }
    
    if (length(covariate) != 0) {
        if (length(group) != 0) {
        dataset[, group] <- dataset[, group][drop = T]
        dataset$group <- as.factor(dataset[, group])
        
    covariate_test <- data.frame(array(dim = c(length(names(taxa)), 
                (1 + length(levels(dataset$group))))))
            names(covariate_test) <- c("taxon", paste(covariate, levels(dataset$group), 
                sep = "_"))
            covariate_test$taxon <- names(taxa)
            rownames(covariate_test) <- names(taxa)
            
            for (i in names(taxa)) {
                for (j in levels(dataset$group)) {
                  tryCatch(covariate_test[i, paste(covariate, j, sep = "_")] <- summary(lm(as.formula(paste(i, 
                    "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = na.omit(dataset[dataset$group == 
                    j, c(covariate, confounders[1], confounders[2], confounders[3], 
                    confounders[4], confounders[5], i, paste("baseline", i, 
                      sep = ""))[c(covariate, confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], i, paste("baseline", 
                      i, sep = "")) != ""]])))$coef[-c(1:(2 + length(confounders[confounders != 
                    ""]))), 4], error = function(e) NULL)
                }
            }
          for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, i])
          covariate_test[, 1] <- 1

          sig <- names(sapply(rownames(covariate_test), function(x) min(covariate_test[x, 
                -1], na.rm = T))[sapply(rownames(covariate_test), function(x) min(covariate_test[x, 
                -1], na.rm = T)) < p.cutoff])
            sig <- sapply(sig, function(x) gsub("_NA", ".", x))
            sig <- sapply(sig, function(x) gsub("_1", ".", x))
            sig <- sapply(sig, function(x) gsub("_2", ".", x))
            sig <- sapply(sig, function(x) gsub("_3", ".", x))
            sig <- sapply(sig, function(x) gsub("_4", ".", x))
            sig <- sapply(sig, function(x) gsub("_5", ".", x))
            sig <- sapply(sig, function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, 
                split = "_", fixed = T)[[1]])])
            sig <- as.character(sig)
            
            names(covariate_test)[-1] <- paste("p",names(covariate_test)[-1],sep="_")
            for (j in levels(dataset$group)) {
                covariate_test[, paste("p", covariate, j, "FDR", sep = "_")] <- p.adjust(covariate_test[, 
                  paste("p", covariate, j, sep = "_")], "fdr")
            }
            covariate_test$taxon <- rownames(covariate_test)
            write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                "_ChangeTest_", covariate, "_", group, "_", select.by, select, 
                ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
            
            if (length(sig) > 0) {
                taxa2 <- deltataxa
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_NA",".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_1", ".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_2", ".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_3", ".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_4",".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) gsub("_5",".", x))
                names(taxa2) <- sapply(names(taxa2), function(x) strsplit(x, 
                  split = "_", fixed = T)[[1]][length(strsplit(x, split = "_", 
                  fixed = T)[[1]])])
                dataset2 <- data.frame(metadata[metadata$time != 1, ], taxa2)
                dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, 
                  ]
                if (length(select.by) != 0) {
                  dataset2$selection <- dataset2[, select.by]
                  dataset2 <- dataset2[dataset2$selection == select & !is.na(dataset2$selection), 
                    ]
                }
                df = na.omit(reshape2::melt(dataset2[, c(covariate, group, 
                  sig)], id = c(covariate, group)))
                names(df) <- c("x", "gr", "variable", "value")
                p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x, color = factor(gr)), 
                  environment = environment()) + ggplot2::stat_smooth(method = "lm", 
                  formula = y ~ x, ggplot2::aes(fill = factor(gr)), se = T) + 
                  ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(sig))), 
                    scales = "free") + ggplot2::theme_bw() + ggplot2::xlab(covariate) + 
                  ggplot2::ylab("Change in relative abundance") + ggplot2::theme(legend.position = "right") + 
                  ggplot2::scale_color_manual(name = group, values = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(levels(factor(df[, "gr"])))]) + 
                  ggplot2::scale_fill_manual(name = group, values = c("skyblue", 
                    "yellowgreen", "pink", "turquoise2", "plum", "darkorange", 
                    "lightyellow", "gray")[1:length(levels(factor(df[, "gr"])))])
                
                if (pdf) {
                  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                    "_change_", covariate, "_", group, "_", select.by, select, 
                    "_Covariateplot.pdf", sep = ""))
                  plot(p)
                  dev.off()
                }
                quartz()
                plot(p)
                
                
            }
        } else {
            covariate_test <- data.frame(taxon = names(taxa), p = rep(NA, length(names(taxa))))
            names(covariate_test)[2] <- covariate
            rownames(covariate_test) <- names(taxa)
            
            for (i in names(taxa)) {
                tryCatch(covariate_test[i, 2] <- summary(lm(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = na.omit(dataset[, c(covariate, 
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], i, paste("baseline", i, sep = ""))[c(covariate, 
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], i, paste("baseline", i, sep = "")) != ""]])))$coef[-c(1:(2 + 
                  length(confounders[confounders != ""]))), 4], error = function(e) NULL)
            }
            for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, 
                i])
            covariate_test[, 1] <- 1
            covariate_test <- na.omit(covariate_test)
            sig <- rownames(covariate_test)[covariate_test[, 2] < p.cutoff & 
                !is.na(covariate_test[, 2])]
            sig <- sapply(sig, function(x) gsub("_NA", ".", x))
            sig <- sapply(sig, function(x) gsub("_1", ".", x))
            sig <- sapply(sig, function(x) gsub("_2", ".", x))
            sig <- sapply(sig, function(x) gsub("_3", ".", x))
            sig <- sapply(sig, function(x) gsub("_4", ".", x))
            sig <- sapply(sig, function(x) gsub("_5", ".", x))
            sig <- sapply(sig, function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, 
                split = "_", fixed = T)[[1]])])
            sig <- as.character(sig)
            names(covariate_test)[-1] <- paste("p",names(covariate_test)[-1],sep="_")
            covariate_test[, paste("p", covariate, "FDR", sep = "_")] <- p.adjust(covariate_test[, 
                2], "fdr")
            covariate_test$taxon <- rownames(covariate_test)
            write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                "_ChangeTest_", covariate, "_", select.by, select, ".txt", sep = ""), 
                quote = F, row.names = F, sep = "\t")
            
            taxa2 <- deltataxa
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_NA", ".", x))
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_1", ".", x))
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_2", ".",  x))
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_3", ".",  x))
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_4", ".", x))
            names(taxa2) <- sapply(names(taxa2), function(x) gsub("_5", ".", x))
            names(taxa2) <- sapply(names(taxa2), function(x) strsplit(x, split = "_", 
                fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])])
            dataset2 <- data.frame(metadata[metadata$time != 1, ], taxa2)
            dataset2 <- dataset2[dataset2$ReadCount > readcount.cutoff, ]
            if (length(select.by) != 0) {
                dataset2$selection <- dataset2[, select.by]
                dataset2 <- dataset2[dataset2$selection == select & !is.na(dataset2$selection), 
                  ]
            }
            if (length(sig) > 0) {
                df = na.omit(reshape2::melt(dataset2[, c(covariate, sig)], 
                  id = covariate))
                names(df) <- c("x", "variable", "value")
                p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
                  ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = T, 
                    fill = "skyblue") + ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(sig))), 
                  scales = "free") + ggplot2::theme_bw() + ggplot2::xlab(covariate) + 
                  ggplot2::ylab("Change in relative abundance") + ggplot2::theme(legend.position = "right")
                
                if (pdf) {
                  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                    "_change_", covariate, "_", select.by, select, "Covariateplot.pdf", 
                    sep = ""))
                  plot(p)
                  dev.off()
                }
                
                quartz()
                plot(p)
            }
            
        }
    }
}
 
