ChangeTest <- function(taxonomic.table, meta, group = NULL, compare.to = NULL, 
    covariate = NULL, readcount.cutoff = 0, confounders = NULL, subject.ID,  time,
    outlier.cutoff = 3, p.cutoff = 0.05, select.by = NULL, select = NULL, pdf = F, 
    consecutive = T, min.prevalence = 0, min.abundance = 0, label.direction = 1, quartz = T, keep.result =F) {
    
  taxa <- read.delim(taxonomic.table)
  colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    
  taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
    
  metadata <- read.delim(meta)
 

    if (length(select.by) != 0) {
        metadata$selection <- metadata[, select.by]
        taxa <- taxa[metadata$selection == select, ]
        metadata <- metadata[metadata$selection == select, ]
    }
    
    metadata$ID <- metadata[, subject.ID]
    metadata$time <- as.numeric(ordered(metadata[, time]))
    
    reltaxa <- (1 + taxa)/metadata$ReadCount
    
    deltataxa <- matrix(nrow = nrow(reltaxa), ncol = ncol(reltaxa))
    rownames(deltataxa) <- rownames(reltaxa)
    colnames(deltataxa) <- colnames(reltaxa)
    deltataxa <- as.data.frame(deltataxa)
    
    if (consecutive) {
        for (k in colnames(deltataxa)) deltataxa[, paste("baseline", k, sep = "")] <- NA
        for (i in names(table(metadata$ID)[table(metadata$ID) > 1])) {
            #metadata$time[metadata$ID == i] <- order(metadata$time[metadata$ID == i])
            for (j in unique(metadata$time[metadata$ID == i][order(metadata$time[metadata$ID == i])])[-1]) {
                for (k in colnames(reltaxa)) {
                  deltataxa[metadata$ID == i & metadata$time == j, k] <- reltaxa[metadata$ID == 
                    i & metadata$time == j, k] - reltaxa[metadata$ID == i &  metadata$time == (j - 1), k]
                  tryCatch(deltataxa[metadata$ID == i & metadata$time == j, 
                    paste("baseline", k, sep = "")] <- reltaxa[metadata$ID == i & metadata$time == (j - 1), k], 
                    error = function(e) NULL)
                }
            }
        }
    } else {
        for (k in colnames(deltataxa)) deltataxa[, paste("baseline", k, sep = "")] <- NA
        for (i in names(table(metadata$ID)[table(metadata$ID) > 1])) {
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
    
    model <- list()
    
    
    if (length(group) != 0 & length(covariate) == 0) {
      
    if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
    dataset$G <- as.factor(dataset[,group])
    dataset[, group] <- as.character(dataset[, group])
    dataset[, group][dataset[, group] == "0" & !is.na(dataset[, group])] <- "group0"
    dataset[, group][dataset[, group] == compare.to & !is.na(dataset[, group])] <- "0"
    dataset[, group] <- as.factor(dataset[, group])
    
        grouptime <- paste(dataset[dataset$time != 1 & dataset[, "G"] != 
            compare.to,"time"], dataset[dataset$time != 1 & dataset[, "G"] != 
            compare.to, "G"], sep = "/")
        
        modeldata <- na.omit(dataset[, c("time","G", subject.ID,time,group, confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], names(taxa), paste("baseline", names(taxa), sep = ""))[c("time", "G",subject.ID,group, time,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = "")) != ""]])
        
        for(i in unique(modeldata[,subject.ID])){ 
 modeldata[modeldata[,subject.ID]==i,names(taxa)[apply(taxa[metadata[,subject.ID]==i,names(taxa)], MARGIN=2,FUN=max)==0]]<-NA
  }
        group_test <- data.frame(array(dim = c(length(names(taxa)), 1 + (length(levels(as.factor(grouptime)))))))
        rownames(group_test) <- names(taxa)
        names(group_test) <- c("taxon", levels(as.factor(grouptime)))
        
        for (i in names(taxa)) {
            for (k in unique(modeldata$time)) {
               tryCatch(model[[paste(i,k)]] <- lm(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+ as.factor(", group, ")", sep = "")), data = modeldata[modeldata$time == k&!is.na(modeldata[,i]),]), 
                  error = function(e) NULL)
               if(length(model[[paste(i,k)]])>0){
                  if(nrow(summary(model[[paste(i,k)]])$coef)>1){
            if(summary(model[[paste(i,k)]])$coef[2,4]!="NaN"){
          tmp <- data.frame(res = resid(model[[paste(i,k)]],type="deviance"), pred = predict(model[[paste(i,k)]]))  
          tmp$group <- modeldata[modeldata$time == k&!is.na(modeldata[,i]),group]
          tmp$resdev <- abs(tmp$res-0) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
               group_test[i, paste(k, levels(as.factor(dataset[, "G"]))[levels(as.factor(dataset[, "G"]))!=compare.to],sep = "/")] <- summary(model[[paste(i,k)]])$coef[-c(1:(2 + 
                  length(confounders[confounders != ""]))), 4]
                } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
          } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
            } else print(paste("Model for",i,"showed violation against independence! You may be missing an important covariate."))
              } else print(paste("Model for",i,"showed violation against independence! You may be missing an important covariate."))
            } else print(paste("Model for",i,"failed!"))
               } else print(paste("Model for",i,"failed!"))
               } else print(paste("Model for",i,"failed!"))
            }
        }
        
     group_test$taxon <- rownames(group_test)
      group_test <- na.omit(group_test)  
      sig <- as.character(rownames(group_test)[apply(group_test[,-1],MARGIN=1,FUN=min,na.rm=T) < p.cutoff])
     
        names(group_test)[-1] <- paste("p",names(group_test)[-1],sep="_")
        for (k in names(group_test)[-1]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,  k], "fdr")

      resids <- modeldata
  
    for(i in rownames(group_test)){
      for(k in unique(modeldata$time)){
      resids[resids$time==k&!is.na(resids[,i]),i] <-  tryCatch(resid(update(model[[paste(i,k)]], as.formula(paste(".~. -as.factor(",group,")",sep="")))),
                               error = function(e) NULL)
    }}
 
    for(i in levels(resids$G)){
      for(j in rownames(group_test)){
        for(k in unique(resids$time)){
        group_test[j,paste("DevianceFromExpected",k,i,sep="_")] <-  mean(resids[resids$G==i&resids$time==k,j],na.rm=T)
      }
    }
    }

 write.table(group_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
            "_", "ChangeTest_", group, compare.to, "_", select.by, select, ".txt", 
            sep = ""), quote = F, row.names = F, sep = "\t")
   
    dataset2 <- data.frame(resids)
    dataset2[,group] <- dataset2$G
    dataset2$grouptime <-paste(dataset2[, time], dataset2[, group], sep="/")
      
        if (length(sig) > 0) {     
       
        
            if (pdf) {
                pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_change_", 
                  group, compare.to, "_", select.by, select, "_Boxplot.pdf", 
                  sep = ""))
                par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                  1), mgp = c(3.5, 0.5, 0), mar = c(3, 5, 1, 1), tck = -0.01, 
                  cex.axis = 1.1, cex.lab = 1.5)
                for (i in sig) {
                  # dataset3 <- na.omit(dataset2[,c(i,"grouptime",group,time)]) 
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                 #bxp(z=list(stats=matrix(unlist(tapply(dataset3[,i],dataset3[,"grouptime"],summary)),nrow=6,ncol=length(unique(dataset3[,"grouptime"])))[-3,],
                #           n=table(dataset3[,"grouptime"]),names=levels(dataset3[,"grouptime"])),
                    boxplot(dataset2[, i] ~ dataset2[, "grouptime"], 
                    ylab = lab, xlab = group, boxfill = c("skyblue", 
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
                  1), mgp = c(2, 0.5, 0), mar = c(3, 3.5, 1, 1), tck = -0.01, 
                  cex.axis = 1.1, cex.lab = 1.5)
                for (i in sig) {
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                  tryCatch(beanplot::beanplot(dataset2[, i] ~ paste(dataset2[, 
                    time], dataset2[, group], sep="/"), ll = 0.1, ylab = lab, xlab = group, 
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
            
            if (quartz) quartz() else x11()
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(3.5, 0.5, 0), mar = c(3, 5, 1, 1), tck = -0.01, 
                cex.axis = 1.1, cex.lab = 1.5)
            
           
            for (i in sig) { 
              #dataset3 <- na.omit(dataset2[,c(i,"grouptime",group,time)]) 
              lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
               # boxplot(dataset2[, i] ~ paste(dataset2[dataset2[, time] !=  1, time], dataset2[, group], sep="/"), 
              boxplot(dataset2[, i] ~ dataset2[, "grouptime"],   
              #bxp(z=list(stats=matrix(unlist(tapply(dataset2[,i],dataset2[,"grouptime"],summary)),nrow=6,ncol=length(unique(dataset2[,"grouptime"])))[-3,],
               #            n=table(dataset2[,"grouptime"]),names=levels(dataset2[,"grouptime"])),
                    ylab = lab, xlab = group, 
                  boxfill = c("skyblue",  "yellowgreen", "pink", "turquoise2", "plum","darkorange", 
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
            if (quartz) quartz() else x11()
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.5, 0), mar = c(3, 3.5, 1, 1), tck = -0.01, 
                cex.axis = 1.1, cex.lab = 1.5)
            for (i in sig) {
              lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                tryCatch(beanplot::beanplot(dataset2[, i] ~ paste(dataset2[, 
                  time], dataset2[, group], sep="/"), ll = 0.1, ylab = lab, xlab = group, 
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
       if(keep.result)   return(group_test)
    }
   
    if (length(covariate) != 0) {
        if (length(group) != 0) {
        dataset[, group] <- dataset[, group][drop = T]
        dataset$group <- as.factor(dataset[, group])
        
   covariate_test <- data.frame(array(dim = c(length(names(taxa)), (1 + 2*length(levels(dataset$group))))))
   names(covariate_test) <- c("taxon", c(paste(covariate, levels(dataset$group),"estimate", sep = "_"),
                                              paste(covariate, levels(dataset$group),"p", sep = "_")))
  covariate_test$taxon <- names(taxa)
  rownames(covariate_test) <- names(taxa)
        
  modeldata <- na.omit(dataset[, c("group",subject.ID, covariate, confounders[1], confounders[2], confounders[3], 
                    confounders[4], confounders[5], names(taxa), paste("baseline", names(taxa), 
                      sep = ""))[c("group", subject.ID, covariate, confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], names(taxa), paste("baseline", 
                      names(taxa), sep = "")) != ""]])
    
 for(i in unique(modeldata[,subject.ID])){ 
 modeldata[modeldata[,subject.ID]==i,names(taxa)[apply(taxa[metadata[,subject.ID]==i,names(taxa)], MARGIN=2,FUN=max)==0]]<-NA
  }  

    for (i in names(taxa)) {
    for (j in levels(dataset$group)) {
    tryCatch(model[[paste(i,j)]] <- lm(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = modeldata[modeldata$group==j&!is.na(modeldata[,i]),]),
                    error = function(e) NULL)
    if(length(model[[paste(i,j)]])>0 ){ 
    if(nrow(summary(model[[paste(i,j)]])$coef)>1){
    if(summary(model[[paste(i,j)]])$coef[covariate,4]!="NaN"&summary(model[[paste(i,j)]])$coef[covariate,4]!="NA"){
    tmp <- data.frame(res = resid(model[[paste(i,j)]],type="deviance"), pred = predict(model[[paste(i,j)]]))  
          tmp$covariate <- modeldata[modeldata$group==j&!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res-0) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
   covariate_test[i, c(paste(covariate, j, "estimate",sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$coef[covariate,c(1,4)]# -c(1:(2 + length(confounders[confounders !=""])))
                 } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
          } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
            } else print(paste("Model for",i,"in group",j,"showed violation against independence! You may be missing an important covariate."))
              } else print(paste("Model for",i,"in group",j,"showed violation against independence! You may be missing an important covariate."))
            } else print(paste("Model for",i,"in group",j,"failed!"))
             } else print(paste("Model for",i,"in group",j,"failed!"))
             } else print(paste("Model for",i,"in group",j,"failed!"))
                }
            }
          for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, i])
        #covariate_test <- na.omit(covariate_test)
         sig <- na.omit(names(apply(covariate_test[,grepl(pattern="_p",x=colnames(covariate_test))],MARGIN = 1,FUN = min)[apply(covariate_test[,grepl(pattern="_p",x=colnames(covariate_test))],MARGIN = 1,FUN = min)<p.cutoff]))
   
            for (j in levels(dataset$group)) {
                covariate_test[, paste(covariate, j, "p","FDR", sep = "_")] <- p.adjust(covariate_test[, 
                  paste(covariate, j,"p",  sep = "_")], "fdr")
            }

            write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                "_ChangeTest_", covariate, "_", group, "_", select.by, select, 
                ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
            
for(h in levels(modeldata[,"group"])){            
covariate_test_summary <- na.omit(covariate_test[,grepl(pattern=paste("_",h,"_",sep=""),x=names(covariate_test))])
covariate_test_summary$class <-  unlist(lapply(rownames(covariate_test_summary), function(x) strsplit(x, split="_")[[1]][2]))
covariate_test_summary <- covariate_test_summary[order(covariate_test_summary[,1]),]
palette(c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                  "darkorange", "lightyellow", "gray","royalblue", "olivedrab4", "red", "turquoise4", "purple", "darkorange3", "lightyellow4", "black"))

if (quartz) quartz() 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(covariate_test_summary[,1],horiz=T,names.arg = paste(rownames(covariate_test_summary),"p =",
      round(covariate_test_summary[,2],digits=2)),las=2,col=as.factor(covariate_test_summary$class),
      xaxt="n",xlab=paste("Estimated effect size of",covariate),main=paste(group, "=", h) )  
abline(v=0)
axis(side=1) 

if (pdf) {
  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_",group, "_", select.by,select, "_", "CovariateTestResult.pdf", sep = "")) 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(covariate_test_summary[,2],horiz=T,names.arg = paste(rownames(covariate_test_summary),"p =",
      round(covariate_test_summary[,3],digits=2)),las=2,col=as.factor(covariate_test_summary$class),
      xaxt="n",xlab=paste("Estimated effect size of",covariate))   
abline(v=0)
axis(side=1) 
dev.off()} 
}

       if (length(sig) > 0) {
      resids <- modeldata
    for(i in sig){
      for (j in levels(resids$group)) {
    tryCatch( 
      resids[resids$group==j&!is.na(resids[,i]),i] <-  resid(update(model[[paste(i,j)]],as.formula(paste(".~. -",covariate,sep="")))),
                              error = function(e) NULL)
    }
           }
    dataset2 <- data.frame(resids)
        
                df = na.omit(reshape2::melt(dataset2[, c(covariate, "group", sig)], id = c(covariate, "group")))
                names(df) <- c("x", "gr", "variable", "value")
                 df$variable <- sapply(df$variable, function(x) gsub(pattern = "_NA",replacement = "", x=x))
                df$variable <- sapply(df$variable, function(x) gsub(pattern = "_incertae_sedis",replacement = "", x=x))
                df$variable <- unlist(sapply(as.character(df$variable), function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])]))

   p<- ggplot2::ggplot(df, ggplot2::aes(y=value, x=x,color=factor(gr)),environment = environment()) +
  ggplot2::stat_smooth(method = "lm", formula = y ~ x, ggplot2::aes(fill=factor(gr)),se=T) +
  ggplot2::facet_wrap(~variable,ncol=floor(sqrt(length(sig))),scales='free')+
  ggplot2::theme_bw()+
  ggplot2::xlab(covariate)+
  ggplot2::geom_point(pch=20) + 
  ggplot2::scale_color_manual(name=group,values=c('gray30','red','turquoise4','olivedrab4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c('cornflowerblue','pink','turquoise2','yellowgreen','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])+
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
      ggplot2::ylab('Deviance from expected') 

                if (pdf) {
                  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                    "_change_", covariate, "_", group, "_", select.by, select, 
                    "_Covariateplot.pdf", sep = ""))
                  plot(p)
                  dev.off()
                }
                if (quartz) quartz() else x11()
                plot(p)
                
                
       }
        } else {
    covariate_test <- data.frame(taxon = names(taxa), estimate = rep(NA, length(names(taxa))), p = rep(NA, length(names(taxa))))
    names(covariate_test)[2:3] <- c(paste(covariate,"estimate",sep="_"),paste(covariate,"p",sep="_"))
    rownames(covariate_test) <- names(taxa)
       
            modeldata <- na.omit(dataset[, c(covariate, subject.ID,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = ""))[c(covariate, subject.ID,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = "")) != ""]])
 
            for(i in unique(modeldata[,subject.ID])){ 
 modeldata[modeldata[,subject.ID]==i,names(taxa)[apply(taxa[metadata[,subject.ID]==i,names(taxa)], MARGIN=2,FUN=max)==0]]<-NA
  } 
            for (i in names(taxa)) {
                tryCatch(model[[i]] <- lm(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = modeldata), 
                  error = function(e) NULL)
          if(length(model[[i]])>0){
          if(nrow(summary(model[[i]])$coef)>1){
         if(summary(model[[i]])$coef[2,4]!="NaN"){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), pred = predict(model[[i]]))  
          tmp$covariate <- modeldata[!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res-0) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                  covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$coef[-c(1:(2 + length(confounders[confounders != ""]))),c(1, 4)]
               } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
          } else print(paste("Model for",i,"showed violation against homoscedasticity!"))
            } else print(paste("Model for",i,"showed violation against independence! You may be missing an important covariate."))
              } else print(paste("Model for",i,"showed violation against independence! You may be missing an important covariate."))
            } else print(paste("Model for",i,"failed!"))
                   } else print(paste("Model for",i,"failed!"))
                   } else print(paste("Model for",i,"failed!"))
            }

            covariate_test <- na.omit(covariate_test)
           
          sig <-  rownames(covariate_test)[covariate_test[,3]<p.cutoff]
   
   
            covariate_test[, paste(covariate,"p",  "FDR", sep = "_")] <- p.adjust(covariate_test[,paste(covariate,"p", sep = "_")], "fdr")
            write.table(covariate_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                "_ChangeTest_", covariate, "_", select.by, select, ".txt", sep = ""), 
                quote = F, row.names = F, sep = "\t")


covariate_test_summary <- covariate_test
covariate_test_summary$class <-  unlist(lapply(rownames(covariate_test_summary), function(x) strsplit(x, split="_")[[1]][2]))
covariate_test_summary <- covariate_test_summary[order(covariate_test_summary[,2]),]
palette(c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                  "darkorange", "lightyellow", "gray","royalblue", "olivedrab4", "red", "turquoise4", "purple", "darkorange3", "lightyellow4", "black"))

if (quartz) quartz() 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(covariate_test_summary[,2],horiz=T,names.arg = paste(rownames(covariate_test_summary),"p =",
      round(covariate_test_summary[,3],digits=2)),las=2,col=as.factor(covariate_test_summary$class),
      xaxt="n",xlab=paste("Estimated effect size of",covariate))   
abline(v=0)
axis(side=1) 

if (pdf) {
  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3],"_",covariate,"_",group, "_", select.by,select, "_", "CovariateTestResult.pdf", sep = "")) 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(covariate_test_summary[,2],horiz=T,names.arg = paste(rownames(covariate_test_summary),"p =",
      round(covariate_test_summary[,3],digits=2)),las=2,col=as.factor(covariate_test_summary$class),
      xaxt="n",xlab=paste("Estimated effect size of",covariate))   
abline(v=0)
axis(side=1) 
dev.off()}

   
    if (length(sig) > 0) { 
    resids<- modeldata

    for(i in sig){
      tryCatch(resids[!is.na(resids[,i]),i] <-  resid(update(model[[i]],as.formula(paste(".~. -",covariate,sep="")))),
                               error = function(e) NULL)
    }
  
    dataset2 <- data.frame(resids)

                df = na.omit(reshape2::melt(dataset2[, c(covariate, sig)], id = covariate))
                names(df) <- c("x", "variable", "value")
                df$variable <- sapply(df$variable, function(x) gsub(pattern = "_NA",replacement = "", x=x))
                df$variable <- sapply(df$variable, function(x) gsub(pattern = "_incertae_sedis",replacement = "", x=x))
                df$variable <- unlist(sapply(as.character(df$variable), function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])]))
                
                p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
                  ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = T, fill = "cornflowerblue", color = "gray30", lwd = 0.5) + 
                  ggplot2::facet_wrap(~variable, ncol = floor(sqrt(length(sig))), scales = "free") + 
                  ggplot2::theme_bw() + 
                  ggplot2::xlab(covariate) + 
                  ggplot2::geom_point(pch=20) + 
  ggplot2::scale_color_manual(name=group,values=c('gray30','red','turquoise4','olivedrab4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c('cornflowerblue','pink','turquoise2','yellowgreen','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])+
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
      ggplot2::ylab('Deviance from expected') 
 
                
                if (pdf) {
                  pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
                    "_change_", covariate, "_", select.by, select, "Covariateplot.pdf", 
                    sep = ""))
                  plot(p)
                  dev.off()
                }
                
                if (quartz) quartz() else x11()
                plot(p)
            }
             if(keep.result) return(covariate_test)
        }
    }
}
 
