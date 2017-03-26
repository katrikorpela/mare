GroupTest <- function(species.table=NULL, genus.table=NULL, family.table=NULL, order.table=NULL, class.table=NULL, phylum.table=NULL, 
                      meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
                      subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, 
                      select.by = NULL, select = NULL, pdf = F,  min.prevalence = 0.1, min.abundance = 0, 
                      label.direction = 1, keep.result = F, nonzero = T) {
  
  if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}
  
  meta <- read.delim(meta)
   if(length(species.table)>0) species <- read.delim(species.table) else species <- NA
  if(length(genus.table)>0) genus <- read.delim(genus.table) else genus <- NA
  if(length(family.table)>0) family <- read.delim(family.table) else family <- NA
  if(length(order.table)>0) order <- read.delim(order.table) else order <- NA
  if(length(class.table)>0) class <- read.delim(class.table) else class <- NA
  if(length(phylum.table)>0) phylum <- read.delim(phylum.table) else phylum <- NA

  
  taxa <- data.frame(cbind(species,genus,family,order,class,phylum))
  taxa <- taxa[,names(colSums(taxa)[!is.na(colSums(taxa))])]

    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))])
   
 taxa <- taxa[, colSums(taxa/ meta$ReadCount > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
       
    if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {
    
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
    dataset[,group] <- as.character(dataset[,group])
    dataset <- dataset[!is.na(dataset[,group]),]
    dataset[,group][dataset[,group]==""] <- "nogroup"
    for(i in unique(dataset[,group])) if(table(dataset[,group])[i] < 3) dataset[,group][dataset[,group]==i] <- "toofewcases"
    dataset <- dataset[dataset[,group]!="toofewcases",]
    dataset[, group] <- dataset[, group][drop = T]
    if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
   dataset$G <- as.factor(dataset[,group])
    if (compare.to != "0") dataset[, group][dataset[, group] == "0" & !is.na(dataset[, group])] <- "group0"
    dataset[, group][dataset[, group] == compare.to & !is.na(dataset[, group])] <- "0"
    dataset[, group] <- as.factor(dataset[, group])
    
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    group_test <- data.frame(array(dim = c(length(names(taxa)),length(c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"),
                                              paste(levels(dataset[,group])[-1],"p", sep = "_"))))))
    rownames(group_test) <- names(taxa)
    names(group_test) <- c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"), paste(levels(dataset[,group])[-1],"p", sep = "_"))
   group_test$taxon <- rownames(group_test)
  group_test_nonzero <- group_test
    model <- list()
    
    if (length(subject.ID) != 0) {
        dataset$ID <- as.factor(dataset[, subject.ID])
        modeldata <- na.omit(dataset[, c(names(taxa), group,"G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], "ID", 
                  "ReadCount")[c(names(taxa), group, "G",confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], "ID", "ReadCount") != ""]])
        modeldata[, group] <- modeldata[, group][drop=T]
       
        for (i in names(taxa)) {
          model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata, 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
        if(length(model[[i]])>0){
          if ((model[[i]]$alpha)<5){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]], type="response"))  
          tmp$group <- modeldata[,group]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
        for(j in levels(modeldata[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
        }
           for (i in rownames(group_test)[is.na(group_test[,2])]) {
           model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]], type="response"))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(modeldata[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)]
            }}}}}}}
if(nonzero){     
for (i in names(taxa)) {
       if(min(table(modeldata[modeldata[,i]>0,group]))>2){
              model[[paste(i,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                " +", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata[modeldata[,i]>0,],
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
            if(length(model[[paste(i,"nonzero")]])>0){
          if ((model[[paste(i,"nonzero")]]$alpha)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$group <- modeldata[modeldata[,i]>0,group]
          tmp$resdev <- abs(tmp$res)
          if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
      for(j in levels(modeldata[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      }}
    
   for (i in rownames(group_test_nonzero)[is.na(group_test_nonzero[,2])]) {
     if(min(table(modeldata[modeldata[,i]>0,group]))>2){
           model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,   data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
           if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]],type="response"))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[modeldata[,i]>0,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
             for(j in levels(modeldata[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)]
            }}}}}}
           }}
}      
    } else {
      modeldata <- na.omit(dataset[, c(names(taxa), group, "G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5],  
                   "ReadCount")[c(names(taxa), group,"G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], "ReadCount") != ""]])
       modeldata[, group] <- modeldata[, group][drop=T]
    
        for (i in names(taxa)) {
           model[[i]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata, control = glm.control(maxit = 500)),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          if ((model[[i]]$deviance/model[[i]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), fitted = predict(model[[i]], type="response"))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(modeldata[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
        }
           for (i in rownames(group_test)[is.na(group_test[, 2])]) {
           model[[i]] <- tryCatch(lm(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),  data = modeldata),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), fitted = predict(model[[i]], type="response"))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(modeldata[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
if(nonzero){          
      for (i in names(taxa)) {
       if(min(table(modeldata[modeldata[,i]>0,group]))>2){
          model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), data = modeldata[modeldata[,i]>0,],
                 control = glm.control(maxit = 500)),
                error = function(e) NULL)
          if(length(model[[paste(i,"nonzero")]])>0){
          if ((model[[paste(i,"nonzero")]]$deviance/model[[paste(i,"nonzero")]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$group <- modeldata[modeldata[,i]>0,group]
          tmp$resdev <- abs(tmp$res)
          if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
             for(j in levels(modeldata[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      }}
    
           for (i in rownames(group_test_nonzero)[is.na(group_test_nonzero[,2])]) {
           if(min(table(modeldata[modeldata[,i]>0,group]))>2){
           model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),  data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
           if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), fitted = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[modeldata[,i]>0,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
          for(j in levels(modeldata[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      }  
  }    
    }
     if(nrow(na.omit(group_test))>0){
       
    for (k in names(group_test)[grepl(pattern="_p$",x=names(group_test))]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,k], "fdr")
    tmp <- group_test
    tmp[, -1][is.na(tmp[, -1])]<-1 
    sig <-    as.character(rownames(tmp)[sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T) < p.cutoff])
    group_test <- na.omit(group_test)
    for (k in names(group_test)[grepl(pattern="_estimate$",x=names(group_test))]) group_test[, paste(k, "FC", sep = "_")] <- exp(group_test[,k])
    
    if(nonzero){
    for (k in names(group_test_nonzero)[grepl(pattern="_p$",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test_nonzero[,k], "fdr")
    tmp2 <- group_test_nonzero
    tmp2[, -1][is.na(tmp2[, -1])]<-1 
    sigNonzero <-   as.character(rownames(tmp2)[sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < p.cutoff & sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T)])
    group_test_nonzero <- na.omit(group_test_nonzero)
    for (k in names(group_test_nonzero)[grepl(pattern="_estimate$",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FC", sep = "_")] <- exp(group_test_nonzero[,k])
    }  
    
    if(length(confounders[confounders!=""])>0){
   resids <- modeldata
  
    library(nlme)
    for(i in rownames(group_test)){
     resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), 
                                   type="pearson"),
                               error = function(e) NA)
   
    }  
    for(i in levels(resids$G)){
      for(j in rownames(group_test)){
        group_test[j,paste("DevianceFromExcpected",i,sep="_")] <-  mean(resids[resids$G==i,j],na.rm=T)
      }
    }
  detach(package:nlme)
  
    if(nonzero){ 
      residsNonzero <- modeldata
   for(i in rownames(group_test_nonzero)){
  
      residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
      residsNonzero[modeldata[,i]==0,i] <-NA
   }
   for(i in levels(resids$G)){
      for(j in rownames(group_test_nonzero)){
      group_test_nonzero[j,paste("DevianceFromExcpected",i,sep="_")] <-  mean(residsNonzero[resids$G==i,j],na.rm=T)
      }
   }
     names(group_test_nonzero)[-1] <-  paste("nonzero",names(group_test_nonzero)[-1],sep="_")
    group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)
    rownames(group_test)<-group_test$taxon
    
resids2 <- data.frame(cbind(resids[,setdiff(colnames(resids),colnames(taxa))],
                            resids[,setdiff(sig,sigNonzero)],
                            residsNonzero[,sigNonzero]))
colnames(resids2)<-c(setdiff(colnames(resids),colnames(taxa)),setdiff(sig,sigNonzero),sigNonzero)

 } else {
   sigNonzero <- NULL  
   resids2 <- data.frame(cbind(resids[,setdiff(colnames(resids),colnames(taxa))],resids[,sig]))
   colnames(resids2)<-c(setdiff(colnames(resids),colnames(taxa)),sig)
 }

    dataset2 <- data.frame(resids2)
    dataset2[,group] <- dataset2$G
    
    } else {    
      for(i in levels(modeldata$G)){
      for(j in group_test$taxon){
        group_test[j,paste("Mean",i,sep="_")] <-   mean(modeldata[modeldata$G==i,j]/modeldata$ReadCount[modeldata$G==i],na.rm=T)
      }}
      for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
      group_test[,paste("FoldChange",i,sep="_")] <- group_test[,paste("Mean",i,sep="_")] / group_test[,paste("Mean",compare.to,sep="_")]
    }
      
      if(nonzero){
      for(i in levels(modeldata$G)){
      for(j in group_test_nonzero$taxon){
        group_test_nonzero[j,paste("Mean",i,sep="_")] <-   mean(modeldata[modeldata$G==i&modeldata[,j]>0,j]/modeldata$ReadCount[modeldata$G==i&modeldata[,j]>0],na.rm=T)
      }}
      
    for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
      group_test_nonzero[,paste("FoldChange",i,sep="_")] <- group_test_nonzero[,paste("Mean",i,sep="_")] / group_test_nonzero[,paste("Mean",compare.to,sep="_")]
    }
      names(group_test_nonzero)[-1] <- paste("nonzero",names(group_test_nonzero)[-1], sep="_") 
      group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)
      } else  sigNonzero <- NULL  
      names(group_test)[grepl(pattern="group0",x=names(group_test))] <-gsub(pattern="group0",replacement="0",x=names(group_test)[grepl(pattern="group0",x=names(group_test))])
      
      
    dataset2 <- modeldata
    dataset2[,names(taxa)]<-log(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
    if(nonzero) for(i in setdiff(sigNonzero,sig)) dataset2[modeldata[,i]==0,i] <- NA
    }
    
    write.table(group_test, paste("GroupTest_", group, compare.to, "_", select.by, select, ".txt", sep = ""), 
        quote = F, row.names = F, sep = "\t")

#------------
   sig <- unique(c(sig,sigNonzero))
   if(length(sig)>0){    
    
library(metacoder)
    
if(length(species.table)!=0){
  sp <- colnames(species)
} else{
 if(length(genus.table)!=0){
   sp <- colnames(genus)
 } else{
  if(length(family.table)!=0){
    sp <- colnames(family)
  } else{
    if(length(order.table)!=0){
      sp <- colnames(order)
    } else{
      if(length(class.table)!=0){
        sp <- colnames(class)
      } else sp <- colnames(phylum)
    }
  }
}
}

for(h in unique(dataset[,"G"])[unique(dataset[,"G"])!=compare.to & !is.na(unique(dataset[,"G"]))]){ 

if(nonzero){  
a <- group_test[!is.na(group_test[,paste(h,"p",sep="_")])&!is.na(group_test[,paste("nonzero",h,"p",sep="_")])&group_test[,paste(h,"p",sep="_")]==group_test[,paste("nonzero",h,"p",sep="_")]|
               !is.na(group_test[,paste(h,"p",sep="_")])&!is.na(group_test[,paste("nonzero",h,"p",sep="_")])&group_test[,paste(h,"p",sep="_")]<group_test[,paste("nonzero",h,"p",sep="_")]|
               !is.na(group_test[,paste(h,"p",sep="_")])&is.na(group_test[,paste("nonzero",h,"p",sep="_")]),c("taxon",paste(h,"estimate",sep="_"),paste(h,"p",sep="_"))]
b <-  group_test[!is.na(group_test[,paste(h,"p",sep="_")])&!is.na(group_test[,paste("nonzero",h,"p",sep="_")])&group_test[,paste(h,"p",sep="_")]>group_test[,paste("nonzero",h,"p",sep="_")]|
               is.na(group_test[,paste(h,"p",sep="_")])&!is.na(group_test[,paste("nonzero",h,"p",sep="_")]),c("taxon",paste("nonzero",h,"estimate",sep="_"),paste("nonzero",h,"p",sep="_"))]
names(b) <- names(a)
group_test_summary <- rbind(a,b)
} else{
group_test_summary <- group_test[!is.na(group_test[,paste(h,"p",sep="_")]),c("taxon",paste(h,"estimate",sep="_"),paste(h,"p",sep="_"))]
}
  
group_test_summary$name <- sapply(group_test_summary$taxon, function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])
rownames(group_test_summary)<-  group_test_summary$taxon
group_test_summary[,paste(h,"estimate",sep="_")][group_test_summary[,paste(h,"p",sep="_")] > p.cutoff] <- 0
group_test_summary2 <- group_test_summary[intersect(sp,rownames(group_test_summary)),]
  
seqs1 <- list()
for(i in rownames(group_test_summary2)){
  if(!is.na(group_test_summary2[i,paste(h,"p",sep="_")]))  if(group_test_summary2[i,paste(h,"p",sep="_")]<p.cutoff)  seqs1[[group_test_summary2[i,"taxon"]]]<- group_test_summary2[i,paste(h,"estimate",sep="_")]
}

for(i in rownames(group_test_summary2)){
  if(!is.na(group_test_summary2[i,paste(h,"p",sep="_")]))  if(group_test_summary2[i,paste(h,"p",sep="_")]>p.cutoff) seqs1[[group_test_summary2[i,"taxon"]]]<- 0
}

abu <- list()
for(i in names(seqs1)) abu[[i]] <-  mean(reltaxa[,i])

newseqs2 <- paste("XX", names(seqs1),seqs1, abu, sep=" ")

tmp2 <- extract_taxonomy(input=newseqs2,regex = "^(.*)\\ (.*)\\ (.*)\\ (.*)",
                         key=c(id = "obs_info","class","taxon_info","taxon_info"),class_sep = "_")

for(k in rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)&tmp2$taxon_data$name!="NA"]){
if(tmp2$taxon_data$name[as.numeric(k)]%in%group_test_summary$name) tmp2$taxon_data$taxon_info_1[as.numeric(k)] <- group_test_summary[group_test_summary$name==tmp2$taxon_data$name[as.numeric(k)],paste(h,"estimate",sep="_")][1]
}
for(k in rev(rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)])){
 tmp2$taxon_data$taxon_info_1[as.numeric(k)] <-  mean(as.numeric(tmp2$taxon_data$taxon_info_1[!is.na(tmp2$taxon_data$supertaxon_ids)&tmp2$taxon_data$supertaxon_ids==tmp2$taxon_data$taxon_ids[as.numeric(k)]]),na.rm=T)
}


treltaxa <- as.data.frame(t(reltaxa))
treltaxa$name <-  sapply(rownames(treltaxa), function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])

for(k in rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_2)]){
if(tmp2$taxon_data$name[as.numeric(k)]%in%treltaxa$name) tmp2$taxon_data$taxon_info_2[as.numeric(k)] <- rowMeans(treltaxa[treltaxa$name==tmp2$taxon_data$name[as.numeric(k)],-ncol(treltaxa)][1,])
}


for(k in rev(rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_2)])){
 tmp2$taxon_data$taxon_info_2[as.numeric(k)] <-  sum(as.numeric(tmp2$taxon_data$taxon_info_2[!is.na(tmp2$taxon_data$supertaxon_ids)&tmp2$taxon_data$supertaxon_ids==tmp2$taxon_data$taxon_ids[as.numeric(k)]]),na.rm=T)
}

heat_tree(tmp2, node_size = as.numeric(taxon_info_2)*100, 
                   node_label = ifelse(name == "NA", NA, name),
                   node_color = as.numeric(taxon_info_1),
                   node_color_range=c("royalblue","cornflowerblue","gray90","hotpink","red"),
                   node_color_interval = c(-max(abs(as.numeric(tmp2$taxon_data$taxon_info_1))), max(abs(as.numeric(tmp2$taxon_data$taxon_info_1)))),
                   edge_color_interval = c(-max(abs(as.numeric(tmp2$taxon_data$taxon_info_1))), max(abs(as.numeric(tmp2$taxon_data$taxon_info_1)))),
                   node_color_axis_label = paste("Group",h,"compared to",compare.to),
                   node_size_axis_label = "Average relative abundance (%)",
          node_label_size_range=c(0.005,0.014),
          node_label_size_trans="ln",
          #node_label_size=0.03,
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          overlap_avoidance = 0.5, node_label_max=150, make_legend=T, 
          node_size_trans="linear",
          node_size_range=c(0.01,0.05),
          node_color_trans="linear",
          title=paste("Differences between groups", h, "and", compare.to), title_size=0.03,
          output_file=paste(group,h,"vs", compare.to, "_", select.by, select,"_HeatTree.pdf",sep=""))


}

#------------
    

   sig <-  na.omit(sig[apply(dataset2[,c("ReadCount",sig)],MARGIN=2,FUN=sum,na.rm=T)!=0])
    if (length(sig) > 0) {
        if (pdf) {
            pdf(paste(group, compare.to, "_", select.by, select, "_Barplot.pdf", sep = ""))
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.5, 0), mar = c(2, 2, 1, 0.5), tck = -0.05, 
                cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
              if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
             boxplot(dataset2[, i] ~ dataset2[, "G"], yaxt="n", ylab="",main =  lab, xlab="", las = label.direction, 
              boxfill= c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                    "G"]))], outpch = 21, outbg = c("skyblue", "yellowgreen", 
                    "pink", "turquoise2", "plum", "darkorange", "lightyellow", 
                    "gray")[1:length(unique(dataset2[, "G"]))], outcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, "G"]))], 
                  boxcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    "G"]))], medcol = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, "G"]))], whiskcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, "G"]))], 
                  staplecol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    "G"]))])
              if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=exp(yaxis))
              } else  axis(side=2)
         }
           dev.off()
            pdf(paste(group, compare.to, "_", select.by, select, "_Beanplot.pdf", sep = ""))
           par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 1), mgp = c(2, 0.5, 0), 
               mar = c(2, 2, 1, 0.5), tck = -0.05, cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
               for(i in sig[order(sig)]){
                if(length(confounders[confounders!=""])==0) {
                lab <- gsub(i,pattern="_NA",replacement = ".")  
                lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- round(log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)))
                } else{
                  lab <- gsub(i,pattern="_NA",replacement = ".")  
                  lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
              tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], xlab="", las = label.direction,yaxt="n",
                  ll = 0.1, ylab = "", main=lab, beanlines="median", col = list(c("skyblue", 
                    "royalblue", "royalblue", "royalblue"), c("yellowgreen", 
                    "olivedrab4", "olivedrab4", "olivedrab4"), c("pink", "red", 
                    "red", "red"), c("turquoise2", "turquoise4", "turquoise4", 
                    "turquoise4"), c("plum", "purple", "purple", "purple"), 
                    c("darkorange", "darkorange3", "darkorange3", "darkorange3"), 
                    c("lightyellow", "lightyellow4", "lightyellow4", "lightyellow4"), 
                    c("gray", "black", "black", "black")), border = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")), error = function(e) NULL)
                 if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=signif(exp(yaxis),digits=1))
              } else  axis(side=2)
        }    
            dev.off()
        }
  quartz()
  par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 1), 
      mgp = c(2, 0.5, 0), mar = c(2, 2, 1, 0.5), tck = -0.05, 
                cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
              if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
             boxplot(dataset2[, i] ~ dataset2[, "G"], yaxt="n", ylab="",main =  lab, xlab="", las = label.direction, 
              boxfill= c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                    "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                    "G"]))], outpch = 21, outbg = c("skyblue", "yellowgreen", 
                    "pink", "turquoise2", "plum", "darkorange", "lightyellow", 
                    "gray")[1:length(unique(dataset2[, "G"]))], outcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, "G"]))], 
                  boxcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    "G"]))], medcol = c("royalblue", "olivedrab4", "red", 
                    "turquoise4", "purple", "darkorange3", "lightyellow4", 
                    "black")[1:length(unique(dataset2[, "G"]))], whiskcol = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")[1:length(unique(dataset2[, "G"]))], 
                  staplecol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                    "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                    "G"]))])
              if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=exp(yaxis))
              } else  axis(side=2)
         }
    quartz()
                  par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 1), mgp = c(2, 0.5, 0), 
               mar = c(2, 2, 1, 0.5), tck = -0.05, cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
               for(i in sig[order(sig)]){
                if(length(confounders[confounders!=""])==0) {
                lab <- gsub(i,pattern="_NA",replacement = ".")  
                lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- round(log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)))
                } else{
                  lab <- gsub(i,pattern="_NA",replacement = ".")  
                  lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
              tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], xlab="", las = label.direction,yaxt="n",
                  ll = 0.1, ylab = "", main=lab, beanlines="median", col = list(c("skyblue", 
                    "royalblue", "royalblue", "royalblue"), c("yellowgreen", 
                    "olivedrab4", "olivedrab4", "olivedrab4"), c("pink", "red", 
                    "red", "red"), c("turquoise2", "turquoise4", "turquoise4", 
                    "turquoise4"), c("plum", "purple", "purple", "purple"), 
                    c("darkorange", "darkorange3", "darkorange3", "darkorange3"), 
                    c("lightyellow", "lightyellow4", "lightyellow4", "lightyellow4"), 
                    c("gray", "black", "black", "black")), border = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")), error = function(e) NULL)
                 if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=signif(exp(yaxis),digits=1))
              } else  axis(side=2)
        } 
   }
    } 
   if(keep.result)   return(group_test)
    }
    }
    palette("default")
} 
