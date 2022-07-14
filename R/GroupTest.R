GroupTest <- function(species.table=NULL, genus.table=NULL, family.table=NULL, order.table=NULL, class.table=NULL, phylum.table=NULL, 
                      meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
                      subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, 
                      select.by = NULL, select = NULL, pdf = F,  min.prevalence = 0.1, min.abundance = 0, 
                      label.direction = 1, keep.result = F, nonzero = T, relative = T) {

  

  
  if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}
library(nlme)  
  meta <- read.delim(meta)
  if(!relative) meta[,"ReadCount"]<-1
  
  if(length(species.table)>0) species <- read.delim(species.table) else species <- NA
  if(length(genus.table)>0) genus <- read.delim(genus.table) else genus <- NA
  if(length(family.table)>0) family <- read.delim(family.table) else family <- NA
  if(length(order.table)>0) order <- read.delim(order.table) else order <- NA
  if(length(class.table)>0) class <- read.delim(class.table) else class <- NA
  if(length(phylum.table)>0) phylum <- read.delim(phylum.table) else phylum <- NA

   if(min.prevalence<0.3) min.prevalence <- 0.3
  taxa <- data.frame(cbind(species,genus,family,order,class,phylum))
  
   if (length(select.by) != 0) {
        meta$selection <- meta[, select.by]
        taxa <- taxa[meta$selection == select, ]
        meta <- meta[meta$selection == select, ]
    }
  
  taxa <- taxa[,names(colSums(taxa)[!is.na(colSums(taxa))])]
  
    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))])

for(i in grep(pattern="_IncertaeSedis",x=colnames(taxa))){
  colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_IncertaeSedis",fixed=T,
replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
}
    
for(i in grep(pattern="_incertaesedis",x=colnames(taxa))){
  colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_incertaesedis",fixed=T,
replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
}
    
for(i in grep(pattern="_uncultured",x=colnames(taxa))){
  colnames(taxa)[i] <- gsub(x=colnames(taxa)[i] ,pattern="_uncultured",fixed=T,
replacement = paste("_",strsplit(x=colnames(taxa)[i],split = "_")[[1]][4],strsplit(x=colnames(taxa)[i],split = "_")[[1]][5],sep=""))
}
    

    taxa <- taxa[meta$ReadCount > readcount.cutoff, ]
    meta <- meta[meta$ReadCount > readcount.cutoff, ]
    
   
    
   reltaxa <- (taxa)/meta$ReadCount
    for (i in names(reltaxa)) {
        for (j in 1:nrow(taxa)) {
            reltaxa[j, i][reltaxa[j, i] > (mean(reltaxa[, i]) + outlier.cutoff*sd(reltaxa[, i]))] <- mean(reltaxa[, i]) + 1.05*outlier.cutoff*sd(reltaxa[, i])
        }
    }

      
    taxa <- round(reltaxa * meta$ReadCount)
    taxa[taxa<0]<-0
    
    sel <- colSums(taxa > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)
    
    if(length(sel[sel==T])==1) {
     taxa <- as.data.frame(taxa[,sel])
      colnames(taxa)[1]<-names(sel[sel==T])
    } else {
      taxa <- taxa[, colSums(taxa > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
    }
    
  
       
    if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {   
      
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
    
     for(i in c(1:ncol(dataset))[-c(1:ncol(meta),ncol(dataset))]){
      if(min(dataset[,i]) == max(dataset[,i]))dataset[,i] <- NA
     }
    
    
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    group_test <- data.frame(array(dim = c(length(names(taxa)),length(c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"),
                                              paste(levels(dataset[,group])[-1],"p", sep = "_"))))))
    rownames(group_test) <- names(taxa)
    names(group_test) <- c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"), 
                           paste(levels(dataset[,group])[-1],"p", sep = "_"))
   group_test$taxon <- rownames(group_test)
   group_test_nonzero <- group_test
    model <- list()
    updatedmodel <- list()
    
    if (length(subject.ID) != 0) {
        dataset$ID <- as.factor(dataset[, subject.ID])
        modeldata <- na.omit(dataset[, c(names(taxa), group,"G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], "ID", 
                  "ReadCount")[c(names(taxa), group, "G",confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], "ID", "ReadCount") != ""]])
        modeldata[, group] <- modeldata[, group][drop=T]
       
       resids <- modeldata
       resids[,names(taxa)]<-NA
       modeldata2 <- modeldata
       modeldata2[,names(taxa)]<- modeldata2[,names(taxa)]+1
       
        for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
           model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata, 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
         
            if(length(model[[i]])==0|max(model[[i]]$b)==0){
               model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata, 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
            }
           
               if(length(model[[i]])==0|max(model[[i]]$b)==0){
               model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "poisson", 
                data = modeldata, 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
            }
             
          if(length(model[[i]])==0|max(model[[i]]$b)==0 | model[[i]]$deviance/model[[i]]$df.residual > 1.4){
            if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.255){
           model[[i]] <- tryCatch(nlme::lme(as.formula(paste("((", i, ")/ReadCount)~", confounders[1], "+", 
                                                             confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata2),
                error = function(e) NULL)
           
            } else{
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                                                             confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata2),
                error = function(e) NULL)
           }}  
     
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
          
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          
          if(length(summary(model[[i]])$tTable)==0){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$coef),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)] }   
          } else {
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)] }
          }
         } else {
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
           if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
             model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varExp(),random = ~1 | ID,  
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
             
           } else {
              model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varExp(),random = ~1 | ID,  
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           }
          } else {
          if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
            model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varIdent(form=as.formula(paste("~1|",group))),
                random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
          } else {
               model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varIdent(form=as.formula(paste("~1|",group))),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)  
              }
          }    
          if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)]  }
          } else {
             if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){  
            model[[i]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, random = ~1 | ID, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
             } else {
             model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, random = ~1 | ID, 
                weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)  
             }
            if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
        group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,5)]  }  
          }
            }
          }
          }
         } 
          
              if(length(model[[i]])>0){
     group_test[i, "model"] <- strsplit(as.character(summary(model[[i]])$call),split="(",fixed=T)[[1]][1] 
     
     if(group_test[i, "model"]== "lme.formula"){
    if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)>0.24999){
    group_test[i, "model"] <- paste("log",group_test[i, "model"])
     } 
     }
     
      for(j in levels(modeldata[,group])[-1]){

       if(group_test[i, "model"]== "glmmADMB::glmmadmb"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
       } 
        
        if(group_test[i, "model"] == "log lme.formula"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
        }
        
      if(group_test[i, "model"] == "lme.formula"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- (summary(model[[i]])$tTable[1,1]+summary(model[[i]])$tTable[paste(group,j,sep=""),1])/summary(model[[i]])$tTable[1,1]
      }
    }
          
    if(length(confounders[confounders!=""])>0){
            resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), 
                                   type="pearson"),
                               error = function(e) NA)
     }
         
        }  
}
         if(nonzero){       
  residsNonzero <- modeldata
    confounders2 <- confounders
   for(h in seq_along(confounders[confounders!=""])) if(min(table(modeldata[modeldata[,i]>0,confounders[h]]))==0) confounders2[h] <- ""

     for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
        if(min(table(modeldata[modeldata[,i]>0,group]))>2){        
          model[[paste(i,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata[modeldata[,i]>0,], 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
          
            if(length(model[[paste(i,"nonzero")]])==0){
               model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata[modeldata[,i]>0,], 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
            }
           
              if(length(model[[paste(i,"nonzero")]])==0){
               model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "poisson", 
                data = modeldata[modeldata[,i]>0,], 
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
              }
          
             if(length(model[[paste(i,"nonzero")]])>0){
                if(model[[paste(i,"nonzero")]]$deviance/model[[paste(i,"nonzero")]]$df.residual > 1.4) model[[paste(i,"nonzero")]] <- NULL
              }
          
          if(length(model[[paste(i,"nonzero")]])==0){
            if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){
           model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("((", i, ")/ReadCount)~", confounders[1], "+", 
                                                             confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
            } else {
              model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                                                             confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
           }}  
     
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
          
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          
          if(length(summary(model[[paste(i,"nonzero")]])$tTable)==0){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$coef),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)] }   
          } else {
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)] }
          }
         } else {
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
            if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){
             model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           } else {
              model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(),random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           }
          } else {
          if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){  
            model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], 
                weights = nlme::varIdent(form=as.formula(paste("~1|",group))),random = ~1 | ID,  
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
          } else {
               model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form=as.formula(paste("~1|",group))),
                random = ~1 | ID,  control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)  
              }
          }    
          if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)]  }
          } else {
              if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<0.25){ 
            model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], random = ~1 | ID, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
             } else {
             model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], random = ~1 | ID, 
                weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form=as.formula(paste("~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)  
             }
            if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
        group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,5)]  }  
          }
            }
          }
          }
         }
          
 if(length(model[[paste(i,"nonzero")]])>0){
     group_test_nonzero[i, "model"] <- strsplit(as.character(summary(model[[paste(i,"nonzero")]])$call),split="(",fixed=T)[[1]][1] 
     
     if(group_test_nonzero[i, "model"]== "lme.formula"){
    if(abs(mean(modeldata[,i]/modeldata$ReadCount)-median(modeldata[,i]/modeldata$ReadCount))/median(modeldata[,i]/modeldata$ReadCount)>0.24999){
    group_test_nonzero[i, "model"] <- paste("log",group_test_nonzero[i, "model"])
     } 
     }
     
      for(j in levels(modeldata[,group])[-1]){
       if(group_test_nonzero[i, "model"]== "glmmADMB::glmmadmb"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
       } 
        
        if(group_test_nonzero[i, "model"] == "log lme.formula"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
        }
        
      if(group_test_nonzero[i, "model"] == "lme.formula"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- (summary(model[[paste(i,"nonzero")]])$tTable[1,1]+summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1])/summary(model[[paste(i,"nonzero")]])$tTable[1,1]
      }
    }
        }  
          
    if(length(confounders[confounders!=""])>0){
      residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), 
                                   type="pearson"),
                               error = function(e) NA)
     }
         
        }  

  }
         }
      } else {
      modeldata <- na.omit(dataset[, c(names(taxa), group, "G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5],  
                   "ReadCount")[c(names(taxa), group,"G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], "ReadCount") != ""]])
       modeldata[, group] <- modeldata[, group][drop=T]
       resids <- modeldata
       modeldata2 <- modeldata
       modeldata2[,names(taxa)]<- modeldata2[,names(taxa)]+1
       
        for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
           model[[i]] <- tryCatch(MASS::glm.nb(formula(paste("round(", i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata, control = glm.control(maxit = 1000,epsilon=1e-12)),
                error = function(e) NULL)
          if(length(model[[i]])==0){
              model[[i]] <- tryCatch(MASS::glm.nb(formula(paste("round(", i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata, control = glm.control(maxit = 1000,epsilon=1e-12)),
                error = function(e) NULL)
          }
            if(length(model[[i]])==0){
              model[[i]] <- tryCatch(glm(formula(paste("round(", i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), family=poisson,
                data = modeldata,control = glm.control(maxit = 1000,epsilon=1e-12)),
                error = function(e) NULL)
           }
           if(length(model[[i]])==0 | model[[i]]$deviance/model[[i]]$df.residual > 1.4){
             if(abs(mean(modeldata[,i]/modeldata$ReadCount)-median(modeldata[,i]/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){
             model[[i]] <- tryCatch(lm(formula(paste("(", i, ")/ReadCount~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata),
                error = function(e) NULL)   
             } else {
               model[[i]] <- tryCatch(lm(formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata),
                error = function(e) NULL)   
             }
           }  
     if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
     
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$coef),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)] }   
           
         } else {
           
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
           if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){   
            model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           } else {
             model[[i]] <- tryCatch(nlme::gls(formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varExp(), control = nlme::glsControl(maxIter=5000)),
                error = function(e) NULL) 
           }
             
          } else {
           if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){    
             model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varIdent(form = formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
            } else(
              model[[i]] <- tryCatch(nlme::gls(formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varIdent(form = formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
            )
          }    
            if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,4)]}
          } else {
         if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)<0.25){    
          model[[i]] <- tryCatch(nlme::gls(formula(paste("(", i, "/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, 
                weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = formula(paste(" ~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           } else{
              model[[i]] <- tryCatch(nlme::gls(formula(paste("log((", i, ")/ReadCount)~", confounders[1], "+", 
                   confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)), 
                data = modeldata2, weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = formula(paste(" ~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           }
          
          if(length( model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), fitted = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          for(j in intersect(levels(modeldata[,group]),gsub(rownames(summary(model[[i]])$tTable),pattern=group,replacement=""))){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$tTable[paste(group,j,sep=""),c(1,4)]
           } 
         }  
          }
          }
            }
         }
    if(length(model[[i]])>0){
     group_test[i, "model"] <- strsplit(as.character(summary(model[[i]])$call),split="(",fixed=T)[[1]][1] 
     if(group_test[i, "model"]== "lm" | group_test[i, "model"]== "nlme::gls"){
    if(abs(mean((modeldata[,i]+1)/modeldata$ReadCount)-median((modeldata[,i]+1)/modeldata$ReadCount))/median((modeldata[,i]+1)/modeldata$ReadCount)>0.24999){
    group_test[i, "model"] <- paste("log",group_test[i, "model"])
     } 
     }
     
      for(j in levels(modeldata[,group])[-1]){

 if(group_test[i, "model"]== "MASS::glm.nb" |group_test[i, "model"]== "glm"){     
  group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test[i, paste(j,"estimate",sep="_")])
         if(length(confounders[confounders!=""])>0){
 resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                      init.theta = 5, control = glm.control(epsilon = 1e-12, maxit = 2500, trace = FALSE)), type="pearson"),
                               error = function(e) NA)
 if(is.na(sum(resids[,i]))){
resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                      control = glm.control(epsilon = 1e-20, maxit = 1000, trace = FALSE)), type="pearson"),
                               error = function(e) NA)
   }
  if(is.na(sum(resids[,i]))){
resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep="")),
                                       init.theta = 10,control = glm.control(epsilon = 1e-20, maxit = 2500, trace = FALSE)), type="pearson"),
                               error = function(e) NA)
   
 }  }  } 
        
   if(group_test[i, "model"] == "nlme::gls"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- mean(fitted(model[[i]])[modeldata$G==j])/mean(fitted(model[[i]])[modeldata$G==compare.to])#(summary(model[[i]])$tTable[1,1]+summary(model[[i]])$tTable[paste(group,j,sep=""),1])/summary(model[[i]])$tTable[1,1]
         if(length(confounders[confounders!=""])>0){
 resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
         }   } 
        
                if(group_test[i, "model"] == "log nlme::gls"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <- mean(exp(fitted(model[[i]])[modeldata$G==j]))/mean(exp(fitted(model[[i]])[modeldata$G==compare.to]))#exp(summary(model[[i]])$tTable[1,1]+(summary(model[[i]])$tTable[paste(group,j,sep=""),1]))/exp(summary(model[[i]])$tTable[1,1])
         if(length(confounders[confounders!=""])>0){
 resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
       }   } 
        
          if(group_test[i, "model"] == "lm"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <-  mean(fitted(model[[i]])[modeldata$G==j])/mean(fitted(model[[i]])[modeldata$G==compare.to])#(summary(model[[i]])$coef[1,1]+summary(model[[i]])$coef[paste(group,j,sep=""),1])/summary(model[[i]])$coef[1,1]
           if(length(confounders[confounders!=""])>0){
 resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
       }   } 
        
        if(group_test[i, "model"] == "log lm"){     
        group_test[i,paste(j,"estimate_FoldChange",sep="_")] <-    mean(exp(fitted(model[[i]])[modeldata$G==j]))/mean(exp(fitted(model[[i]])[modeldata$G==compare.to]))#(exp(summary(model[[i]])$coef[1,1]+summary(model[[i]])$coef[paste(group,j,sep=""),1]))/exp(summary(model[[i]])$coef[1,1])
             if(length(confounders[confounders!=""])>0){
 resids[,i] <-  tryCatch(resid(update(model[[i]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
             }  } 
        
     
      }
    
        }  
        }
        }        
if(nonzero){       
  residsNonzero <- modeldata
  
      for (i in names(colSums(modeldata[,names(taxa)]>0)[colSums(modeldata[,names(taxa)]>0)>5])) {
    if(min(table(modeldata[modeldata[,i]>0,group]))>3){    
    confounders2 <- confounders
   for(h in seq_along(confounders[confounders!=""])) if(min(table(modeldata[modeldata[,i]>0,confounders[h]]))==0) confounders2[h] <- ""
      
             model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, ")~", confounders2[1], "+", confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata[modeldata[,i]>0,], control = glm.control(maxit = 1000)),
                error = function(e) NULL)
          if(length(model[[paste(i,"nonzero")]])==0){
              model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, ")~", confounders2[1], "+", confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata[modeldata[,i]>0,], control = glm.control(maxit = 1000)),
                error = function(e) NULL)
          }
              if(length(model[[paste(i,"nonzero")]])==0){
              model[[i]] <- tryCatch(glm(formula(paste("round(", i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), family=poisson,
                data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
              }
             
              if(length(model[[paste(i,"nonzero")]])>0){
                if(model[[paste(i,"nonzero")]]$deviance/model[[paste(i,"nonzero")]]$df.residual > 1.4) model[[paste(i,"nonzero")]] <- NULL
              }
             
           if(length(model[[paste(i,"nonzero")]])==0){
             if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){
             model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("(", i, ")/ReadCount~", confounders2[1], "+", confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)   
             } else {
               model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("log((", i, ")/ReadCount)~", confounders2[1], "+", confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)   
             }
           }  
           
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
          
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]!="NaN" & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]!="NaN" &
              summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 &
             anova(lm(res ~ group, data=tmp))$Pr[1]!="NaN" & anova(lm(resdev ~ group, data=tmp))$Pr[1]!="NaN" &
              anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
       
            for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$coef),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)] }   
          
         } else {
           
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
          if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){
            model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           } else {
             model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log((", i, ")/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varExp(), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL) 
           }
             
          } else if(anova(lm(res ~ group, data=tmp))$Pr[1]<0.01 | anova(lm(resdev ~ group, data=tmp))$Pr[1]<0.01){
            if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){  
             model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form = as.formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
            } else(
              model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log((", i, ")/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varIdent(form = as.formula(paste(" ~1|",group))), control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
            )
          }    
            if(length( model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
            for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,4)]}
          } else {
           if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])<1){   
          model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("(", i, "/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], 
                weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = as.formula(paste(" ~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           } else{
              model[[paste(i,"nonzero")]] <- tryCatch(nlme::gls(as.formula(paste("log((", i, ")/ReadCount)~", confounders2[1], "+", 
                   confounders2[2], "+", confounders2[3], 
                "+", confounders2[4], "+", confounders2[5], "+ ",group)), 
                data = modeldata[modeldata[,i]>0,], weights = nlme::varComb(nlme::varExp(),nlme::varIdent(form = as.formula(paste(" ~1|",group)))), 
                control = nlme::glsControl(maxIter=1000)),
                error = function(e) NULL)
           }
          
          if(length( model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), fitted = predict(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[modeldata[,i]>0,group]
           if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.1 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.1 |
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.1 | anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.1){
          for(j in intersect(levels(modeldata[modeldata[,i]>0,group]),gsub(rownames(summary(model[[paste(i,"nonzero")]])$tTable),pattern=group,replacement=""))){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),c(1,4)]
          } 
         }  
          }
          }
            }
         }    
          
              if(length(model[[paste(i,"nonzero")]])>0){
     group_test_nonzero[i, "model"] <- strsplit(as.character(summary(model[[paste(i,"nonzero")]])$call),split="(",fixed=T)[[1]][1] 
     if(group_test_nonzero[i, "model"]== "lm" | group_test_nonzero[i, "model"]== "nlme::gls"){
    if(abs(mean(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])-median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0]))/median(modeldata[modeldata[,i]>0,i]/modeldata$ReadCount[modeldata[,i]>0])>0.24999){
    group_test_nonzero[i, "model"] <- paste("log",group_test_nonzero[i, "model"])
     } 
     }    
     
          for(j in levels(modeldata[,group])[-1]){
            
 if(group_test_nonzero[i, "model"]== "MASS::glm.nb" |group_test_nonzero[i, "model"]== "glm"){     
  group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- exp(group_test_nonzero[i, paste(j,"estimate",sep="_")])
       
    if(length(confounders[confounders!=""])>0){
try( residsNonzero[residsNonzero[,i]>0,i] <-  resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                      init.theta = 5, control = glm.control(epsilon = 1e-12, maxit = 2500, trace = FALSE)), type="pearson"))#,
                               #error = function(e) NA)
 if(is.na(sum(residsNonzero[,i]))){
residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                      control = glm.control(epsilon = 1e-20, maxit = 1000, trace = FALSE)), type="pearson"),
                               error = function(e) NA)
   }
  if(is.na(sum(residsNonzero[,i]))){
residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep="")),
                                       init.theta = 10,control = glm.control(epsilon = 1e-20, maxit = 2500, trace = FALSE)), type="pearson"),
                               error = function(e) NA)
   
 }  }  } 
        
   if(group_test_nonzero[i, "model"] == "nlme::gls"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0])/mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0])#(summary(model[[paste(i,"nonzero")]])$tTable[1,1]+summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1])/summary(model[[paste(i,"nonzero")]])$tTable[1,1]
         if(length(confounders[confounders!=""])>0){
 residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
         }   } 
        
                if(group_test_nonzero[i, "model"] == "log nlme::gls"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <- mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0]))/mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0]))#exp(summary(model[[paste(i,"nonzero")]])$tTable[1,1]+(summary(model[[paste(i,"nonzero")]])$tTable[paste(group,j,sep=""),1]))/exp(summary(model[[paste(i,"nonzero")]])$tTable[1,1])
         if(length(confounders[confounders!=""])>0){
 residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
       }   } 
        
          if(group_test_nonzero[i, "model"] == "lm"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <-  mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0])/mean(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0])#(summary(model[[paste(i,"nonzero")]])$coef[1,1]+summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),1])/summary(model[[paste(i,"nonzero")]])$coef[1,1]
           if(length(confounders[confounders!=""])>0){
 residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
       }   } 
        
        if(group_test_nonzero[i, "model"] == "log lm"){     
        group_test_nonzero[i,paste(j,"estimate_FoldChange",sep="_")] <-    mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==j&modeldata[,i]>0]))/mean(exp(fitted(model[[paste(i,"nonzero")]])[modeldata$G==compare.to&modeldata[,i]>0]))#(exp(summary(model[[paste(i,"nonzero")]])$coef[1,1]+summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),1]))/exp(summary(model[[paste(i,"nonzero")]])$coef[1,1])
             if(length(confounders[confounders!=""])>0){
 residsNonzero[residsNonzero[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), type="pearson"),
                               error = function(e) NA)
             }  } 
        }

     
     
     
    }
          
  if(length(confounders2[confounders2!=""])>0){
    
            residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update( model[[paste(i,"nonzero")]], as.formula(paste(".~. -",group,sep=""))), 
                                   type="pearson"),
                               error = function(e) NA)
           }      
}     
        
}
}
      }
    if(nrow(na.omit(group_test))>0){
       
    for (k in names(group_test)[grepl(pattern="_p$",x=names(group_test))]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,k], "fdr")
    tmp <- group_test
    tmp[, -1][is.na(tmp[, -1])]<-1 
    sig <-    as.character(rownames(tmp)[sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T) < p.cutoff])
 
    
    if(nonzero){
    for (k in names(group_test_nonzero)[grepl(pattern="_p$",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test_nonzero[,k], "fdr")
    tmp2 <- group_test_nonzero
    tmp2[, -1][is.na(tmp2[, -1])]<-1 
    sigNonzero <-   as.character(rownames(tmp2)[sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < p.cutoff & sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T)])
    }  
    
    if(length(confounders[confounders!=""])>0){
  
    if(nonzero){ 
   
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
   if(relative) dataset2[,names(taxa)]<-(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
    if(nonzero) for(i in setdiff(sigNonzero,sig)) dataset2[modeldata[,i]==0,i] <- NA
    }
    
    group_test <- group_test[order(group_test$taxon),]
    
  write.table(group_test, paste("GroupTest_", group, compare.to, "_", select.by, select, ".txt", sep = ""), 
        quote = F, row.names = F, sep = "\t")
    
if(length(confounders[confounders!=""])>0){  
  write.table(resids, paste("GroupTest_residuals_", group, compare.to, "_", select.by, select, ".txt", sep = ""), 
        quote = F, row.names = F, sep = "\t")
}
#------------
   sig <- unique(c(sig,sigNonzero))
   if(length(sig)>0){    
    
#------------
  if (length(sig) > 1) {  sig <-  na.omit(sig[apply(dataset2[,sig],MARGIN=2,FUN=sum,na.rm=T)!=0])}
   
lphy <- length(sig)

if(lphy < 5) {
  ncols <- lphy
  nrows <- 1
} else {
ncols <-  round(sqrt(lphy)) + 1
nrows <-  floor(sqrt(lphy))
}

   if(relative) yaxis <- round((100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))) else yaxis <- round((seq(0,max(dataset[,sig]),max(dataset[,sig])/10)))

        if (pdf) {
            pdf(width=ncols*3,height=nrows*3, paste(group, compare.to, "_", select.by, select, "_Boxplot.pdf", sep = ""))
   par(mfrow = c(nrows, ncols),  mgp = c(2, 0.5, 0), mar = c(6, 2, 1, 0.5), tck = -0.01,  cex.axis = 0.8, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
              if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                } else{
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
           boxplot(dataset2[, i] ~ dataset2[, "G"], ylab="", main =  lab, xlab="", las = label.direction, #yaxt="n", 
              col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
       "#A65628", "#F781BF", "#999999","blue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
       outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
       "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
       }
          dev.off()
          
           pdf(width=ncols*3,height=nrows*3, paste(group, compare.to, "_", select.by, select, "_Beanplot.pdf", sep = ""))
                    par(mfrow = c(nrows, ncols), mgp = c(2, 0.5, 0), 
               mar = c(6, 2, 1, 0.5), tck = -0.01, cex.axis = 0.8, cex.lab = 1,cex.main=0.8)
               for(i in sig[order(sig)]){
                if(length(confounders[confounders!=""])==0) {
                lab <- gsub(i,pattern="_NA",replacement = ".")  
                lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                } else {
                  lab <- gsub(i,pattern="_NA",replacement = ".")  
                  lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
              tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], xlab="", las = label.direction,#yaxt="n",
                  ll = 0.1, ylab = "", main=lab, beanlines="median",   col=list(c('#E41A1C','black','black','black'),
         c('orange','black','black','black'),
          c('#377EB8','black','black','black'),
          c('skyblue','black','black','black'),
          c("#4DAF4A",'black','black','black'),
          c('#984EA3','black','black','black'),
          c('#FFFF33','black','black','black'),
         c('#A65628','black','black','black'),
         c('#F781BF','black','black','black'),
         c('#999999','black','black','black'),
          c('blue','black','black','black'),
          c('firebrick4','black','black','black'),
          c('yellowgreen','black','black','black'),
          c('pink','black','black','black'),
          c('turquoise2','black','black','black'),
          c('plum','black','black','black'),
          c('darkorange','black','black','black'),
          c('lightyellow','black','black','black'),
          c('gray','black','black','black')),
                  border = "black"), error = function(e) NULL)
        } 
                    dev.off()
        }
 
       quartz()
  par(mfrow = c(nrows, ncols), 
      mgp = c(2, 0.5, 0), mar = c(6, 2, 1, 0.5), tck = -0.01, 
                cex.axis = 0.8, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
              if(length(confounders[confounders!=""])==0) {
              lab <- gsub(i,pattern="_NA",replacement = ".")  
              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                } else{
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
           boxplot(dataset2[, i] ~ dataset2[, "G"], ylab="", main =  lab, xlab="", las = label.direction, #yaxt="n", 
              col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
       "#A65628", "#F781BF", "#999999","blue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
       outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33", 
       "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
            #  if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=exp(yaxis)) } else  axis(side=2)
         }
  #  quartz()
  #                par(mfrow = c(nrows, ncols), mgp = c(2, 0.5, 0), 
  #             mar = c(6, 2, 1, 0.5), tck = -0.01, cex.axis = 0.8, cex.lab = 1,cex.main=0.8)
  #             for(i in sig[order(sig)]){
  #              if(length(confounders[confounders!=""])==0) {
  #              lab <- gsub(i,pattern="_NA",replacement = ".")  
  #              lab <- strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])]
  #              if(length(lab)==0) lab <- "Unassigned taxa"
  #              } else {
  #                lab <- gsub(i,pattern="_NA",replacement = ".")  
  #                lab <- paste(strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])],"deviance")
  #                if(lab==" deviance") lab <- "Unassigned taxa deviance"
  #                }
  #            tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, "G"], xlab="", las = label.direction,#yaxt="n",
  #                ll = 0.1, ylab = "", main=lab, beanlines="median",   col=list(c('#E41A1C','black','black','black'),
  #       c('orange','black','black','black'),
  #        c('#377EB8','black','black','black'),
  #        c('skyblue','black','black','black'),
  #        c("#4DAF4A",'black','black','black'),
  #        c('#984EA3','black','black','black'),
  #        c('#FFFF33','black','black','black'),
  #       c('#A65628','black','black','black'),
  #       c('#F781BF','black','black','black'),
  #       c('#999999','black','black','black'),
  #        c('blue','black','black','black'),
  #        c('firebrick4','black','black','black'),
  #        c('yellowgreen','black','black','black'),
  #        c('pink','black','black','black'),
  #        c('turquoise2','black','black','black'),
  #        c('plum','black','black','black'),
  #        c('darkorange','black','black','black'),
  #        c('lightyellow','black','black','black'),
  #        c('gray','black','black','black')),
  #                border = "black"), error = function(e) NULL)
  #            #   if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=signif(exp(yaxis),digits=1)) } else  axis(side=2)
  #      } 
    


     
     
    } 
   if(keep.result)   return(group_test)

    }
    }
    palette("default")
    detach(package:nlme)
} 
