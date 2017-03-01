GroupTest <- function(taxonomic.table, meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
    subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, #zinf.cutoff = 0, 
    select.by = NULL, select = NULL, pdf = F,  min.prevalence = 0, min.abundance = 0, 
    label.direction = 1, quartz = T, keep.result = F) {
 
    taxa <- read.delim(taxonomic.table)
    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
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
    if(compare.to != "0") { dataset[, group][dataset[, group] == "0" & !is.na(dataset[, group])] <- "group0"}
    dataset[, group][dataset[, group] == compare.to & !is.na(dataset[, group])] <- "0"
    dataset[, group] <- as.factor(dataset[, group])
    
    confounders <- c(confounders, rep("", 5 - length(confounders)))
    #group_test <- data.frame(array(dim = c(length(names(taxa)), (length(levels(as.factor(dataset[,group])))))))
    group_test <- data.frame(array(dim = c(length(names(taxa)),length(c("taxon", paste(levels(dataset[,group])[-1],"estimate", sep = "_"),
                                              paste(levels(dataset[,group])[-1],"p", sep = "_"))))))
    rownames(group_test) <- names(taxa)
    #names(group_test) <- levels(as.factor(dataset[, group]))
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
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = fitted.values(model[[i]]))  
          tmp$group <- modeldata[,group]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
          # group_test[i, -1] <- summary(model[[i]])$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4]
             for(j in levels(dataset[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
        }
           for (i in rownames(group_test)[is.na(group_test[,paste(j,"estimate",sep="_")])]) {
           model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,  data = modeldata),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), fitted = fitted.values(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(dataset[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
     
for (i in names(taxa)) {
              model[[paste(i,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ as.factor(", group, 
                ")+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                data = modeldata[modeldata[,i]>0,],
                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)),
                error = function(e) NULL)
            if(length(model[[paste(i,"nonzero")]])>0){
          if ((model[[paste(i,"nonzero")]]$alpha)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), pred = fitted.values(model[[paste(i,"nonzero")]]))  
          tmp$group <- modeldata[modeldata[,i]>0,group]
          tmp$resdev <- abs(tmp$res)
          if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
          # group_test_nonzero[i, -1] <-  summary(model[[paste(i,"nonzero")]])$coef[paste("as.factor(",group,")",names(group_test),sep="")[-1],4]
           for(j in levels(dataset[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      }
    
           for (i in rownames(group_test_nonzero)[is.na(group_test_nonzero[,2])]) {
           model[[i]] <- tryCatch(nlme::lme(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),random = ~1 | ID,   data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), fitted = fitted.values(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[modeldata[,i]>0,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(dataset[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
        
    } else {
      modeldata <- na.omit(dataset[, c(names(taxa), group, "G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5],  
                   "ReadCount")[c(names(taxa), group,"G", confounders[1], confounders[2], 
                  confounders[3], confounders[4], confounders[5], "ReadCount") != ""]])
    
        for (i in names(taxa)) {
           model[[i]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", i, "+1)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group,"+", "offset(log(ReadCount))")), 
                data = modeldata, control = glm.control(maxit = 500)),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          if ((model[[i]]$deviance/model[[i]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), fitted = fitted.values(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
          tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(dataset[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
        }
           for (i in rownames(group_test)[is.na(group_test[,2])]) {
           model[[i]] <- tryCatch(lm(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),  data = modeldata),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), fitted = fitted.values(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(dataset[,group])[-1]){
            group_test[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
         
      for (i in names(taxa)) {
          model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste("round(", 
                i, ")~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ", group, 
                "+", "offset(log(ReadCount))")), data = modeldata[modeldata[,i]>0,],
                 control = glm.control(maxit = 500)),
                error = function(e) NULL)
          if(length(model[[paste(i,"nonzero")]])>0){
          if ((model[[paste(i,"nonzero")]]$deviance/model[[paste(i,"nonzero")]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), pred = fitted.values(model[[paste(i,"nonzero")]]))  
          tmp$group <- modeldata[modeldata[,i]>0,group]
          tmp$resdev <- abs(tmp$res)
          if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
         for(j in levels(dataset[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[paste(i,"nonzero")]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      }
    
           for (i in rownames(group_test_nonzero)[is.na(group_test_nonzero[,2])]) {
           model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("log((", i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                "+", confounders[4], "+", confounders[5], "+ ",group)),  data = modeldata[modeldata[,i]>0,]),
                error = function(e) NULL)
           if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), fitted = fitted.values(model[[paste(i,"nonzero")]]))  
          tmp$resdev <- abs(tmp$res) 
           tmp$group <- modeldata[modeldata[,i]>0,group]
         if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(res ~ group, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.05){
            for(j in levels(dataset[,group])[-1]){
            group_test_nonzero[i, c(paste(j,"estimate",sep="_"),paste(j,"p",sep="_"))] <- summary(model[[i]])$coef[paste(group,j,sep=""),c(1,4)]
            }}}}}}}
      
    }
    
   # names(group_test)[1] <- "taxon"
  #  group_test$taxon <- rownames(group_test)
  #  names(group_test)[-1] <- paste("p",names(group_test)[-1],sep="_")
    for (k in names(group_test)[grepl(pattern="_p$",x=names(group_test))]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,k], "fdr")
     
   # names(group_test_nonzero)[1] <- "taxon"
  #  group_test_nonzero$taxon <- rownames(group_test_nonzero)
  #  names(group_test_nonzero)[-1] <- paste("p",names(group_test_nonzero)[-1],sep="_")
    for (k in names(group_test_nonzero)[grepl(pattern="_p$",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test_nonzero[,k], "fdr")
    
    tmp <- group_test
    tmp[, -1][is.na(tmp[, -1])]<-1 
    tmp2 <- group_test_nonzero
    tmp2[, -1][is.na(tmp2[, -1])]<-1 
    
    sig <-    as.character(rownames(tmp)[sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T) < p.cutoff])
    sigNonzero <-   as.character(rownames(tmp2)[sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < p.cutoff & sapply(data.frame(t(tmp2[, grepl(pattern="_p$",x=names(group_test_nonzero))])), min,na.rm=T) < sapply(data.frame(t(tmp[, grepl(pattern="_p$",x=names(group_test))])), min,na.rm=T)])
    group_test <- na.omit(group_test)
    group_test_nonzero <- na.omit(group_test_nonzero)
    
    if(nrow(group_test)>0){
    #rownames(group_test)<-group_test$taxon
    #rownames(group_test_nonzero)<-group_test_nonzero$taxon
    
    for (k in names(group_test)[grepl(pattern="_estimate$",x=names(group_test))]) group_test[, paste(k, "FC", sep = "_")] <- exp(group_test[,k])
     for (k in names(group_test_nonzero)[grepl(pattern="_estimate$",x=names(group_test_nonzero))]) group_test_nonzero[, paste(k, "FC", sep = "_")] <- exp(group_test_nonzero[,k])
      
    
    if(length(confounders[confounders!=""])>0){
   resids <- modeldata
   residsNonzero <- modeldata
   
   for(i in rownames(group_test_nonzero)){
      residsNonzero[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]], 
                              as.formula(paste(".~. -",group,sep="")))),
                               error = function(e) NA)
      residsNonzero[modeldata[,i]==0,i] <-NA
    }
    for(i in rownames(group_test)){
     resids[,i] <-  tryCatch(resid(update(model[[i]], 
                              as.formula(paste(".~. -",group,sep="")))),
                               error = function(e) NA)
    }
   # resids[,setdiff(rownames(group_test),colnames(resids))] <- NA
  #resids[,setdiff(rownames(group_test_nonzero),colnames(residsNonzero))] <- NA
    
    for(i in levels(resids$G)){
      for(j in rownames(group_test)){
        group_test[j,paste("DevianceFromExcpected",i,sep="_")] <-  mean(resids[resids$G==i,j],na.rm=T)
      }
    }

    for(i in levels(resids$G)){
      for(j in rownames(group_test_nonzero)){
      group_test_nonzero[j,paste("DevianceFromExcpected",i,sep="_")] <-  mean(residsNonzero[resids$G==i,j],na.rm=T)
      }
    }
   
names(group_test_nonzero)[-1] <-  paste("nonzero",names(group_test_nonzero)[-1],sep="_")
group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)

resids2 <- data.frame(cbind(resids[,setdiff(colnames(resids),colnames(taxa))],
                            resids[,setdiff(sig,sigNonzero)],
                            residsNonzero[,sigNonzero]))
colnames(resids2)<-c(setdiff(colnames(resids),colnames(taxa)),setdiff(sig,sigNonzero),sigNonzero)

    dataset2 <- data.frame(resids2)
    dataset2[,group] <- dataset2$G
    } else{    
      for(i in levels(modeldata$G)){
      for(j in rownames(group_test)){
        group_test[j,paste("Mean",i,sep="_")] <-  mean(modeldata[modeldata$G==i,j]/modeldata$ReadCount[modeldata$G==i],na.rm=T)
      }}
      for(i in levels(modeldata$G)){
      for(j in rownames(group_test_nonzero)){
        group_test_nonzero[j,paste("Mean",i,sep="_")] <-   mean(modeldata[modeldata$G==i&modeldata[,j]>0,j]/modeldata$ReadCount[modeldata$G==i&modeldata[,j]>0],na.rm=T)
      }}
      
    for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
      group_test[,paste("FoldChange",i,sep="_")] <- group_test[,paste("Mean",i,sep="_")] / group_test[,paste("Mean",compare.to,sep="_")]
    }
      
    for(i in levels(modeldata$G)[levels(modeldata$G)!= compare.to]) {
      group_test_nonzero[,paste("FoldChange",i,sep="_")] <- group_test_nonzero[,paste("Mean",i,sep="_")] / group_test_nonzero[,paste("Mean",compare.to,sep="_")]
    }
    
      names(group_test_nonzero)[-1] <- paste("nonzero",names(group_test_nonzero)[-1], sep="_") 
      group_test <- merge(group_test,group_test_nonzero,by="taxon",all=T)
      
      names(group_test)[grepl(pattern="group0",x=names(group_test))] <-gsub(pattern="group0",replacement="0",x=names(group_test)[grepl(pattern="group0",x=names(group_test))])
      
      
    dataset2 <- modeldata
    dataset2[,names(taxa)]<-log(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
    for(i in setdiff(sigNonzero,sig)) dataset2[modeldata[,i]==0,i]<-NA
      }
    write.table(group_test, paste(strsplit(taxonomic.table, split = "_")[[1]][3], 
        "_GroupTest_", group, compare.to, "_", select.by, select, ".txt", sep = ""), 
        quote = F, row.names = F, sep = "\t")
 
    
group_test$class <-  unlist(lapply(group_test$taxon, function(x) strsplit(x, split="_")[[1]][2]))
group_test$color <- as.numeric(ordered(group_test$class))        
               
for(i in unique(dataset[,"G"])[unique(dataset[,"G"])!=compare.to & !is.na(unique(dataset[,"G"]))]){ 
group_test_group <-  group_test[,c("taxon","color",colnames(group_test)[grepl(pattern=paste(i,"_",sep=""),x=colnames(group_test))])]  
a <- group_test_group[!is.na(group_test_group[,paste(i,"p",sep="_")])&!is.na(group_test_group[,paste("nonzero",i,"p",sep="_")])&group_test_group[,paste(i,"p",sep="_")]==group_test_group[,paste("nonzero",i,"p",sep="_")]|!is.na(group_test_group[,paste(i,"p",sep="_")])&!is.na(group_test_group[,paste(i,"p",sep="_")])&group_test_group[,paste("nonzero",i,"p",sep="_")]<group_test_group[,paste("nonzero",i,"p",sep="_")]|!is.na(group_test_group[,paste(i,"p",sep="_")])&is.na(group_test_group[,paste("nonzero",i,"p",sep="_")]),c("taxon","color",paste(i,"estimate",sep="_"),paste(i,"p",sep="_"))]
b <- group_test_group[!is.na(group_test_group[,paste(i,"p",sep="_")])&!is.na(group_test_group[,paste("nonzero",i,"p",sep="_")])&group_test_group[,paste(i,"p",sep="_")]>group_test_group[,paste("nonzero",i,"p",sep="_")]|is.na(group_test_group[,paste(i,"p",sep="_")])&!is.na(group_test_group[,paste("nonzero",i,"p",sep="_")]),c("taxon","color",paste("nonzero",i,"estimate",sep="_"),paste("nonzero",i,"p",sep="_"))]
names(b) <- names(a)
group_test_summary <- rbind(a,b)
                             
group_test_summary <- group_test_summary[order(group_test_summary[,3]),]
palette(c("purple","yellow", "yellowgreen", "pink","skyblue","darkorange", "turquoise2","plum", "red", "gray","royalblue","lightyellow",
                  "darkorange",   "olivedrab4", "red", "turquoise4",  "lightyellow4", "black"))

if (quartz)  quartz() 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(group_test_summary[,3],horiz=T,names.arg = paste(group_test_summary$taxon,"p =",
      round(group_test_summary[,4],digits=2)),las=2,col=group_test_summary$color,
      xaxt="n",xlab=paste("Estimated effect size of",group),main=paste(group,"=",i))   
abline(v=0)
axis(side=1) 
}
    
       
   sig <- unique(c(sig,sigNonzero))
   if(length(sig)>0){
   sig <-  na.omit(sig[apply(dataset2[,c("ReadCount",sig)],MARGIN=2,FUN=sum,na.rm=T)!=0])
    if (length(sig) > 0) {
        if (pdf) {
            pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_", 
                group, compare.to, "_", select.by, select, "_Barplot.pdf", 
                sep = ""))
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, 
                cex.axis = 1.3, cex.lab = 1.5)
        # if(length(dataset2[,i][!is.na(dataset2[,i])])>0){
            for(i in sig){
              if(length(confounders[confounders!=""])==0) {
              lab <- strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
             boxplot(dataset2[, i] ~ dataset2[, "G"], 
             #  bxp(z=list(stats=matrix(unlist(tapply(dataset2[,i],dataset2[,"G"],summary)),nrow=6,ncol=length(levels(dataset2[,"G"])))[-3,],
              #             n=table(dataset2[,"G"]),names=levels(dataset2[,"G"])),
                          yaxt="n",ylab =  lab, xlab="", las = label.direction, 
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
          #  }
         }
           dev.off()
            pdf(paste(strsplit(taxonomic.table, split = "_")[[1]][3], "_", 
                group, compare.to, "_", select.by, select, "_Beanplot.pdf", 
                sep = ""))
            par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
                1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, 
                cex.axis = 1.3, cex.lab = 1.5)
           #  if(length(dataset2[,i][!is.na(dataset2[,i])])>0){
               for(i in sig){
                if(length(confounders[confounders!=""])==0) {
                lab <- strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
              tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, group], xlab="", las = label.direction,yaxt="n",
                  ll = 0.1, ylab = lab, beanlines="median",col = list(c("skyblue", 
                    "royalblue", "royalblue", "royalblue"), c("yellowgreen", 
                    "olivedrab4", "olivedrab4", "olivedrab4"), c("pink", "red", 
                    "red", "red"), c("turquoise2", "turquoise4", "turquoise4", 
                    "turquoise4"), c("plum", "purple", "purple", "purple"), 
                    c("darkorange", "darkorange3", "darkorange3", "darkorange3"), 
                    c("lightyellow", "lightyellow4", "lightyellow4", "lightyellow4"), 
                    c("gray", "black", "black", "black")), border = c("royalblue", 
                    "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                    "lightyellow4", "black")), error = function(e) NULL)
                 if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=(yaxis))
              } else  axis(side=2)
          #  }
        }    
            dev.off()
        }
         if (quartz) quartz() else x11()
        par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
            1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, cex.axis = 1.3, 
            cex.lab = 1.5)
           for(i in sig){
            # if(length(dataset2[,i][!is.na(dataset2[,i])])>0){
          #dataset3 <- na.omit(dataset2[,c(i,group)])
           if(length(confounders[confounders!=""])==0) {
                lab <- strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
                  boxplot(dataset2[, i] ~ dataset2[, "G"], 
               #bxp(z=list(stats=matrix(unlist(tapply(dataset2[,i],dataset2[,"G"],summary)),
              #                          nrow=6,ncol=length(levels(dataset2[,"G"])))[-3,],
              #             n=table(dataset2[,"G"]),names=levels(dataset2[,"G"])),
                          yaxt="n",ylab =  lab, xlab="", las = label.direction, 
              boxfill= c("skyblue", "yellowgreen", "pink", "turquoise2", "plum", 
                  "darkorange", "lightyellow", "gray")[1:length(unique(dataset2[, 
                  "G"]))], outpch = 21, outbg = c("skyblue", "yellowgreen", 
                  "pink", "turquoise2", "plum", "darkorange", "lightyellow", 
                  "gray")[1:length(unique(dataset2[, "G"]))], outcol = c("royalblue", 
                  "olivedrab4", "red", "turquoise4", "purple", "darkorange3", 
                  "lightyellow4", "black")[1:length(unique(dataset2[, "G"]))], 
                boxcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                  "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  "G"]))], medcol = c("royalblue", "olivedrab4", "red", "turquoise4", 
                  "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  "G"]))], whiskcol = c("royalblue", "olivedrab4", "red", 
                  "turquoise4", "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  "G"]))], staplecol = c("royalblue", "olivedrab4", "red", 
                  "turquoise4", "purple", "darkorange3", "lightyellow4", "black")[1:length(unique(dataset2[, 
                  "G"]))])
              if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=exp(yaxis))
              } else  axis(side=2)
             }
         #  }
       
        if (quartz) quartz() else x11()
        par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 
            1), mgp = c(2, 0.3, 0), mar = c(7, 3.5, 1, 1), tck = -0.01, cex.axis = 1.3, 
            cex.lab = 1.5)
        
        #for (i in names(colSums(dataset2[,sig],na.rm=T))[colSums(dataset2[,sig],na.rm=T)!=0]) {
           for(i in sig){
          #      if(length(dataset2[,i][!is.na(dataset2[,i])])>0){
            if(length(confounders[confounders!=""])==0) {
                lab <- strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])]
                if(length(lab)==0) lab <- "Unassigned taxa"
                yaxis <- log(100*c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1))
                } else{
                  lab <- paste(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"][length(strsplit(i, split = "_", fixed = T)[[1]][strsplit(i, split = "_", fixed = T)[[1]]!="NA"])],"deviance")
                  if(lab==" deviance") lab <- "Unassigned taxa deviance"
                  }
            tryCatch(beanplot::beanplot(dataset2[, i] ~ dataset2[, group],
                                        xlab="", las = label.direction,yaxt="n",
                ll = 0.1, ylab = lab,beanlines="median", col = list(c("skyblue", "royalblue", 
                  "royalblue", "royalblue"), c("yellowgreen", "olivedrab4", 
                  "olivedrab4", "olivedrab4"), c("pink", "red", "red", "red"), 
                  c("turquoise2", "turquoise4", "turquoise4", "turquoise4"), 
                  c("plum", "purple", "purple", "purple"), c("darkorange", 
                    "darkorange3", "darkorange3", "darkorange3"), c("lightyellow", 
                    "lightyellow4", "lightyellow4", "lightyellow4"), c("gray", 
                    "black", "black", "black")), border = c("royalblue", "olivedrab4", 
                  "red", "turquoise4", "purple", "darkorange3", "lightyellow4", 
                  "black")), error = function(e) NULL)
                 if(length(confounders[confounders!=""])==0) { axis(side=2,at=yaxis,labels=(yaxis))
              } else  axis(side=2)
                }
           #}
    }
   }
   if(keep.result)   return(group_test)
    }
    }
    palette("default")
} 
