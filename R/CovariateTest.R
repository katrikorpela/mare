CovariateTest <- function(species.table=NULL, genus.table=NULL, family.table=NULL, order.table=NULL, class.table=NULL, phylum.table=NULL,
    meta, covariate, readcount.cutoff = 0, 
    confounders = NULL, subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05,
    group = NULL, select.by = NULL, select = NULL, pdf = F, 
    min.prevalence = 0.1, min.abundance = 0, keep.result = F, nonzero = T) {
    
    if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}
  
  metadata <- read.delim(meta)
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

  taxa <- taxa[, colSums(taxa/ metadata$ReadCount > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
   
  if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {
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
    if(nonzero){
    covariate_test[,paste(covariate,"nonzero_estimate", sep = "_")] <- NA
    covariate_test[,paste(covariate,"nonzero_p", sep = "_")] <- NA 
    }
     model <- list()
    
    if (length(group) != 0) {
        dataset$group <- as.factor(dataset[, group])
        dataset$group <- dataset$group[drop = T]
        covariate_test <- data.frame(array(dim = c(length(names(taxa)), (1 + 2*length(levels(dataset$group))))))
        names(covariate_test) <- c("taxon", c(paste(covariate, levels(dataset$group),"estimate", sep = "_"),
                                              paste(covariate, levels(dataset$group),"p", sep = "_")))
        if(nonzero){ covariate_test[,paste(covariate,levels(dataset$group),"nonzero_estimate", sep = "_")] <- NA
        covariate_test[,paste(covariate,levels(dataset$group),"nonzero_p", sep = "_")] <- NA }
        covariate_test$taxon <- names(taxa)
        rownames(covariate_test) <- names(taxa)
        
        if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])
           modeldata <- na.omit(dataset[, c("group", group,names(taxa), confounders[1], confounders[2], confounders[3], 
                      confounders[4], confounders[5], "ID", covariate, "ReadCount")[c("group",group, names(taxa),confounders[1], 
                      confounders[2], confounders[3], confounders[4], confounders[5], 
                      "ID", covariate, "ReadCount") != ""]])
           modeldata[,group]<-modeldata[,group][drop=T]
             if (max(table(dataset[, group], dataset[, subject.ID])) > 1) {
            for (j in levels(dataset$group)) { 
              
               for (i in names(taxa)) { 
                     model[[paste(i,j)]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste(i, "~", 
                    confounders[1], "+", confounders[2], "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                    data = modeldata[modeldata$group==j,], admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                      error = function(e) NULL)
              if(length(model[[paste(i,j)]])>0){
              if (model[[paste(i,j)]]$alpha<5){
          tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
            covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$coef[covariate, c(1,4)]
                }}}}}}
                  }
               for(i in rownames(covariate_test)[is.na(covariate_test[,paste(covariate, j,"estimate", sep = "_")])]){
                  model[[paste(i,j)]] <- tryCatch(nlme::lme(as.formula(paste(i, "~", 
                    confounders[1], "+", confounders[2], "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1|ID,  data = modeldata[modeldata$group==j,]), 
                      error = function(e) NULL)
              if(length(model[[paste(i,j)]])>0){
         tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                   covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$tTable[covariate, c(1,5)]
          }}}}}
                  }}
             
        if(nonzero){         
        for (j in levels(dataset$group)) { 
          for(i in names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3])){     
                 
            model[[paste(i,j,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste(i, "~",
                       confounders[1], "+", confounders[2], "+", 
                      confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1 | 
                      ID, family = "nbinom", data = modeldata[modeldata$group==j&modeldata[,i]>0,],
                      admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                      error = function(e) NULL)
              if(length(model[[paste(i,j,"nonzero")]])>0){
              if (model[[paste(i,j,"nonzero")]]$alpha<5){
               tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){       
            covariate_test[i, c(paste(covariate, j,"nonzero_estimate", sep = "_"),paste(covariate, j,"nonzero_p", sep = "_"))] <- summary(model[[paste(i,j,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}
                  }
        
          for(i in intersect(rownames(covariate_test)[is.na(covariate_test[,paste(covariate, j,"nonzero_estimate", sep = "_")])], 
            names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3]))){     
                    model[[paste(i,j,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste(i, "~",
                       confounders[1], "+", confounders[2], "+", 
                      confounders[3], "+", confounders[4], "+", confounders[5], 
                      "+", covariate, "+", "offset(log(ReadCount))")), random = ~1 | ID,  
                      data = modeldata[modeldata$group==j&modeldata[,i]>0,]), 
                      error = function(e) NULL)
              if(length(model[[paste(i,j,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
         covariate_test[i, c(paste(covariate, j,"nonzero_estimate", sep = "_"),paste(covariate, j,"nonzero_p", sep = "_"))] <- summary(model[[paste(i,j,"nonzero")]])$tTable[covariate, c(1,5)]
               }}}}}
        }}
        }
         }
           } else {
    modeldata <- na.omit(dataset[ , c("group",group, names(taxa), confounders[1], confounders[2], confounders[3], confounders[4], 
                    confounders[5],  covariate, "ReadCount")[c("group",group, names(taxa), confounders[1], 
                    confounders[2], confounders[3], confounders[4], confounders[5], 
                    covariate, "ReadCount") != ""]])
    modeldata[,group]<-modeldata[,group][drop=T]
      for (j in levels(dataset$group)) {      
        
        for (i in names(taxa)){
                  model[[paste(i,j)]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate, 
                    "+", "offset(log(ReadCount))")), data = modeldata[modeldata$group==j,], 
                    control = glm.control(maxit = 500)), 
                    error = function(e) NULL)
          if(length(model[[paste(i,j)]])>0){
          if ((model[[paste(i,j)]]$deviance/model[[paste(i,j)]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[paste(i,j)]],type="deviance"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                  covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$coef[covariate, c(1,4)]
          }}}}}}}
         
        for(i in rownames(covariate_test)[is.na(covariate_test[,paste(covariate,j,"p",sep="_")])]){
                  model[[paste(i,j)]] <- tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", 
                                                                      confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate)), data = modeldata[modeldata$group==j,]),
                    error = function(e) NULL)
          if(length(model[[paste(i,j)]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,j)]],type="deviance"), pred = predict(model[[paste(i,j)]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                  covariate_test[i, c(paste(covariate, j,"estimate", sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$coef[covariate, c(1,4)]
          }}}}}}
        
        }
          
        
          if(nonzero){       
            for (j in levels(dataset$group)) {  
              
          for(i in names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3])){     
          model[[paste(i,j,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate, 
                    "+", "offset(log(ReadCount))")),  data = modeldata[modeldata$group==j&modeldata[,i]>0,],
                    control = glm.control(maxit = 500)), 
                    error = function(e) NULL)
                   if(length(model[[paste(i,j,"nonzero")]])>0){
          if ((model[[paste(i,j,"nonzero")]]$deviance/model[[paste(i,j,"nonzero")]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[paste(i,j,"nonzero")]],type="deviance"), pred = predict(model[[paste(i,j,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j&modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                covariate_test[i, c(paste(covariate, j,"nonzero_estimate", sep = "_"),paste(covariate, j,"nonzero_p", sep = "_"))] <- summary(model[[paste(i,j,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}}
      
    for(i in intersect(rownames(covariate_test)[is.na(covariate_test[,paste(covariate,j,"nonzero_p",sep="_")])], 
                       names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3]))){     
          model[[paste(i,j,"nonzero")]] <- tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                    "+", confounders[4], "+", confounders[5], "+", covariate)),  
                    data = modeldata[modeldata$group==j&modeldata[,i]>0,]), 
                    error = function(e) NULL)
                   if(length(model[[paste(i,j,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,j,"nonzero")]],type="deviance"), pred = predict(model[[paste(i,j,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata$group==j&modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                covariate_test[i, c(paste(covariate, j,"nonzero_estimate", sep = "_"),paste(covariate, j,"nonzero_p", sep = "_"))] <- summary(model[[paste(i,j,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}   
              }
          }
 }
          
  for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, i])
  
      sig <- na.omit(names(apply(covariate_test[,grepl(pattern="_p$",x=colnames(covariate_test))],MARGIN = 1,FUN = min,na.rm=T)[apply(covariate_test[,grepl(pattern="_p$",x=colnames(covariate_test))],MARGIN = 1,FUN = min,na.rm=T)<p.cutoff]))
      if(nonzero) {
        sigNonzero <- na.omit(names(apply(covariate_test[,grepl(pattern="_nonzero_p$",x=colnames(covariate_test))],MARGIN = 1,FUN = min,na.rm=T)[apply(covariate_test[,grepl(pattern="_nonzero_p$",x=colnames(covariate_test))],MARGIN = 1,FUN = min,na.rm=T)<p.cutoff]))
      } else  sigNonzero <-  NULL
      
       for (j in levels(dataset$group)) {
            covariate_test[, paste(covariate, j, "FDR", sep = "_")] <- p.adjust(covariate_test[, paste(covariate, j, "p", sep = "_")], "fdr")
         if(nonzero)   covariate_test[, paste(covariate, j, "nonzero_FDR", sep = "_")] <- p.adjust(covariate_test[, paste(covariate, j, "nonzero_p", sep = "_")], "fdr")
            }
        covariate_test$taxon <- rownames(covariate_test)
        write.table(covariate_test, paste("CovariateTest_", covariate, "_", group, "_", select.by, select, ".txt", sep = ""), 
            quote = F, row.names = F, sep = "\t")
      
        if (length(sig) > 0) {
          
          if(length(confounders[confounders!=""])>0){
      resids <- modeldata
       resids[,rownames(covariate_test)] <- NA
            
       for (j in levels(modeldata$group)) {

    for(i in setdiff(sig,sigNonzero)){
    resids[resids$group==j,i] <-  tryCatch(resid(update(model[[paste(i,j)]],as.formula(paste(".~. -",covariate,sep="")))),
                              error = function(e) NA)
    }
 if(nonzero) {
      for(i in sigNonzero){
     resids[resids$group==j&modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,j,"nonzero")]],
                                                    as.formula(paste(".~. -",covariate,sep="")))),
                               error = function(e) NA)
     }
 }
           }
   dataset2 <- data.frame(resids)
            } else {
              dataset2 <- modeldata
               dataset2[,names(taxa)]<-log(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
              if(nonzero) for(i in sigNonzero){ dataset2[modeldata[,i]==0,i] <- NA  } 
             
            }
 
                
  df = na.omit(reshape2::melt(dataset2[,c(covariate,group,sig)], id=c(covariate,group)))
  names(df) <- c("x","gr","variable","value") 
  df$variable <- sapply(df$variable, function(x) gsub(pattern = "_NA",replacement = "", x=x))
  df$variable <- sapply(df$variable, function(x) gsub(pattern = "_incertae_sedis",replacement = "", x=x))
  df$variable <- unlist(sapply(as.character(df$variable), function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])]))
  
  if(length(confounders[confounders!=""])==0) { 
  p<- ggplot2::ggplot(df, ggplot2::aes(y=value, x=x,color=factor(gr)),environment = environment()) +
  ggplot2::stat_smooth(method = "lm", formula = y ~ x, ggplot2::aes(fill=factor(gr)),se=T) +
  ggplot2::facet_wrap(~variable,ncol=floor(sqrt(length(sig))),scales='free')+
  ggplot2::theme_bw()+
  ggplot2::xlab(covariate)+
  ggplot2::geom_point(pch=20) + 
  ggplot2::scale_color_manual(name=group,values=c('gray30','red','turquoise4','olivedrab4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c('cornflowerblue','pink','turquoise2','yellowgreen','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])+
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
  ggplot2::ylab('Relative abundance') +
  ggplot2::scale_y_continuous(breaks=log(c(0.0001,0.001,0.01,0.1,1,10,100)),labels=c(0,0.001,0.01,0.1,1,10,100))
    } else {
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
    }
  if(pdf){
pdf(paste(covariate,"_",group, "_", select.by,select, "_", "Covariateplot.pdf", sep = ""));
plot(p)
dev.off() 
}       
 quartz() 
  plot(p)

#------------
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

sp <- intersect(sp,colnames(taxa))  
 
for(h in unique(modeldata[,group])){ 
covariate_test_group <-  cbind(covariate_test[,colnames(covariate_test)[grepl(pattern=paste("_",h,"_",sep=""),x=colnames(covariate_test))]])  
if(nonzero){
a <- covariate_test_group[!is.na(covariate_test[,paste(covariate,h,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,h,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,h,"p",sep="_")]==covariate_test[,paste(covariate,h,"nonzero_p",sep="_")]|!is.na(covariate_test[,paste(covariate,h,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,h,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,h,"p",sep="_")]<covariate_test[,paste(covariate,h,"nonzero_p",sep="_")]|!is.na(covariate_test[,paste(covariate,h,"p",sep="_")])&is.na(covariate_test[,paste(covariate,h,"nonzero_p",sep="_")]),c(paste(covariate,h,"estimate",sep="_"),paste(covariate,h,"p",sep="_"))]
b <- covariate_test_group[!is.na(covariate_test[,paste(covariate,h,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,h,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,h,"p",sep="_")]>covariate_test[,paste(covariate,h,"nonzero_p",sep="_")]|is.na(covariate_test[,paste(covariate,h,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,h,"nonzero_p",sep="_")]),c(paste(covariate,h,"nonzero_estimate",sep="_"),paste(covariate,h,"nonzero_p",sep="_"))]
names(b) <- names(a)
covariate_test_summary <- rbind(a,b)
} else covariate_test_summary <- covariate_test_group 
if(nrow(na.omit(covariate_test_summary))>0){
covariate_test_summary$name <- sapply(rownames(covariate_test_summary), function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])
covariate_test_summary[,paste(covariate,h,"estimate",sep="_")][covariate_test_summary[,paste(covariate,h,"p",sep="_")] > p.cutoff] <- 0
covariate_test_summary2 <- covariate_test_summary[intersect(sp,rownames(covariate_test_summary)),]
  
seqs1 <- list()
for(i in rownames(covariate_test_summary2)){
  if(!is.na(covariate_test_summary2[i,paste(covariate,h,"p",sep="_")]))  if(covariate_test_summary2[i,paste(covariate,h,"p",sep="_")]<p.cutoff)  seqs1[[i]]<- covariate_test_summary2[i,paste(covariate,h,"estimate",sep="_")]
}

for(i in rownames(covariate_test_summary2)){
  if(!is.na(covariate_test_summary2[i,paste(covariate,h,"p",sep="_")]))  if(covariate_test_summary2[i,paste(covariate,h,"p",sep="_")]>p.cutoff) seqs1[[i]]<- 0
}

abu <- list()
for(i in names(seqs1)) abu[[i]] <-  mean(reltaxa[,i])

newseqs2 <- paste("XX", names(seqs1),seqs1, abu, sep=" ")

tmp2 <-  metacoder::extract_taxonomy(input=newseqs2,regex = "^(.*)\\ (.*)\\ (.*)\\ (.*)",
                         key=c(id = "obs_info","class","taxon_info","taxon_info"),class_sep = "_")

for(k in rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)&tmp2$taxon_data$name!="NA"]){
if(tmp2$taxon_data$name[as.numeric(k)]%in%covariate_test_summary$name) tmp2$taxon_data$taxon_info_1[as.numeric(k)] <- covariate_test_summary[covariate_test_summary$name==tmp2$taxon_data$name[as.numeric(k)],paste(covariate,h,"estimate",sep="_")][1]
}
for(k in rev(rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)])){
 tmp2$taxon_data$taxon_info_1[as.numeric(k)] <-  mean(as.numeric(tmp2$taxon_data$taxon_info_1[!is.na(tmp2$taxon_data$supertaxon_ids)&tmp2$taxon_data$supertaxon_ids==tmp2$taxon_data$taxon_ids[as.numeric(k)]]),na.rm=T)
}


treltaxa <- as.data.frame(t(reltaxa))
treltaxa$name <-   sapply(rownames(treltaxa), function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])

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
                   node_color_axis_label = "Effect size",
                   node_size_axis_label = "Average relative abundance (%)",
          node_label_size_range=c(0.01,0.014),
          node_label_size_trans="ln",
          #node_label_size=0.03,
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          overlap_avoidance = 0.5, node_label_max=150, make_legend=T, 
          node_size_trans="linear",
          node_size_range=c(0.012,0.05),
          node_color_trans="linear",
          title=paste("Associations with", covariate,"in group", h), title_size=0.03,
          output_file=paste(covariate,"_",h, "_", select.by, select,"_HeatTree.pdf",sep=""))


}
}
#------------

        }
      
    } else {
              
        if (length(subject.ID) != 0) {
            dataset$ID <- as.factor(dataset[, subject.ID])
            modeldata <- na.omit(dataset[, c(names(taxa), confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5],  "ID", 
                    covariate, "ReadCount")[c(names(taxa), confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5],  "ID", 
                    covariate, "ReadCount") != ""]])
          
            for (i in names(taxa)){
               model[[i]] <-  tryCatch(glmmADMB::glmmadmb(as.formula(paste(i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                  data = modeldata, admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                  error = function(e) NULL)
          if(length(model[[i]])>0){
          if ((model[[i]]$alpha)<5){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]], type="response"))  
          tmp$covariate <- modeldata[,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
      covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$coef[covariate, c(1,4)]
          }}}}}}}
            
      for(i in rownames(covariate_test)[is.na(covariate_test[,paste(covariate, "estimate", sep = "_")])]){
           model[[i]] <-  tryCatch(nlme::lme(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate)), random = ~1 | ID,
                  data = modeldata), 
                  error = function(e) NULL)
          if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]], type="response"))  
          tmp$covariate <- modeldata[,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
           covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$tTable[covariate, c(1,5)]
          }}}}}}
            
        if(nonzero){     
       for (i in names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3])){      
               model[[paste(i,"nonzero")]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste( i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), random = ~1 | ID, family = "nbinom", 
                  data = modeldata[modeldata[,i]>0,], 
                  admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                  error = function(e) NULL)
              if(length(model[[paste(i,"nonzero")]])>0){
          if ((model[[paste(i,"nonzero")]]$alpha)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
             covariate_test[i, c(paste(covariate, "nonzero_estimate", sep = "_"),paste(covariate, "nonzero_p", sep = "_"))] <-  summary(model[[paste(i,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}}
              
          for(i in intersect(rownames(covariate_test)[is.na(covariate_test[,paste(covariate, "nonzero_estimate", sep = "_")])], 
                             names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3]))){     
               model[[paste(i,"nonzero")]] <- tryCatch(nlme::lme(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate)), random = ~1 | ID,  
                  data = modeldata[modeldata[,i]>0,]), 
                  error = function(e) NULL)
              if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
             covariate_test[i, c(paste(covariate, "nonzero_estimate", sep = "_"),paste(covariate, "nonzero_p", sep = "_"))] <-  summary(model[[paste(i,"nonzero")]])$tTable[covariate, c(1,5)]
          }}}}}}
        }     
        
       
        } else {
          modeldata <- na.omit(dataset[, c(names(taxa), confounders[1], confounders[2], confounders[3], confounders[4], 
                    confounders[5], covariate, "ReadCount")[c(names(taxa), confounders[1], 
                    confounders[2], confounders[3], confounders[4], confounders[5], 
                     covariate, "ReadCount") != ""]])
         
           for (i in names(taxa)){
                model[[i]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), data = modeldata, control = glm.control(maxit = 1000)), 
                  error = function(e) NULL)
               if(length(model[[i]])>0){
               if ((model[[i]]$deviance/model[[i]]$df.residual)<5){
              tmp <- data.frame(res = resid(model[[i]],type="deviance"), pred = predict(model[[i]], type="response"))  
          tmp$covariate <- modeldata[,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
                 covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]] )$coef[covariate, c(1,4)]
          }}}}}}}
          
          for(i in rownames(covariate_test)[is.na(covariate_test[,paste(covariate, "estimate", sep = "_")])]){
           model[[i]] <-  tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate)),data = modeldata), 
                  error = function(e) NULL)
          if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]], type="response"))  
          tmp$covariate <- modeldata[,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
           covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$coef[covariate, c(1,4)]
          }}}}}}
          
        if(nonzero){ 
          for (i in names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3])){    
            model[[paste(i,"nonzero")]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate, 
                  "+", "offset(log(ReadCount))")), data = modeldata[modeldata[,i]>0,], 
                  control = glm.control(maxit = 1000)), 
                  error = function(e) NULL)
           if(length(model[[paste(i,"nonzero")]])>0){
              if ((model[[paste(i,"nonzero")]]$deviance/model[[paste(i,"nonzero")]]$df.residual)<5){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="deviance"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
           covariate_test[i, c(paste(covariate, "nonzero_estimate", sep = "_"),paste(covariate, "nonzero_p", sep = "_"))] <- summary(model[[paste(i,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}}
          for(i in intersect(rownames(covariate_test)[is.na(covariate_test[,paste(covariate, "nonzero_estimate", sep = "_")])], 
                            names(colSums(modeldata[,colnames(taxa)]>0)[colSums(modeldata[,colnames(taxa)]>0)>3]))){     
               model[[paste(i,"nonzero")]] <- tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~", confounders[1], "+", confounders[2], "+", confounders[3], 
                  "+", confounders[4], "+", confounders[5], "+", covariate)),  
                  data = modeldata[modeldata[,i]>0,]), 
                  error = function(e) NULL)
              if(length(model[[paste(i,"nonzero")]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,"nonzero")]],type="pearson"), pred = predict(model[[paste(i,"nonzero")]], type="response"))  
          tmp$covariate <- modeldata[modeldata[,i]>0,covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
          if(anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.05){
             covariate_test[i, c(paste(covariate, "nonzero_estimate", sep = "_"),paste(covariate, "nonzero_p", sep = "_"))] <-  summary(model[[paste(i,"nonzero")]])$coef[covariate, c(1,4)]
          }}}}}}
        }    
        }
        for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[,i])
  
    sig <-  na.omit(rownames(covariate_test)[covariate_test[,grepl(pattern=paste(covariate,"p$",sep="_"),x=colnames(covariate_test))]<p.cutoff])
    if(nonzero){   sigNonzero <-  na.omit(rownames(covariate_test)[covariate_test[,grepl(pattern="nonzero_p$",x=colnames(covariate_test))]<p.cutoff]) 
        } else  sigNonzero <- NULL
         
         sig <- unique(c(sig, sigNonzero)) 
 
          covariate_test[, paste(covariate, "FDR", sep = "_")] <- p.adjust(covariate_test[,paste(covariate, "p", sep = "_")], "fdr")
      if(nonzero) covariate_test[, paste(covariate, "nonzero_FDR", sep = "_")] <- p.adjust(covariate_test[, paste(covariate, "nonzero_p", sep = "_")], "fdr")
        
        covariate_test$taxon <- rownames(covariate_test)
        write.table(covariate_test, paste("CovariateTest_", covariate, "_", select.by, select, ".txt", sep = ""), quote = F, 
            row.names = F, sep = "\t")       
              
if (length(sig) > 0) {
          
    if(length(confounders[confounders!=""])>0){
    resids<- modeldata
    resids[,rownames(covariate_test)] <- NA
   
   if(nonzero){
    for(i in sigNonzero){
    resids[modeldata[,i]>0,i] <-  tryCatch(resid(update(model[[paste(i,"nonzero")]],
                                                    as.formula(paste(".~. -",covariate,sep="")))),
                               error = function(e) NA)
    }
   }
    for(i in setdiff(sig,c(sigNonzero))){
     resids[,i] <-   tryCatch(resid(update(model[[i]],as.formula(paste(".~. -",covariate,sep="")))),
                               error = function(e) NA)
    }
    
    dataset2 <- data.frame(resids)
          } else {
              dataset2 <- modeldata
               dataset2[,names(taxa)]<-log(100*((dataset2[,names(taxa)]+1)/dataset2$ReadCount))
              if(nonzero) for(i in sigNonzero){ dataset2[modeldata[,i]==0,i] <- NA  } 
           
          }
 
 df = na.omit(reshape2::melt(dataset2[,c(covariate,sig)], id=c(covariate)))
  names(df) <- c("x","variable","value") 
  df$variable <- sapply(df$variable, function(x) gsub(pattern = "_NA",replacement = "", x=x))
  df$variable <- sapply(df$variable, function(x) gsub(pattern = "_incertae_sedis",replacement = "", x=x))
  df$variable <- unlist(sapply(as.character(df$variable), function(x) strsplit(x, split = "_", fixed = T)[[1]][length(strsplit(x, split = "_", fixed = T)[[1]])]))
 
  if(length(confounders[confounders!=""])==0) { 
  p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
  ggplot2::stat_smooth(data=df,method = "lm", formula = y ~ x, se = T, fill = "cornflowerblue", color = "gray30", lwd = 0.5)+
  ggplot2::facet_wrap(~variable,ncol=floor(sqrt(length(sig))),scales='free')+
  ggplot2::theme_bw()+
  ggplot2::xlab(covariate)+
  ggplot2::geom_point(pch=20,color = "gray30") +  
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
      ggplot2::ylab('Relative abundance') +
      ggplot2::scale_y_continuous(breaks=log(c(0.0001,0.001,0.01,0.1,1,10,100)),labels=c(0,0.001,0.01,0.1,1,10,100))
   } else{
      p <- ggplot2::ggplot(df, ggplot2::aes(y = value, x = x), environment = environment()) + 
     ggplot2::stat_smooth(data=df,method = "lm", formula = y ~ x, se = T, fill = "cornflowerblue", color = "gray30", lwd = 0.5)+
   ggplot2::facet_wrap(~variable,ncol=floor(sqrt(length(sig))),scales='free')+
  ggplot2::theme_bw()+
  ggplot2::xlab(covariate)+
          ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
  ggplot2::geom_point(pch=20,color = "gray30") +  
 ggplot2::ylab('Deviance from expected') 
   }

  if(pdf){
pdf(paste(covariate,"_",group, "_", select.by,select, "_", "Covariateplot.pdf", sep = ""));
plot(p)
dev.off() 
}       
quartz() 
  plot(p)
         
      
if (pdf) {
  pdf(paste(covariate,"_",group, "_", select.by,select, "_", "CovariateTestResult.pdf", sep = "")) 
par(mar=c(4,20,2,2),cex.axis=0.5)
barplot(covariate_test_summary[,2],horiz=T,names.arg = paste(rownames(covariate_test_summary),"p =",
      round(covariate_test_summary[,3],digits=2)),las=2,col=as.factor(covariate_test_summary$class),
      xaxt="n",xlab=paste("Estimated effect size of",covariate))   
abline(v=0)
axis(side=1) 
dev.off()}


#------------
#------------
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

sp <- intersect(sp,colnames(taxa))  
  
if(nonzero){
a <- covariate_test[!is.na(covariate_test[,paste(covariate,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,"p",sep="_")]==covariate_test[,paste(covariate,"nonzero_p",sep="_")]|!is.na(covariate_test[,paste(covariate,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,"p",sep="_")]<covariate_test[,paste(covariate,"nonzero_p",sep="_")]|!is.na(covariate_test[,paste(covariate,"p",sep="_")])&is.na(covariate_test[,paste(covariate,"nonzero_p",sep="_")]),c(paste(covariate,"estimate",sep="_"),paste(covariate,"p",sep="_"))]
b <- covariate_test[!is.na(covariate_test[,paste(covariate,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,"nonzero_p",sep="_")])&covariate_test[,paste(covariate,"p",sep="_")]>covariate_test[,paste(covariate,"nonzero_p",sep="_")]|is.na(covariate_test[,paste(covariate,"p",sep="_")])&!is.na(covariate_test[,paste(covariate,"nonzero_p",sep="_")]),c(paste(covariate,"nonzero_estimate",sep="_"),paste(covariate,"nonzero_p",sep="_"))]
names(b) <- names(a)
covariate_test_summary <- rbind(a,b)
} else covariate_test_summary <- covariate_test 
covariate_test_summary$name <- sapply(rownames(covariate_test_summary), function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])
covariate_test_summary[,paste(covariate,"estimate",sep="_")][covariate_test_summary[,paste(covariate,"p",sep="_")] > p.cutoff] <- 0
covariate_test_summary2 <- covariate_test_summary[intersect(sp,rownames(covariate_test_summary)),]
  
seqs1 <- list()
for(i in rownames(covariate_test_summary2)){
  if(!is.na(covariate_test_summary2[i,paste(covariate,"p",sep="_")]))  if(covariate_test_summary2[i,paste(covariate,"p",sep="_")]<p.cutoff)  seqs1[[i]]<- covariate_test_summary2[i,paste(covariate,"estimate",sep="_")]
}

for(i in rownames(covariate_test_summary2)){
  if(!is.na(covariate_test_summary2[i,paste(covariate,"p",sep="_")]))  if(covariate_test_summary2[i,paste(covariate,"p",sep="_")]>p.cutoff) seqs1[[i]]<- 0
}

abu <- list()
for(i in names(seqs1)) abu[[i]] <-  mean(reltaxa[,i])

newseqs2 <- paste("XX", names(seqs1),seqs1, abu, sep=" ")

tmp2 <-  metacoder::extract_taxonomy(input=newseqs2,regex = "^(.*)\\ (.*)\\ (.*)\\ (.*)",
                         key=c(id = "obs_info","class","taxon_info","taxon_info"),class_sep = "_")

for(k in rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)&tmp2$taxon_data$name!="NA"]){
if(tmp2$taxon_data$name[as.numeric(k)]%in%covariate_test_summary$name) tmp2$taxon_data$taxon_info_1[as.numeric(k)] <- covariate_test_summary[covariate_test_summary$name==tmp2$taxon_data$name[as.numeric(k)],paste(covariate,"estimate",sep="_")][1]
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
                   node_color_axis_label = "Effect size",
                   node_size_axis_label = "Average relative abundance (%)",
          node_label_size_range=c(0.01,0.014),
          node_label_size_trans="ln",
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          overlap_avoidance = 0.5, node_label_max=150, make_legend=T, 
          node_size_trans="linear",
          node_size_range=c(0.012,0.05),
          node_color_trans="linear",
          title=paste("Associations with", covariate), title_size=0.03,
          output_file=paste(covariate, "_", select.by, select,"_HeatTree.pdf",sep=""))




#------------
}



    }         
      if(keep.result)        return(covariate_test)
    } 
    palette("default")
}