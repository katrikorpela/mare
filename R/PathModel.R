PathModel <- function(taxonomic.table, meta, response = NULL, binomial.response = F, variables, readcount.cutoff = 0, 
    subject.ID = NULL,  outlier.cutoff, 
    select.by = NULL, select = NULL, pdf = F, 
    min.prevalence = 0, min.abundance = 0, keep.result = F) {
    
    if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}


    metadata <- read.delim(meta)
    taxa <- read.delim(taxonomic.table)
    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]

    if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {
    taxa <- taxa[metadata$ReadCount > readcount.cutoff, ]
    metadata <- metadata[metadata$ReadCount > readcount.cutoff, ]
    background.variables <- names(metadata[,variables])
    
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
    
    dataset <- data.frame(metadata, reltaxa)
    taxa <- round(reltaxa * metadata$ReadCount - 1)
    taxa[taxa<0]<-0
    dataset2 <- data.frame(metadata, taxa)
   
    modeldata <- na.omit(dataset[, c(names(taxa), response, background.variables)])
    modeldata2 <- na.omit(dataset2[, c(names(taxa), response, background.variables)])
    
    modeldata[,c(names(taxa),background.variables)] <- scale(modeldata[,c(names(taxa),background.variables)])
    modeldata2[,background.variables] <- scale(modeldata2[,background.variables])
   
   if(length(response)>0){
     
   covariate_test <- data.frame(var= c(names(taxa),background.variables), 
                                p = rep(NA, length(c(names(taxa),background.variables))))
   names(covariate_test)[2] <- "p"
   rownames(covariate_test) <- c(names(taxa),background.variables)

     model <- list()
     
    if(binomial.response){
      for (i in c(names(taxa),background.variables)){
      if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])  
      model[[i]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste(response, "~",i)), random = ~1|ID, 
                                                admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE),
                                                data = modeldata, family="binomial"), 
                  error = function(e) NULL)
      } else {model[[i]] <- tryCatch(glm(as.formula(paste(response, "~",i)), data = modeldata, family="binomial"), 
                  error = function(e) NULL)
      }
               if(length(model[[i]])>0){
              if(length(na.omit(coef(model[[i]])))>1){  
              if ((model[[i]]$deviance/model[[i]]$df.residual)<5){
     tmp <- data.frame(res = resid(model[[i]],type="deviance"), pred = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res-0) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05 ){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
            covariate_test[i,"p"] <- summary(model[[i]])$coef[2,4]
          }} }  }}
      }
    } else{
       for (i in c(names(taxa),background.variables)){
        
       if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])  
      model[[i]] <- tryCatch(nlme:lme(as.formula(paste(response, "~",i)), random = ~1|ID, data = modeldata), 
                  error = function(e) NULL)
       } else{
        model[[i]] <- tryCatch(lm(as.formula(paste(response, "~",i)), data = modeldata), 
                  error = function(e) NULL)
      }    
         
               if(length(model[[i]])>0){
              if(length(na.omit(coef(model[[i]])))>1){  
     tmp <- data.frame(res = resid(model[[i]],type="deviance"), pred = predict(model[[i]]))  
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05 ){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
            covariate_test[i,"p"] <- summary(model[[i]])$coef[2,4]
          }} }  }
    }
    }  
  
  
    sig <- rownames(covariate_test)[order(covariate_test$p)][1:(nrow(modeldata)/3)]
    
    if(length(sig)>0){
    if (length(sig) > 1) sig2 <- paste(sig[1],sig[2],sep=" + ") else sig2 <- sig
    if (length(sig) > 2) for(i in 3:length(sig)) sig2 <-  paste(sig2,sig[i],sep=" + ")
    if(binomial.response){
      combmodel <- step(glm(as.formula(paste(response, "~",sig2)), data = modeldata, family="binomial"),k=3)
    } else{
       combmodel <- step(lm(as.formula(paste(response, "~",sig2)), data = modeldata),k=3)
    }
    sigbac <- intersect(names(taxa),names(combmodel$coefficients)[-1])
    sigvars <- setdiff(names(combmodel$coefficients)[-1],names(taxa))
      bacmodel <- list()
      bac_test <- list()
      combbacmodel <- list()
      
      if(length(sigbac)>0){
    for(i in sigbac){    
   bac_test[[i]] <- data.frame(var= background.variables, p = rep(NA, length(background.variables)))
   names(bac_test[[i]])[2] <- "p"
   rownames(bac_test[[i]]) <- background.variables

        for(j in background.variables){
           if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])  
            bacmodel[[paste(i,j)]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste(i, "~",j)),random = ~1 | ID, family = "nbinom",
                                                                  data = modeldata2, admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                  error = function(e) NULL)
           }else {
           bacmodel[[paste(i,j)]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~",j)), data = modeldata2, 
                                               control = glm.control(maxit = 1000)), 
                  error = function(e) NULL)  
           }
               if(length(bacmodel[[paste(i,j)]])>0){
                  if(length(na.omit(coef(bacmodel[[paste(i,j)]])))>1){ 
               if ((bacmodel[[paste(i,j)]]$deviance/bacmodel[[paste(i,j)]]$df.residual)<5){
              tmp <- data.frame(res = resid(bacmodel[[paste(i,j)]],type="deviance"), pred = predict(bacmodel[[paste(i,j)]]))  
          tmp$resdev <- abs(tmp$res-0) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
                 bac_test[[i]][j,"p"] <- summary(bacmodel[[paste(i,j)]])$coef[2,4]
          }}}}}
        }

    sigi <- rownames(bac_test[[i]])[order(bac_test[[i]]$p)][1:(nrow(modeldata2)/3)]
    if(length(sigi)>0){
    if (length(sigi) > 1) sigi2 <- paste(sigi[1],sigi[2],sep=" + ") else sigi2 <- sigi
    if (length(sigi) > 2) for(h in 3:length(sigi)) sigi2 <-  paste(sigi2,sigi[h],sep=" + ")
    combbacmodel[[i]] <- tryCatch(MASS::stepAIC(MASS::glm.nb(as.formula(paste(i, "~",sigi2)), data = modeldata2),k=3),
                                   error = function(e) NULL)
   
      }
    }
    }
  
    } else print(paste("No variables remained in the final model"))
    } else {
  sigbac <- colnames(taxa)
  bacmodel <- list()
  bac_test <- list()
  combbacmodel <- list()
  
  for(i in sigbac){    
   bac_test[[i]] <- data.frame(var= background.variables, p = rep(NA, length(background.variables)))
   names(bac_test[[i]])[2] <- "p"
   rownames(bac_test[[i]]) <- background.variables

        for(j in background.variables){
    if (length(subject.ID) != 0) {
           dataset$ID <- as.factor(dataset[, subject.ID])  
            bacmodel[[paste(i,j)]] <- tryCatch(glmmADMB::glmmadmb(as.formula(paste(i, "~",j)),random = ~1 | ID, family = "nbinom",
                                                                  data = modeldata2, admb.opts = glmmADMB::admbControl(shess = F, noinit = FALSE)), 
                  error = function(e) NULL)
           }else {
           bacmodel[[paste(i,j)]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~",j)), data = modeldata2, 
                                               control = glm.control(maxit = 1000)), 
                  error = function(e) NULL)  
           }
               if(length(bacmodel[[paste(i,j)]])>0){
                  if(length(na.omit(coef(bacmodel[[paste(i,j)]])))>1){ 
               if ((bacmodel[[paste(i,j)]]$deviance/bacmodel[[paste(i,j)]]$df.residual)<5){
              tmp <- data.frame(res = resid(bacmodel[[paste(i,j)]],type="deviance"), pred = predict(bacmodel[[paste(i,j)]]))  
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05){
          if(summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05){
                 bac_test[[i]][j,"p"] <- summary(bacmodel[[paste(i,j)]])$coef[2,4]
          }}}}}
        }

    sigi <- na.omit(rownames(bac_test[[i]])[order(bac_test[[i]]$p)][1:round(nrow(modeldata2)/3)])
    if(length(sigi)>0){
    if (length(sigi) > 1) sigi2 <- paste(sigi[1],sigi[2],sep=" + ") else sigi2 <- sigi
    if (length(sigi) > 2) for(h in 3:length(sigi)) sigi2 <-  paste(sigi2,sigi[h],sep=" + ")
    combbacmodel[[i]] <- tryCatch(MASS::stepAIC(MASS::glm.nb(as.formula(paste(i, "~",sigi2)), data = modeldata2),k=3),
                                   error = function(e) NULL)
   
      }
    }    
      
      
    }  
  
    result <- matrix(nrow=length(c(response,sigbac,background.variables)),ncol=length(c(response,sigbac,background.variables)))
    rownames(result) <- c(response,sigbac,background.variables)
    colnames(result) <- c(response,sigbac,background.variables)

    if(length(response)>0){
   result[names(combmodel$coefficients)[-1],response] <- combmodel$coef[-1]
    }
 if(length(sigbac)>0){
   for(i in sigbac){
   result[names(combbacmodel[[i]]$coefficients[-1]),i] <-  combbacmodel[[i]]$coefficients[-1] 
   } 
    }
  
res <- result[unique(c(rownames(result)[rowSums(result,na.rm=T)!=0],colnames(result)[colSums(result,na.rm=T)!=0])),
               unique(c(rownames(result)[rowSums(result,na.rm=T)!=0],colnames(result)[colSums(result,na.rm=T)!=0]))] 
  
namesres <- vector("list", ncol(res))
names(namesres)<-colnames(res)
namesres[intersect(colnames(res),names(taxa))] <-  lapply(intersect(colnames(res),names(taxa)), function(x) gsub(pattern="_NA",replacement = ".",x))
namesres[intersect(colnames(res),names(taxa))] <- unlist(lapply(unlist(namesres[intersect(colnames(res),names(taxa))]) , function(x) strsplit(x, split = "_")[[1]][length(strsplit(x, split = "_")[[1]])]))
namesres[setdiff(colnames(res),names(taxa))] <- setdiff(colnames(res),names(taxa))

class <-  vector("list", ncol(res))
names(class)<-colnames(res)
class[intersect(colnames(res),names(taxa))] <-  unlist(lapply(intersect(colnames(res),names(taxa)), function(x) strsplit(x, split="_")[[1]][1]))
class[setdiff(colnames(res),names(taxa))] <-  "1"
if (length(response)>0) class[[response]] <- "2"
color <- as.numeric(ordered(unlist(class))) 


curves <- res
curves[!is.na(curves)]<-0.25

if (length(response)>0) {
palette(c("gray","white","yellow", "yellowgreen", "pink","skyblue","darkorange", "turquoise2","plum", "red", "aliceblue","royalblue","lightyellow",
                  "darkorange",   "olivedrab4", "firebrick4", "turquoise4",  "lightyellow4"))
} else {
palette(c("gray","yellow", "yellowgreen", "pink","skyblue","darkorange", "turquoise2","plum", "red", "aliceblue","royalblue","lightyellow",
                  "darkorange",   "olivedrab4", "firebrick4", "turquoise4",  "lightyellow4"))
  
}

if(pdf) {
  pdf(paste("PathModel_",response,".pdf",sep=""))
qgraph::qgraph(input=res,asize=5,vsize=10,esize=10,directed=T,colFactor=0,arrowAngle=0.6,
           labels=substr(namesres,start=1,stop=15),groups=class,color=color,legend=F,
           curve.all=T,borders=F,negCol = "red",posCol="yellowgreen",
           edge.label.cex=1,edge.labels=T)
dev.off()
}
quartz()
qgraph::qgraph(input=res,asize=5,vsize=8,esize=10,directed=T,colFactor=0,arrowAngle=0.6,
           labels=substr(namesres,start=1,stop=15),groups=class,color=color,legend=F,
          curveShape=-0.5,curveScale=T,
          curve=curves, edge.label.cex=1,edge.labels=T, borders=F,negCol = "red",posCol="yellowgreen")
  

    if (keep.result) return(list(combmodel,combbacmodel))
}
}    
    
    
    

     
     
     
     
     
     
     
     
