PathModel <- function(taxonomic.table, meta, response = NULL, 
                      binomial.response = F, count.response = F, variables = NULL, readcount.cutoff = 0, 
                      outlier.cutoff = 3, 
                      select.by = NULL, select = NULL, pdf = F, 
                      min.prevalence = 0, min.abundance = 0, keep.result = F, relative = T, k.value = 3.5,
                      logtrans = F) {
  
    if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}

    metadata <- read.delim(meta)
    if(!relative) metadata[,"ReadCount"]<-1
    taxa <- read.delim(taxonomic.table)
    if(logtrans) taxa <- log(taxa + min(taxa[taxa>0],na.rm=T))
    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    taxa <- taxa[, colSums(taxa/rowSums(taxa) > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]

    if(ncol(taxa)==0) print("No taxa that fullfill the abundance and prevalence criteria!")
    if(ncol(taxa)>0) {
      
    taxa <- taxa[metadata$ReadCount > readcount.cutoff, ]
    metadata <- metadata[metadata$ReadCount > readcount.cutoff, ]
   if(length(variables)>1) background.variables <- names(metadata[,variables]) else background.variables <- variables
    
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
    modeldata2 <- na.omit(dataset2[, c(names(taxa), response, background.variables, "ReadCount")])
    
    modeldata[,c(names(taxa), background.variables)] <- scale(modeldata[,c(names(taxa), background.variables)])
    modeldata2[,c(background.variables)] <- scale(modeldata2[,c(background.variables)])
    
    
   if(length(response)>0){
     
   covariate_test <- data.frame(var= c(names(taxa),background.variables), 
                                p = rep(NA, length(c(names(taxa),background.variables))))
   names(covariate_test)[2] <- "p"
   rownames(covariate_test) <- c(names(taxa),background.variables)

     model <- list()
     
    if(binomial.response){
      for (i in c(names(taxa),background.variables)){
        model[[i]] <- tryCatch(glm(as.formula(paste(response, "~",i)), data = modeldata, family="binomial"), 
                  error = function(e) NULL)
      
               if(length(model[[i]])>0){
              if(length(na.omit(coef(model[[i]])))>1){  
            covariate_test[i,"p"] <- summary(model[[i]])$coef[2,4]
                 }}
      }
    } else if (count.response){
      
    for (i in c(names(taxa),background.variables)){
        model[[i]] <- tryCatch(MASS::glm.nb(as.formula(paste(response, "~",i)), data = modeldata), 
                  error = function(e) NULL)
               if(length(model[[i]])>0){
              if(length(na.omit(coef(model[[i]])))>1){  
            covariate_test[i,"p"] <- summary(model[[i]])$coef[2,4]
                }}
      }
      
    } else {
       for (i in c(names(taxa),background.variables)){
        
        model[[i]] <- tryCatch(lm(as.formula(paste(response, "~",i)), data = modeldata), 
                  error = function(e) NULL)
               if(length(model[[i]])>0){
              if(length(na.omit(coef(model[[i]])))>1){  
         covariate_test[i,"p"] <- summary(model[[i]])$coef[2,4] 
                 } }  
    }
    }  
  
  
    sig <- na.omit(rownames(covariate_test)[order(covariate_test$p)][1:(nrow(modeldata)/10)])
    
    if(length(sig)>0){
    if (length(sig) > 1) sig2 <- paste(sig[1],sig[2],sep=" + ") else sig2 <- sig
    if (length(sig) > 2) for(i in 3:length(sig)) sig2 <-  paste(sig2,sig[i],sep=" + ")
    if(binomial.response){
 combmodel <- step(glm(as.formula(paste(response, "~",sig2)), data = modeldata, family="binomial"),k=k.value)
    } else if(count.response){
      combmodel <- step(MASS::glm.nb(as.formula(paste(response, "~",sig2)), data = modeldata),k=k.value)
    } else {
      combmodel <- step(lm(as.formula(paste(response, "~",sig2)), data = modeldata),k=k.value)
    }
    
    
    if(binomial.response == F){
    if(length(combmodel)>0){
      if(length(na.omit(coef(combmodel)))>1){  
     tmp <- data.frame(res = resid(combmodel,type="pearson"), pred = predict(combmodel))  
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]<0.01){
          combmodel <-  tryCatch(nlme::gls(formula(combmodel), weights = nlme::varExp(), 
                    data = modeldata, control = nlme::glsControl(maxIter=5000)), 
                      error = function(e) NULL) 
         }
              }} else print(paste("No variables remained in the final model"))
    }
    
    if(length(combmodel)>0){

     sigbac <- intersect(names(taxa),names(combmodel$coefficients)[-1])
    sigvars <- setdiff(names(combmodel$coefficients)[-1],names(taxa)) 
    
      bacmodel <- list()
      bac_test <- list()
      combbacmodel <- list()
      
      if(length(sigbac)>0){
    for(i in sigbac){    
   bac_test[[i]] <- data.frame(var = background.variables, p = rep(NA, length(background.variables)))
   names(bac_test[[i]])[2] <- "p"
   rownames(bac_test[[i]]) <- background.variables

        for(j in background.variables){
         
           bacmodel[[paste(i,j)]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~",j, "+ offset(log(ReadCount))")), data = modeldata2, 
                                               control = glm.control(maxit = 1000)), 
                  error = function(e) NULL)  
           
            if(length(bacmodel[[paste(i,j)]])==0){
           bacmodel[[paste(i,j)]] <- tryCatch(lm(as.formula(paste("log((", i, "+1)/ReadCount)~",j)), data = modeldata2), 
                  error = function(e) NULL)
              }
          
           bac_test[[i]][j,"p"] <- tryCatch(summary(bacmodel[[paste(i,j)]])$coef[2,4],
                                               error = function(e) NA)
        }


    sigi <- na.omit(rownames(bac_test[[i]])[order(bac_test[[i]]$p)][1:(nrow(modeldata2)/10)])
    if(length(sigi)>0){
    if (length(sigi) > 1) sigi2 <- paste(sigi[1],sigi[2],sep=" + ") else sigi2 <- sigi
    if (length(sigi) > 2) for(h in 3:length(sigi)) sigi2 <-  paste(sigi2,sigi[h],sep=" + ")
    

      combbacmodel[[i]] <- tryCatch(step(MASS::glm.nb(as.formula(paste(i, "~",sigi2, "+ offset(log(ReadCount))")), 
                                                      data = modeldata2)),
                                   error = function(e) NULL)
      
       if (length(combbacmodel[[i]])==0){
             combbacmodel[[i]]<- tryCatch(step(lm(as.formula(paste("log((",i, "+1)/ReadCount)~",sigi2)),
                                                data = modeldata2)), 
                  error = function(e) NULL) 
            
      
       }
      }
    }
    }
  
    } else print("Model fit could not be found. Try adding/removing background variables.")
}
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

           bacmodel[[paste(i,j)]] <- tryCatch(MASS::glm.nb(as.formula(paste(i, "~",j, "+", "offset(log(ReadCount))")), 
                                                           data = modeldata2, 
                                               control = glm.control(maxit = 1000)), 
                  error = function(e) NULL) 
           
             if (length(bacmodel[[paste(i,j)]])==0){
             bacmodel[[paste(i,j)]] <- tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~",j)),
                                                data = modeldata2), 
                  error = function(e) NULL) 
            
           
           }
               if(length(bacmodel[[paste(i,j)]])>0){
                  if(length(na.omit(coef(bacmodel[[paste(i,j)]])))>1){ 
                  if(length(summary(bacmodel[[paste(i,j)]])$tTable)==0)  bac_test[[i]][j,"p"] <- summary(bacmodel[[paste(i,j)]])$coef[2,4]
                 if(length(summary(bacmodel[[paste(i,j)]])$tTable)!=0)  bac_test[[i]][j,"p"] <- summary(bacmodel[[paste(i,j)]])$tTable[2,5]

                 }}
        }

    sigi <- na.omit(rownames(bac_test[[i]])[order(bac_test[[i]]$p)][1:round(nrow(modeldata2)/10)])
    if(length(sigi)>0){
    if (length(sigi) > 1) sigi2 <- paste(sigi[1],sigi[2],sep=" + ") else sigi2 <- sigi
    if (length(sigi) > 2) for(h in 3:length(sigi)) sigi2 <-  paste(sigi2,sigi[h],sep=" + ")
    
    combbacmodel[[i]] <- tryCatch(MASS::stepAIC(MASS::glm.nb(as.formula(paste(i, "~",sigi2, "+ offset(log(ReadCount))")),
                                                             data = modeldata2),k=k.value),
                                   error = function(e) NULL)
    if (length(combbacmodel[[i]])==0){
            combbacmodel[[i]] <- tryCatch(lm(as.formula(paste("log((",i, "+1)/ReadCount)~",sigi2)),
                                                data = modeldata2), 
                  error = function(e) NULL) 
    }
       if (length(combbacmodel[[i]])>0){
    tmp <- data.frame(res = resid(combbacmodel[[i]],type="deviance"), pred = predict(combbacmodel[[i]]))  
    tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]<0.01){
                combbacmodel[[i]] <- NULL
          }
       
    }
      }
    }    
    }  
  
  
    
   
   result <- matrix(nrow=length(c(response,sigbac,background.variables)),
                     ncol=length(c(response,sigbac,background.variables)))
    rownames(result) <- c(response,sigbac,background.variables)
    colnames(result) <- c(response,sigbac,background.variables)

    if(length(combmodel)>0 ){
      if(length(response)>0){
result[names(combmodel$coefficients)[-1],response] <- combmodel$coef[-1]        
   }   

    }
    
    
  
    
if(length(combbacmodel)>0 ){ 
 if(length(sigbac)>0){
   for(i in sigbac){
   result[names(combbacmodel[[i]]$coefficients[-1]),i] <-  combbacmodel[[i]]$coefficients[-1] 

    }
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
class[intersect(colnames(res),names(taxa))] <-  unlist(lapply(intersect(colnames(res),names(taxa)), function(x) strsplit(x, split="_")[[1]][2]))
class[setdiff(colnames(res),names(taxa))] <-  "1"
if (length(response)>0) class[[response]] <- "2"
color <- as.numeric(ordered(unlist(class))) 

curves <- res
curves[!is.na(curves)]<-0.25

alpha = floor(255*0.5)  
if (length(response)>0) {
newColor = col2rgb(col=c("gray90","black","#E41A1C","#FFA500","#377EB8","#87CEFA","#4DAF4A" ,'#9ACD32',"#984EA3",'#DA70D6', "#999999","gainsboro",
      "#008080","#00CED1","#F781BF","thistle1","#8DA0CB","lightsteelblue1","#FFD92F","#FFFFB3",
"#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
"#CCEBC5","#FFED6F","#C71585","#EE82EE","#66C2A5","#FC8D62","#A65628"), alpha=FALSE)
} else{
newColor = col2rgb(col=c("gray90","#E41A1C","#FFA500","#377EB8","#87CEFA","#4DAF4A" ,'#9ACD32',"#984EA3",'#DA70D6', "#999999","gainsboro",
      "#008080","#00CED1","#F781BF","thistle1","#8DA0CB","lightsteelblue1","#FFD92F","#FFFFB3",
"#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
"#CCEBC5","#FFED6F","#C71585","#EE82EE","#66C2A5","#FC8D62","#A65628"), alpha=FALSE)
}
 makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newCol = apply(newColor, 2, makeTransparent, alpha=alpha)
  palette(newCol)

  
if(pdf) {
  pdf(paste("PathModel_",select.by,select,"_",response,".pdf",sep=""))
qgraph::qgraph(input=res,asize=5,vsize=8,esize=10,directed=T,colFactor=0,arrowAngle=0.6,
           labels=substr(namesres,start=1,stop=15),groups=class,color=color,legend=F,
          curveShape=-0.5,curveScale=T,
          curve=curves, 
          edge.label.cex=1,edge.labels=F, 
          borders=F,negCol = "red",posCol="yellowgreen")
dev.off()
}
quartz()
qgraph::qgraph(input=res,asize=5,vsize=8,esize=10,directed=T,colFactor=0,arrowAngle=0.6,
           labels=substr(namesres,start=1,stop=15),groups=class,color=color,legend=F,
          curveShape=-0.5,curveScale=T,
          curve=curves, 
          edge.label.cex=1,edge.labels=F, 
          borders=F,negCol = "red",posCol="yellowgreen")
  

    if (keep.result) return(list(combmodel,combbacmodel))

}
}    
    
    
    

     
     
     
     
     
     
     
     
