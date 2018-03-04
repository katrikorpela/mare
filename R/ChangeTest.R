ChangeTest <- function(species.table=NULL, genus.table=NULL, family.table=NULL, order.table=NULL, class.table=NULL, phylum.table=NULL, 
                      meta, group = NULL, compare.to = NULL, 
    covariate = NULL, readcount.cutoff = 0, confounders = NULL, subject.ID,  time,
    outlier.cutoff = 3, p.cutoff = 0.05, select.by = NULL, select = NULL, pdf = F, 
    consecutive = T, min.prevalence = 0, min.abundance = 0, label.direction = 1, keep.result =F) {
    
      
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
  
    colnames(genus)[grepl(pattern="incertae_sedis",x=colnames(genus))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(genus)[grepl(pattern="incertae_sedis",x=colnames(genus))])
    colnames(genus)[grepl(pattern="Incertae_Sedis_",x=colnames(genus))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(genus)[grepl(pattern="Incertae_Sedis_",x=colnames(genus))])
    colnames(genus)[grepl(pattern="Erysipelotrichi_",x=colnames(genus))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(genus)[grepl(pattern="Erysipelotrichi_",x=colnames(genus))])

    colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))] <- gsub(pattern="_incertae_sedis",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="incertae_sedis",x=colnames(taxa))])
    colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))] <-  gsub(pattern="_Incertae_Sedis_",replacement = "incertaesedis",x= colnames(taxa)[grepl(pattern="Incertae_Sedis_",x=colnames(taxa))])
    colnames(taxa)[grepl(pattern="Erysipelotrichi_",x=colnames(taxa))] <-  gsub(pattern="Erysipelotrichi_",replacement = "Erysipelotrichia_",x= colnames(taxa)[grepl(pattern="Erysipelotrichi_",x=colnames(taxa))])
  
 taxa <- taxa[, colSums(taxa/ metadata$ReadCount > min.abundance, na.rm = T) > min.prevalence * nrow(taxa)]
       
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
            metadata$time[metadata$ID == i] <- order(metadata$time[metadata$ID == i])
            for (j in unique(metadata$time[metadata$ID == i][order(metadata$time[metadata$ID == i])])[-1]) {
                for (k in colnames(reltaxa)) {
                  deltataxa[metadata$ID == i & metadata$time == j, k] <- reltaxa[metadata$ID == i & metadata$time == j, k] - 
                    reltaxa[metadata$ID == i &  metadata$time == (j - 1), k]
                  tryCatch(deltataxa[metadata$ID == i & metadata$time == j, 
                    paste("baseline", k, sep = "")] <- reltaxa[metadata$ID == i & metadata$time == (j - 1), k], 
                    error = function(e) NULL)
                }
            }
        }
    } else {
        for (k in colnames(deltataxa)) deltataxa[, paste("baseline", k, sep = "")] <- NA
        for (i in names(table(metadata$ID)[table(metadata$ID) > 1])) {
          if(1%in%metadata$time[metadata$ID==i]){
            for (j in unique(metadata$time[metadata$ID == i])[unique(metadata$time[metadata$ID == i])!=1]) {
                for (k in colnames(reltaxa)) {
                  deltataxa[metadata$ID == i & metadata$time == j, k] <- reltaxa[metadata$ID ==  i & metadata$time == j, k] - 
                    reltaxa[metadata$ID == i & metadata$time == 1, k]
                  tryCatch(deltataxa[metadata$ID == i & metadata$time == j, 
                    paste("baseline", k, sep = "")] <- reltaxa[metadata$ID ==  i & metadata$time == 1, k], 
                    error = function(e) NULL)
                }
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
      
    dataset[,group] <- as.character(dataset[,group])
    dataset <- dataset[!is.na(dataset[,group]),]
    dataset[,group][dataset[,group]==""] <- "nogroup"
    for(i in unique(dataset[,group])) if(table(dataset[,group])[i] < 2) dataset[,group][dataset[,group]==i] <- "toofewcases"
    dataset <- dataset[dataset[,group]!="toofewcases",]
    dataset[, group] <- dataset[, group][drop = T]
    if(length(compare.to)==0) compare.to = levels(as.factor(dataset[,group]))[1]
   dataset$G <- as.factor(dataset[,group])
    if (compare.to != "0") dataset[, group][dataset[, group] == "0" & !is.na(dataset[, group])] <- "group0"
    dataset[, group][dataset[, group] == compare.to & !is.na(dataset[, group])] <- "0"
    dataset[, group] <- as.factor(dataset[, group])
   
    for(i in c(1:ncol(dataset))[-c(1:ncol(metadata),ncol(dataset))]){
      if(min(dataset[,i], na.rm = T) == max(dataset[,i], na.rm = T)) dataset[,i] <- NA
     }
    
    grouptime <- paste(dataset[dataset$time != 1 & dataset[, group] != "0" & dataset[,group]!="other","time"], 
                       dataset[dataset$time != 1 & dataset[, group] != "0" & dataset[,group]!="other", group], sep = "/")
        
  modeldata <- na.omit(dataset[dataset[,group]!="other", c("time","G", subject.ID,time,group, confounders[1], confounders[2], confounders[3], 
                  confounders[4], confounders[5], names(taxa), paste("baseline", names(taxa), sep = ""))[c("time", "G",subject.ID,group, time,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = "")) != ""]])
  modeldata[, group] <- modeldata[, group][drop=T]
        

        group_test <- data.frame(array(dim = c(length(names(taxa)), 1 + 2*(length(levels(as.factor(grouptime)))))))
        rownames(group_test) <- names(taxa)
        names(group_test) <- c("taxon", paste("estimate_",levels(as.factor(grouptime)),sep=""),paste("p_",levels(as.factor(grouptime)),sep=""))
        
        for (i in names(taxa)) {
            for (k in names(table(modeldata$time)[table(modeldata$time)>3])) {
              model[[paste(i,k)]] <-  tryCatch(lm(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], "+ ", group, sep = "")), 
                  data = modeldata[modeldata$time == k&!is.na(modeldata[,i]),]), 
                  error = function(e) NULL)
               if(length(model[[paste(i,k)]])>0 & nrow(summary(model[[paste(i,k)]])$coef)>1 & summary(model[[paste(i,k)]])$coef[2,4]!="NaN"){
          tmp <- data.frame(res = resid(model[[paste(i,k)]],type="deviance"), fitted = predict(model[[paste(i,k)]]))  
          tmp$group <- modeldata[modeldata$time == k&!is.na(modeldata[,i]),group]
          tmp$resdev <- abs(tmp$res)
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.01 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.01 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.01 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.01){
          for(j in levels(modeldata[,group])[-1]){
            tryCatch(group_test[i, c(paste("estimate_",k,"/",j,sep=""),paste("p_",k,"/",j,sep=""))] <- summary(model[[paste(i,k)]])$coef[paste(group,j,sep=""),c(1,4)],
               error = function(e) NULL  )
            }   
          } else {
            if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]<0.01){
             model[[paste(i,k)]] <-   tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], "+ ", group, sep = "")), 
                  data = modeldata[modeldata$time == k&!is.na(modeldata[,i]),], weights = nlme::varExp(), 
                  control = nlme::glsControl(maxIter=5000)), 
                  error = function(e) NULL)
          } else if(anova(lm(res ~ group, data=tmp))$Pr[1]<0.01 | anova(lm(resdev ~ group, data=tmp))$Pr[1]<0.01){
             model[[paste(i,k)]] <-  tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], "+ ", group, sep = "")), 
                  data = modeldata[modeldata$time == k&!is.na(modeldata[,i]),], weights = nlme::varIdent(~group),
                  control = nlme::glsControl(maxIter=5000)), 
                  error = function(e) NULL)
          }
          if(length( model[[paste(i,k)]])>0){
          tmp <- data.frame(res = resid(model[[paste(i,k)]],type="pearson"), fitted = predict(model[[paste(i,k)]]))  
          tmp$resdev <- abs(tmp$res) 
         tmp$group <- modeldata[modeldata$time == k&!is.na(modeldata[,i]),group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.01 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.01 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.01 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.01){
            for(j in levels(modeldata[,group])[-1]){
            tryCatch(group_test[i, c(paste("estimate_",k,"/",j,sep=""),paste("p_",k,"/",j,sep=""))] <- summary(model[[paste(i,k)]])$tTable[paste(group,j,sep=""),c(1,4)],
            error = function(e) NULL)
              } 
          } else {
            model[[paste(i,k)]] <-  tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], "+ ", group, sep = "")), 
                  data = modeldata[modeldata$time == k&!is.na(modeldata[,i]),], 
                  weights = nlme::varComb(nlme::varExp(),nlme::varIdent(~group)),
                  control = nlme::glsControl(maxIter=5000)), 
                  error = function(e) NULL)
            if(length(model[[paste(i,k)]])>0){
           tmp <- data.frame(res = resid(model[[paste(i,k)]],type="pearson"), fitted = predict(model[[paste(i,k)]]))  
          tmp$resdev <- abs(tmp$res) 
         tmp$group <- modeldata[modeldata$time == k&!is.na(modeldata[,i]),group]
          if(summary(lm(tmp$res ~ tmp$fitted))$coef[2,4]>0.01 & summary(lm(tmp$resdev ~ tmp$fitted))$coef[2,4]>0.01 &
             anova(lm(res ~ group, data=tmp))$Pr[1]>0.01 & anova(lm(resdev ~ group, data=tmp))$Pr[1]>0.01){
           for(j in levels(modeldata[,group])[-1]){
            tryCatch(group_test[i, c(paste("estimate_",k,"/",j,sep=""),paste("p_",k,"/",j,sep=""))] <- summary(model[[paste(i,k)]])$tTable[paste(group,j,sep=""),c(1,4)],
            error = function(e) NULL)
             }  
          }
            }
          }
          }
          }
          }
            }
        }
        
    colnames(group_test) <- gsub(x=colnames(group_test),pattern="group0",replacement = "0")
            
     group_test$taxon <- rownames(group_test)
   
      sig <- as.character(rownames(group_test)[sapply(data.frame(t(group_test[, grepl(pattern="p_",x=names(group_test))])), min,na.rm=T) < p.cutoff])
     for (k in names(group_test)[grepl(pattern="p_",x=names(group_test))]) group_test[, paste(k, "FDR", sep = "_")] <- p.adjust(group_test[,k], "fdr")
   
      resids <- modeldata
  
    for(i in rownames(group_test)){
      for(k in unique(modeldata$time)){
     tryCatch( resids[resids$time==k&!is.na(resids[,i]),i] <-  resid(update(model[[paste(i,k)]], as.formula(paste(".~. -",group,sep="")))),
                               error = function(e) NULL)
    }}
 

 write.table(group_test, paste("ChangeTest_", group, compare.to, "_", select.by, select, ".txt", 
            sep = ""), quote = F, row.names = F, sep = "\t")
   
    dataset2 <- data.frame(resids[,c(sig,"G",time)])
    dataset2[,group] <- dataset2$G
    dataset2$grouptime <-paste(dataset2[, time], dataset2[, group], sep="/")
      
    
        if (length(sig) > 0) {
        if (pdf) {
            pdf(paste(group, compare.to, "_", select.by, select, "_Boxplot.pdf", sep = ""))
   par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 1), 
      mgp = c(2, 0.5, 0), mar = c(3, 3, 2, 0.5), tck = -0.01, 
                cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste("Delta ",strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])])
                  if(lab=="Delta ") lab <- "Delta Unassigned taxon"
                  
           boxplot(dataset2[, i] ~ dataset2[, "grouptime"], yaxt="n", ylab =  lab, xlab="", las = label.direction, 
              col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,
       "#A65628", "#F781BF", "#999999","blue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
       outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,
       "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
               axis(side=2)
         }
         
           dev.off()
            
        }
  quartz()
  par(mfrow = c(floor(sqrt(length(sig))), round(sqrt(length(sig))) + 1), 
      mgp = c(2, 0.5, 0), mar = c(3, 3, 2, 0.5), tck = -0.01, 
                cex.axis = 1.2, cex.lab = 1,cex.main=0.8)
            for(i in sig[order(sig)]){
                 lab <- gsub(i,pattern="_NA",replacement = ".")  
                 lab <- paste("Delta ",strsplit(lab, split = "_", fixed = T)[[1]][length(strsplit(lab, split = "_", fixed = T)[[1]])])
                  if(lab=="Delta ") lab <- "Delta Unassigned taxon"
                  
           boxplot(dataset2[, i] ~ dataset2[, "grouptime"], yaxt="n", ylab =  lab, xlab="", las = label.direction, 
              col = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,
       "#A65628", "#F781BF", "#999999","blue","firebrick4",
       'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[,"G"]))], 
       outpch = 21, outbg = c("#E41A1C","orange","#377EB8","skyblue","#4DAF4A" ,"#984EA3", "#FF7F00" ,
       "#A65628", "#F781BF", "#999999","blue","firebrick4",'yellowgreen','pink','turquoise2','plum','darkorange','lightyellow','gray',
       'royalblue','olivedrab4','red','turquoise4','purple','darkorange3','lightyellow4','black')[1:length(unique(dataset2[, "G"]))])
               axis(side=2)
         }

 
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
  for(t in unique(dataset[,"time"])[unique(dataset[,"time"])>1 & !is.na(unique(dataset[,"time"]))]){
group_test_summary <- tryCatch(group_test[,c("taxon",paste("estimate_",t,"/",h,sep=""),paste("p_",t,"/",h,sep=""))],
         error = function(e) NULL)
if(length(group_test_summary)>0){
if(length(na.omit(group_test_summary[,2]))>0){
group_test_summary$name <- sapply(group_test_summary$taxon, function(x) strsplit(x, split="_")[[1]][length(strsplit(x, split="_")[[1]])])
rownames(group_test_summary)<-  group_test_summary$taxon
group_test_summary[,paste("estimate_",t,"/",h,sep="")][group_test_summary[,paste("p_",t,"/",h,sep="")] > p.cutoff] <- 0
group_test_summary2 <- group_test_summary[intersect(sp,rownames(group_test_summary)),]

seqs1 <- list()
for(i in rownames(group_test_summary2)){
  if(!is.na(group_test_summary2[i,paste("p_",t,"/",h,sep="")]))  if(group_test_summary2[i,paste("p_",t,"/",h,sep="")]<p.cutoff)  seqs1[[group_test_summary2[i,"taxon"]]]<- group_test_summary2[i,paste("estimate_",t,"/",h,sep="")]
}

for(i in rownames(group_test_summary2)){
  if(!is.na(group_test_summary2[i,paste("p_",t,"/",h,sep="")]))  if(group_test_summary2[i,paste("p_",t,"/",h,sep="")]>p.cutoff) seqs1[[group_test_summary2[i,"taxon"]]]<- 0
}

abu <- list()
for(i in names(seqs1)) abu[[i]] <-  mean(reltaxa[,i])

#newseqs2 <- paste("XX", names(seqs1),seqs1, abu, sep=" ")
 if(length(intersect(colnames(phylum),sp))==length(sp)){newseqs2 <- paste("XX ", names(seqs1)," ",seqs1," ", abu, sep="")
#} else if( length(sp[grepl("GH_",sp)])==0) {newseqs2 <- paste("XX Bacteria_", names(seqs1)," ",seqs1," ", abu, sep="")
 } else newseqs2 <- paste("XX Bacteria_", names(seqs1)," ",seqs1," ", abu, sep="")

tmp2 <- extract_taxonomy(input=newseqs2,regex = "^(.*)\\ (.*)\\ (.*)\\ (.*)",
                         key=c(id = "obs_info","class","taxon_info","taxon_info"),class_sep = "_")

for(k in rownames(tmp2$taxon_data)[is.na(tmp2$taxon_data$taxon_info_1)&tmp2$taxon_data$name!="NA"]){
if(tmp2$taxon_data$name[as.numeric(k)]%in%group_test_summary$name) tmp2$taxon_data$taxon_info_1[as.numeric(k)] <- group_test_summary[group_test_summary$name==tmp2$taxon_data$name[as.numeric(k)],paste("estimate_",t,"/",h,sep="")][1]
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
                   node_color_range=c("royalblue","gray95","red"),
                   node_color_interval = c(-3,3),#c(-max(abs(as.numeric(tmp2$taxon_data$taxon_info_1))), max(abs(as.numeric(tmp2$taxon_data$taxon_info_1)))),
                   edge_color_interval = c(-3,3),#c(-max(abs(as.numeric(tmp2$taxon_data$taxon_info_1))), max(abs(as.numeric(tmp2$taxon_data$taxon_info_1)))),
                   node_color_axis_label = paste(group,h,"compared to",compare.to),
                   node_size_axis_label = "Average relative abundance (%)",
          node_label_size_range=c(0.01,0.013),
          node_label_size_trans="ln",
          #node_label_size=0.03,
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          overlap_avoidance = 0.5, node_label_max=150, make_legend=T, 
          node_size_trans="linear",
          node_size_range=c(0.012,0.05),
          node_color_trans="linear",
          title=paste("Differences between groups", h, "and", compare.to), title_size=0.03,
          output_file=paste("Time",t,group,h,"vs", compare.to, "_", select.by, select,"_HeatTree.pdf",sep=""))
}
}
  }  
}
#------------
            
            
            
            
        }
       if(keep.result)   return(group_test)
    
#------------
 
    }
   
    if (length(covariate) != 0) {
        if (length(group) != 0) {

    dataset[,group]<-as.character(dataset[,group])
    for(i in unique(dataset[,group])) if(table(dataset[,group])[i] < 3) dataset[,group][dataset[,group]==i] <- "toofewcases"
    dataset <- dataset[dataset[,group]!="toofewcases",]
    dataset[, group] <- dataset[, group][drop = T]
   dataset$G <- as.factor(dataset[,group])
    dataset[, group] <- as.factor(dataset[, group])

   covariate_test <- data.frame(array(dim = c(length(names(taxa)), (1 + 2*length(levels(dataset[,group]))))))
   names(covariate_test) <- c("taxon", c(paste(covariate, levels(dataset[,group]),"estimate", sep = "_"),
                                              paste(covariate, levels(dataset[,group]),"p", sep = "_")))
  covariate_test$taxon <- names(taxa)
  rownames(covariate_test) <- names(taxa)
        
  modeldata <- na.omit(dataset[, c(group,subject.ID, covariate, confounders[1], confounders[2], confounders[3], 
                    confounders[4], confounders[5], names(taxa), paste("baseline", names(taxa), 
                      sep = ""))[c(group, subject.ID, covariate, confounders[1], confounders[2], 
                    confounders[3], confounders[4], confounders[5], names(taxa), paste("baseline", 
                      names(taxa), sep = "")) != ""]])
    

    for (i in names(taxa)) {
    for (j in unique(modeldata[,group])) {
   model[[paste(i,j)]] <-  tryCatch(lm(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),]),
                    error = function(e) NULL)
    if(length(model[[paste(i,j)]])>0 & nrow(summary(model[[paste(i,j)]])$coef)>1 & summary(model[[paste(i,j)]])$coef[covariate,4]!="NaN" & 
       summary(model[[paste(i,j)]])$coef[covariate,4]!="NA"){ 
   
    tmp <- data.frame(res = resid(model[[paste(i,j)]],type="deviance"), pred = predict(model[[paste(i,j)]]))  
          tmp$covariate <- modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res) 
         
    if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.01 & 
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.01){

   covariate_test[i, c(paste(covariate, j, "estimate",sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$coef[covariate,c(1,4)]
    } else {
      if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]<0.01 | summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]<0.01){
        model[[paste(i,j)]] <-  tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),],
                     weights = nlme::varExp(),  control = nlme::glsControl(maxIter=5000)),
                    error = function(e) NULL)   
    } else {
        model[[paste(i,j)]] <-  tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),],
                     weights = nlme::varExp(form= as.formula(paste("~",covariate))),  
                    control = nlme::glsControl(maxIter=5000)),
                    error = function(e) NULL)   
    }
          
     if(length(model[[paste(i,j)]])>0){ 
        tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]]))  
          tmp$covariate <- modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res)
      if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.01 & 
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.01){
         covariate_test[i, c(paste(covariate, j, "estimate",sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$tTable[covariate,c(1,4)]
      } else {
         model[[paste(i,j)]] <-  tryCatch(nlme::gls(as.formula(paste(i, "~baseline", i, "+", confounders[1], "+", confounders[2], 
                    "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                    "+", covariate, sep = "")), data = modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),],
                     weights = nlme::varComb(nlme::varExp(),nlme::varExp(form= as.formula(paste("~",covariate)))),  
                    control = nlme::glsControl(maxIter=5000)),
                    error = function(e) NULL)  
         
    if(length(model[[paste(i,j)]])>0){ 
        tmp <- data.frame(res = resid(model[[paste(i,j)]],type="pearson"), pred = predict(model[[paste(i,j)]]))  
          tmp$covariate <- modeldata[modeldata[,group]==j&!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res)
      if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(res ~ covariate, data=tmp))$Pr[1]>0.01 & 
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & anova(lm(resdev ~ covariate, data=tmp))$Pr[1]>0.01){
         covariate_test[i, c(paste(covariate, j, "estimate",sep = "_"),paste(covariate, j,"p", sep = "_"))] <- summary(model[[paste(i,j)]])$tTable[covariate,c(1,4)]
      }
      }}
     } 
      } }}}
            
          for (i in names(covariate_test)[-1]) covariate_test[, i] <- as.numeric(covariate_test[, i])
         sig <- na.omit(names(apply(covariate_test[,grepl(pattern="_p",x=colnames(covariate_test))],MARGIN = 1,FUN = min)[apply(covariate_test[,grepl(pattern="_p",x=colnames(covariate_test))],MARGIN = 1,FUN = min)<p.cutoff]))
   
            for (j in levels(dataset[,group])) {
                covariate_test[, paste(covariate, j, "p","FDR", sep = "_")] <- p.adjust(covariate_test[, 
                  paste(covariate, j,"p",  sep = "_")], "fdr")
            }

            write.table(covariate_test, paste("ChangeTest_", covariate, "_", group, "_", select.by, select, 
                ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
            
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
        
    df = na.omit(reshape2::melt(dataset2[, c(covariate, group, sig)], id = c(covariate, group)))
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
  ggplot2::scale_color_manual(name=group,values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'gray30','red','turquoise4','olivedrab4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'cornflowerblue','pink','turquoise2','yellowgreen','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])+
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
  ggplot2::ylab('Deviance from expected') 

                if (pdf) {
                  pdf(paste("ChangeTest", covariate, "_", group, "_", select.by, select, 
                    "_Covariateplot.pdf", sep = ""))
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
covariate_test_summary <- covariate_test_group 
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

#newseqs2 <- paste("XX", names(seqs1),seqs1, abu, sep=" ")
 if(length(intersect(colnames(phylum),sp))==length(sp)){newseqs2 <- paste("XX ", names(seqs1)," ",seqs1," ", abu, sep="")
#} else if( length(sp[grepl("GH_",sp)])==0) {newseqs2 <- paste("XX Bacteria_", names(seqs1)," ",seqs1," ", abu, sep="")
 } else newseqs2 <- paste("XX Bacteria_", names(seqs1)," ",seqs1," ", abu, sep="")

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
          node_label_size_range=c(0.01,0.013),
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
    covariate_test <- data.frame(taxon = names(taxa), estimate = rep(NA, length(names(taxa))), p = rep(NA, length(names(taxa))))
    names(covariate_test)[2:3] <- c(paste(covariate,"estimate",sep="_"),paste(covariate,"p",sep="_"))
    rownames(covariate_test) <- names(taxa)
       
            modeldata <- na.omit(dataset[, c(covariate, subject.ID,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = ""))[c(covariate, subject.ID,
                  confounders[1], confounders[2], confounders[3], confounders[4], 
                  confounders[5], names(taxa), paste("baseline", names(taxa), sep = "")) != ""]])
 
            for (i in names(taxa)) {
                model[[i]] <- tryCatch(lm(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = modeldata), 
                  error = function(e) NULL)
          if(length(model[[i]])>0){
          tmp <- data.frame(res = resid(model[[i]],type="deviance"), pred = predict(model[[i]]))  
          tmp$covariate <- modeldata[!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(res ~ covariate, data=tmp))$coef[2,4]>0.01 &
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(resdev ~ covariate, data=tmp))$coef[2,4]>0.01){
                  covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$coef[covariate,c(1, 4)]
         } else {
           if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.05 | summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.05 ){
            model[[i]] <- tryCatch(nlme::gls(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = modeldata,
                   weights = nlme::varExp()), 
                  error = function(e) NULL)  
           } else {
           model[[i]] <- tryCatch(nlme::gls(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = modeldata,
                  weights = nlme::varExp(form = as.formula(paste("~",covariate)))), 
                  error = function(e) NULL) 
           }
           
          if(length(model[[i]])>0){
           tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]]))  
          tmp$covariate <- modeldata[!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res) 
         if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(res ~ covariate, data=tmp))$coef[2,4]>0.01 &
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(resdev ~ covariate, data=tmp))$coef[2,4]>0.01){
                  covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$tTable[covariate,c(1, 4)]
         } else {
          model[[i]] <- tryCatch(nlme::gls(as.formula(paste(i, 
                  "~baseline", i, "+", confounders[1], "+", confounders[2], 
                  "+", confounders[3], "+", confounders[4], "+", confounders[5], 
                  "+", covariate, sep = "")), data = modeldata,
                weights = nlme::varComb(nlme::varExp(),nlme::varExp(form= as.formula(paste("~",covariate))))), 
                  error = function(e) NULL)  
          if(length(model[[i]])>0){
           tmp <- data.frame(res = resid(model[[i]],type="pearson"), pred = predict(model[[i]]))  
          tmp$covariate <- modeldata[!is.na(modeldata[,i]),covariate]
          tmp$resdev <- abs(tmp$res) 
          if(summary(lm(tmp$res ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(res ~ covariate, data=tmp))$coef[2,4]>0.01 &
       summary(lm(tmp$resdev ~ tmp$pred))$coef[2,4]>0.01 & summary(lm(resdev ~ covariate, data=tmp))$coef[2,4]>0.01){
                  covariate_test[i, c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))] <- summary(model[[i]])$tTable[covariate,c(1, 4)]
                
         }
          } 
         }  
           
           
         }
          
          }}}
            

            covariate_test <- na.omit(covariate_test)
           
          sig <-  rownames(covariate_test)[covariate_test[,paste(covariate, "p", sep = "_")]<p.cutoff]
    covariate_test[, paste(covariate,"p",  "FDR", sep = "_")] <- p.adjust(covariate_test[,paste(covariate,"p", sep = "_")], "fdr")
           
     write.table(covariate_test, paste("ChangeTest_", covariate, "_", select.by, select, ".txt", sep = ""), 
                quote = F, row.names = F, sep = "\t")

   
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
  ggplot2::scale_color_manual(name=group,values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'gray30','red','turquoise4','olivedrab4','purple','darkorange3','lightyellow4','black')[1:length(levels(factor(df[,'gr'])))])+ 
  ggplot2::scale_fill_manual(name=group,values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'cornflowerblue','pink','turquoise2','yellowgreen','plum','darkorange','lightyellow','gray')[1:length(levels(factor(df[,'gr'])))])+
  ggplot2::theme(legend.position = "right",strip.background =  ggplot2::element_rect(color = "white",fill="white"))+
      ggplot2::ylab('Deviance from expected') 
 
                
                if (pdf) {
                  pdf(paste("ChangeTest", covariate, "_", select.by, select, "Covariateplot.pdf", 
                    sep = ""))
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


covariate_test_summary <- covariate_test[,c(paste(covariate, "estimate", sep = "_"),paste(covariate, "p", sep = "_"))]
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

 if(length(intersect(colnames(phylum),sp))==length(sp)){newseqs2 <- paste("XX ", names(seqs1)," ",seqs1," ", abu, sep="")
 } else newseqs2 <- paste("XX Bacteria_", names(seqs1)," ",seqs1," ", abu, sep="")

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
          node_label_size_range=c(0.01,0.013),
          node_label_size_trans="ln",
          #node_label_size=0.03,
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          overlap_avoidance = 0.5, node_label_max=100, make_legend=T, 
          node_size_trans="linear",
          node_size_range=c(0.012,0.05),
          node_color_trans="linear",
          title=paste("Associations with", covariate), title_size=0.03,
          output_file=paste(covariate, "_", select.by, select,"_HeatTree.pdf",sep=""))




#------------

            }
             if(keep.result) return(covariate_test)
        }
    }
}
 
