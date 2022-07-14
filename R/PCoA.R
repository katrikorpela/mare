PCoA <- function(taxonomic.table, meta, readcount.cutoff = 0, group = NULL, group2 = NULL,
    components = c(1, 2), background.variable = NULL, colour.scheme = "terrain.colors", 
    ellipse = F, hull = F, spider = F, dots = T, legendplace = "topright", legend =T,
    select.by = NULL, select = NULL, relative = F, pdf = F, keep.result = F,
    subjectID = NULL, time = NULL, distance = "bray", constrain = NULL) {
    
   if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
}

if(colour.scheme[1]=="terrain.colors") colour.scheme <-terrain.colors(n=50)
if(colour.scheme[1]=="topo.colors") colour.scheme <-topo.colors(n=50)
if(colour.scheme[1]=="cm.colors") colour.scheme <- cm.colors(n=50)
if(colour.scheme[1]=="heat.colors") colour.scheme <- heat.colors(n=50)

    peardist <- function(x,distance){
      as.dist(1-cor(t(log(x + min(x[x>0])))))
    }
    
    meta <- read.delim(meta)
    taxa <- read.delim(taxonomic.table)
    
    taxa <- taxa[meta$ReadCount > readcount.cutoff, ]
    meta <- meta[meta$ReadCount > readcount.cutoff, ]
    
    
    if (length(select.by) != 0) {
        meta$selection <- meta[, select.by]
        taxa <- taxa[meta$selection == select, ]
        meta <- meta[meta$selection == select, ]
    }
    
    meta <- meta[rowSums(taxa)>0&!is.na(rowSums(taxa)),]
    taxa <- taxa[rowSums(taxa)>0&!is.na(rowSums(taxa)),]
    
    taxa <- taxa[,names(colSums(taxa)[colSums(taxa)>0])]
    
    if (relative) {
        taxa <- ((taxa)/rowSums(taxa))
    }
    
    if (length(group) != 0) {
      taxa <- taxa[!is.na(meta[,group]),]
      meta <- meta[!is.na(meta[,group]),]
      meta[, group] <- meta[, group][drop = T]
      groupvar <- meta[, group]
    }
    
    if (length(group2) != 0) {
      taxa <- taxa[!is.na(meta[,group2]),]
      meta <- meta[!is.na(meta[,group2]),]
      meta[, group2] <- meta[, group2][drop = T]
      groupvar2 <- as.numeric(ordered(meta[, group2]))-1
    } else groupvar2 <- 0 
  
    if(length(constrain)>0) {
      meta$constrain <- meta[,constrain]
    taxa <- taxa[!is.na(meta$constrain),]
    meta <- meta[!is.na(meta$constrain),]
  if(distance == "pearson") pcoa <- vegan::capscale(taxa ~ constrain,dfun = "peardist", data=meta)
  if(distance == "euc") pcoa <- vegan::capscale(log(taxa+ min(taxa[taxa>0])/10) ~ constrain, distance = "euc",data=meta)
 if(distance == "bray") pcoa <- vegan::capscale(taxa ~ constrain, distance = "bray",data=meta)
    } else {  
      if(distance == "pearson") pcoa <- vegan::capscale(taxa ~ 1,dfun = "peardist")
     if(distance == "euc")  pcoa <- vegan::capscale(log(taxa+ min(taxa[taxa>0])/10) ~ 1, distance = "euc")
      if(distance == "bray") pcoa <- vegan::capscale(taxa ~ 1, distance = "bray")
    }
    
    spcoa <- summary(pcoa)
  
    if (length(group) != 0) {
      if (length(group2) != 0) {
    if(distance == "pearson") ado <- vegan::adonis(peardist(taxa) ~ groupvar + groupvar2, method = "euc")$aov.tab
    if(distance == "euc") ado <- vegan::adonis(log(taxa+ min(taxa[taxa>0])/10) ~ groupvar + groupvar2, method = "euc")$aov.tab
    if(distance == "bray") ado <- vegan::adonis(taxa ~ groupvar + groupvar2, method = "bray")$aov.tab
    } else {
     if(distance == "pearson")  ado <- vegan::adonis(peardist(taxa) ~ groupvar, method = "euc")$aov.tab
     if(distance == "euc")  ado <- vegan::adonis(log(taxa+ min(taxa[taxa>0])/10) ~ groupvar, method = "euc")$aov.tab
     if(distance == "bray")  ado <- vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab
    }
      }
    
    pc <- data.frame(taxa,spcoa$sites)
   
    
    if(keep.result){
      pcscores <- data.frame(spcoa$species) 
      pctest <- as.data.frame(matrix(nrow=ncol(taxa),ncol=4))
    names(pctest) <- c(paste("cor_MDS",components[1],sep=""),
                        paste("cor_MDS",components[2],sep=""),
                         paste("cor_MDS",components[1],"p",sep=""),
                         paste("cor_MDS",components[2],"p",sep=""))
    rownames(pctest)<-colnames(taxa)
    for(i in rownames(pctest)) pctest[i,paste("cor_MDS",components[1],sep="")] <-cor(spcoa$sites[,components[1]],log(pc[,i]+1))
    for(i in rownames(pctest)) pctest[i,paste("cor_MDS",components[2],sep="")]<-cor(spcoa$sites[,components[2]],log(pc[,i]+1))
    for(i in rownames(pctest)) pctest[i,paste("cor_MDS",components[1],"p",sep="")]<-cor.test(spcoa$sites[,components[1]],log(pc[,i]+1))$p.value
    for(i in rownames(pctest)) pctest[i,paste("cor_MDS",components[2],"p",sep="")]<-cor.test(spcoa$sites[,components[2]],log(pc[,i]+1))$p.value
pcscores <- data.frame(pcscores,pctest)
    }
    temp <- data.frame(spcoa$sites[,  components], meta, taxa)
    names(temp)[c(1, 2)] <- c("MDS1", "MDS2")
    
    #if(nrow(meta)<30) {pointsize <- 30/nrow(meta)} else pointsize <- 1
    pointsize <- 100/nrow(meta)
    if(pointsize<1 & length(group)>0) pointsize <- 1

    
    #op <- par(mar = c(3, 3, 1, 1), xpd = T, cex.lab = 1.5, cex.axis = 1.5, mgp = c(1.5, 0.3, 0), tck = -0.01)
   
    xlimit <- seq(1.5 * min(spcoa$sites[, components[1]]), 1.5 * max(spcoa$sites[, components[1]]), by = 0.01)
    ylimit <- seq(1.5 * min(spcoa$sites[, components[2]]), 1.5 * max(spcoa$sites[, components[2]]), by = 0.01)
        
    if (length(background.variable) != 0) {
        
        grd <- expand.grid(x = xlimit*1.5, y = ylimit*1.5)
        sp::coordinates(grd) <- ~x + y
        sp::gridded(grd) <- TRUE
        
        IDW <- gstat::krige(as.formula(paste("log(", background.variable, "+1)~1")), 
            locations = ~MDS1 + MDS2, data = na.omit(temp[, c(background.variable,  "MDS1", "MDS2")]), newdata = grd)
        
        if (pdf) {
             pdf(paste("PCoA",distance, group, group2,".pdf", sep="_"))
           plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.7), 
            xlim = (range(xlimit) * 0.7), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            image(IDW, col = colour.scheme, use.raster = T, add = T,alpha=c(1:0))
                    
          if(spider & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordispider(pcoa,groups=as.factor(meta[,group]),
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,lwd=3)
            }
          }
          if(hull & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordihull(pcoa,groups=as.factor(meta[,group]),draw="polygon",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          }
          if(ellipse & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordiellipse(pcoa,groups=as.factor(meta[,group]),draw="polygon",kind="sd",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          } 
            if (length(group) != 0) {
             if(dots)   points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
              if(legend) legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
               if(legend) legend(legendplace,inset=c(0.5,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
              if(legend)  mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
              if(legend)   mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                   } else {
                if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
                   }

        
                   if(length(subjectID)>0){  
           temp$subject <- temp[,subjectID]
           temp$subject <- temp$subject[drop=T]
           temp$time <- temp[,time]
            for(i in  names(table(temp$subject)[table(temp$subject)>1])) {
                temp$time[temp$subject==i] <- order(temp$time[temp$subject==i])
              for(j in temp$time[temp$subject==i][order(temp$time[temp$subject==i])][-1]){
              arrows(x0=temp[temp$time==(j-1)&temp$subject==i,"MDS1"],x1=temp[temp$time==j&temp$subject==i,"MDS1"],
             y0=temp[temp$time==(j-1)&temp$subject==i,"MDS2"],y1=temp[temp$time==j&temp$subject==i,"MDS2"],length=0.1)
              }
            }
}
           LL <- pretty(temp[,background.variable],n=length(colour.scheme))
            par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0),  mgp = c(0.75, 0, 0), xpd = F, tck = 0.01,cex.lab=0.8)
            plot(seq(1, length(colour.scheme), 1), rep(1, length(colour.scheme)), col = colour.scheme, pch = 15, cex = 10, axes = F, ylab = "", xlab =background.variable)
            axis(side = 1, at = seq(1,length(colour.scheme)+0.5,length(colour.scheme)/length(LL)), labels=LL)#labels = c(-1*round(exp(seq(5,1,-1))),0,round(exp(seq(1,5,1))))) 
        
            dev.off()
        }
   quartz() 
           plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.7), 
            xlim = (range(xlimit) * 0.7), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            image(IDW, col = colour.scheme, use.raster = T, add = T)
              
                     if(spider & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordispider(pcoa,groups=as.factor(meta[,group]),
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,lwd=3)
            }
          }
          if(hull & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordihull(pcoa,groups=as.factor(meta[,group]),draw="polygon",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          }
          if(ellipse & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordiellipse(pcoa,groups=as.factor(meta[,group]),draw="polygon",kind="sd",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          } 
            if (length(group) != 0) {
                if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
             if(legend) legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
              if(legend)  legend(legendplace,inset=c(0.5,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
              if(legend)  mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
              if(legend)   mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                   } else {
                 if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
                   }
            

        
                   if(length(subjectID)>0){  
           temp$subject <- temp[,subjectID]
           temp$subject <- temp$subject[drop=T]
           temp$time <- temp[,time]
            for(i in  names(table(temp$subject)[table(temp$subject)>1])) {
                temp$time[temp$subject==i] <- order(temp$time[temp$subject==i])
              for(j in temp$time[temp$subject==i][order(temp$time[temp$subject==i])][-1]){
              arrows(x0=temp[temp$time==(j-1)&temp$subject==i,"MDS1"],x1=temp[temp$time==j&temp$subject==i,"MDS1"],
             y0=temp[temp$time==(j-1)&temp$subject==i,"MDS2"],y1=temp[temp$time==j&temp$subject==i,"MDS2"],length=0.1)
              }
            }
}
           LL <- pretty(temp[,background.variable],n=length(colour.scheme))
            par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0),  mgp = c(0.75, 0, 0), xpd = F, tck = 0.01,cex.lab=0.8)
            plot(seq(1, length(colour.scheme), 1), rep(1, length(colour.scheme)), col = colour.scheme, pch = 15, cex = 10, axes = F, ylab = "", xlab =background.variable)
            axis(side = 1, at = seq(1,length(colour.scheme)+0.5,length(colour.scheme)/length(LL)), labels=LL)#labels = c(-1*round(exp(seq(5,1,-1))),0,round(exp(seq(1,5,1))))) 
        
    } else {
      
    if (pdf) {
         pdf(paste("PCoA", distance, group, group2,".pdf", sep="_"))
       plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.7), 
            xlim = (range(xlimit) * 0.7), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
                  if(spider & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordispider(pcoa,groups=as.factor(meta[,group]),
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,lwd=3)
            }
          }
          if(hull & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordihull(pcoa,groups=as.factor(meta[,group]),draw="polygon",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          }
          if(ellipse & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordiellipse(pcoa,groups=as.factor(meta[,group]),draw="polygon",kind="sd",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          }   
       if (length(group) != 0) {
                if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                if(legend)legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
               if(legend) legend(legendplace,inset=c(0.5,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
              if(legend) mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
               if(legend)  mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                 } else {
                if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }

          
                   if(length(subjectID)>0){  
           temp$subject <- temp[,subjectID]
           temp$subject <- temp$subject[drop=T]
           temp$time <- temp[,time]
            for(i in  names(table(temp$subject)[table(temp$subject)>1])) {
                temp$time[temp$subject==i] <- order(temp$time[temp$subject==i])
              for(j in temp$time[temp$subject==i][order(temp$time[temp$subject==i])][-1]){
              arrows(x0=temp[temp$time==(j-1)&temp$subject==i,"MDS1"],x1=temp[temp$time==j&temp$subject==i,"MDS1"],
             y0=temp[temp$time==(j-1)&temp$subject==i,"MDS2"],y1=temp[temp$time==j&temp$subject==i,"MDS2"],length=0.1)
              }
            }
}
      dev.off()
    }
      
 quartz()
       plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.7), 
            xlim = (range(xlimit) * 0.7), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
              if(spider & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordispider(pcoa,groups=as.factor(meta[,group]),
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,lwd=3)
            }
          }
          if(hull & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordihull(pcoa,groups=as.factor(meta[,group]),draw="polygon",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          }
          if(ellipse & length(group) != 0){
            for(i in seq_along(levels(as.factor(meta[,group])))){  
            vegan::ordiellipse(pcoa,groups=as.factor(meta[,group]),draw="polygon",kind="sd",
                                 show.groups=levels(as.factor(meta[,group]))[i],choices=components,col=i,border=i,lwd=3)
            }
          } 
       if (length(group) != 0) {
               if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                if(legend) legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
                if(legend) legend(legendplace,inset=c(0.5,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
              if(legend)    mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
              if(legend)   mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               } } else {
               if(dots) points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }
  
   
                   if(length(subjectID)>0){  
           temp$subject <- temp[,subjectID]
           temp$subject <- temp$subject[drop=T]
           temp$time <- temp[,time]
            for(i in  names(table(temp$subject)[table(temp$subject)>1])) {
              temp$time[temp$subject==i][order(temp$time[temp$subject==i])] <- order(temp$time[temp$subject==i][order(temp$time[temp$subject==i])])
              for(j in temp$time[temp$subject==i][order(temp$time[temp$subject==i])][-1]){
              arrows(x0=temp[temp$time==(j-1)&temp$subject==i,"MDS1"],x1=temp[temp$time==j&temp$subject==i,"MDS1"],
             y0=temp[temp$time==(j-1)&temp$subject==i,"MDS2"],y1=temp[temp$time==j&temp$subject==i,"MDS2"],length=0.1)
              }
            }
}
    }
    
   if(keep.result)     return(list(pcscores,pc))
} 
