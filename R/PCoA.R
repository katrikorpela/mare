PCoA <- function(taxonomic.table, meta, readcount.cutoff = 0, group = NULL, group2 = NULL,
    components = c(1, 2), background.variable = NULL, colour.scheme = terrain.colors, legendplace = "topright", 
    select.by = NULL, select = NULL, relative = F, pdf = F, keep.result = F) {
    
   if(Sys.info()[['sysname']] == "Linux") {
  quartz <- function() {X11()}
  }
    if(Sys.info()[['sysname']] == "Windows") {
  quartz <- function() {X11()}
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
    if (relative) {
        taxa <- (taxa/rowSums(taxa))
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
       # groupvar2 <- meta[, group2]
      groupvar2 <- as.numeric(ordered(meta[, group2]))-1
      groupvar <- meta[, group]
    } else groupvar2 <- 0 
    
    pcoa <- vegan::capscale(taxa ~ 1, distance = "bray")
    spcoa <- summary(pcoa)
    
    if (length(group) != 0) {
      if (length(group2) != 0) {
    ado <- vegan::adonis(taxa ~ groupvar + groupvar2, method = "bray")$aov.tab
    } else  ado <- vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab
    }
    
    pc <- data.frame(taxa,spcoa$sites)
    pctest <- as.data.frame(matrix(nrow=ncol(taxa),ncol=4))
    names(pctest) <- c(paste("MDS",components[1],sep=""),
                        paste("MDS",components[2],sep=""),
                         paste("MDS",components[1],"p",sep=""),
                         paste("MDS",components[2],"p",sep=""))
    rownames(pctest)<-colnames(taxa)
    for(i in rownames(pctest)) pctest[i,paste("MDS",components[1],sep="")] <-cor(pc[,paste("MDS",components[1],sep="")],log(pc[,i]+1))
    for(i in rownames(pctest)) pctest[i,paste("MDS",components[2],sep="")]<-cor(pc[,paste("MDS",components[2],sep="")],log(pc[,i]+1))
    for(i in rownames(pctest)) pctest[i,paste("MDS",components[1],"p",sep="")]<-cor.test(pc[,paste("MDS",components[1],sep="")],log(pc[,i]+1))$p.value
    for(i in rownames(pctest)) pctest[i,paste("MDS",components[2],"p",sep="")]<-cor.test(pc[,paste("MDS",components[2],sep="")],log(pc[,i]+1))$p.value

    
    if(nrow(meta)<50) pointsize <- 50/nrow(meta) else pointsize <- 1
    
    op <- par(mar = c(3, 3, 1, 1), xpd = T, cex.lab = 1.5, cex.axis = 1.5, mgp = c(1.5, 0.3, 0), tck = -0.01)
   
    xlimit <- seq(1.5 * min(spcoa$sites[, components[1]]), 1.5 * max(spcoa$sites[, components[1]]), by = 0.01)
    ylimit <- seq(1.5 * min(spcoa$sites[, components[2]]), 1.5 * max(spcoa$sites[, components[2]]), by = 0.01)
        
    if (length(background.variable) != 0) {
        
        grd <- expand.grid(x = xlimit*1.5, y = ylimit*1.5)
        sp::coordinates(grd) <- ~x + y
        sp::gridded(grd) <- TRUE
        temp <- data.frame(spcoa$sites[,  components], meta, taxa)
        names(temp)[c(1, 2)] <- c("MDS1", "MDS2")
        IDW <- gstat::krige(as.formula(paste("log(", background.variable, "+1)~1")), 
            locations = ~MDS1 + MDS2, data = na.omit(temp[, c(background.variable,  "MDS1", "MDS2")]), newdata = grd)
        
        if (pdf) {
            pdf("PCoA.pdf")
            plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.94), 
            xlim = (range(xlimit) * 0.94), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            image(IDW, col = colour.scheme(50), use.raster = T, add = T)
            if (length(group) != 0) {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
                legend(legendplace,inset=c(0.2,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
                mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
                 mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                 } else {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }
            par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0),  mgp = c(0.5, 0, 0), xpd = F, tck = 0.01)
            plot(seq(1, 50, 1), rep(1, 50), col = colour.scheme(50), pch = 15, cex = 10, axes = F, ylab = "", xlab = background.variable)
            axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
            dev.off()
        }
   quartz() 
           plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.94), 
            xlim = (range(xlimit) * 0.94), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            image(IDW, col = colour.scheme(50), use.raster = T, add = T)
            if (length(group) != 0) {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
                legend(legendplace,inset=c(0.2,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
                mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
                 mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                   } else {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }
            par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0),  mgp = c(0.5, 0, 0), xpd = F, tck = 0.01)
            plot(seq(1, 50, 1), rep(1, 50), col = colour.scheme(50), pch = 15, cex = 10, axes = F, ylab = "", xlab = background.variable)
            axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
        
    } else {
      
    if (pdf) {
         pdf("PCoA.pdf")
       plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.94), 
            xlim = (range(xlimit) * 0.94), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            if (length(group) != 0) {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
                legend(legendplace,inset=c(0.2,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
               mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
                 mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               }
                 } else {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }
          axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
      dev.off()
    }
      
 quartz()
       plot(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], type = "n", ylim = (range(ylimit) * 0.94), 
            xlim = (range(xlimit) * 0.94), xlab = paste("Component ", components[1],  " (", 100 * round(spcoa$cont$importance[2, 
            components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", components[2], 
              " (", 100 * round(spcoa$cont$importance[2, components[2]], digits = 2), "%)", sep = ""))
            if (length(group) != 0) {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[, components[1]], pch = 21 + groupvar2,  bg = as.factor(meta[,group]), col="white", cex = pointsize)
                legend(legendplace, legend = levels(as.factor(meta[,group])), title = group, pch = 21, pt.bg = c(1:length(levels(as.factor(meta[,group])))),  bty = "n",col="white")
                if(length(group2)!=0){
                legend(legendplace,inset=c(0.2,0),legend = levels(as.factor(meta[,group2])),  title = group2, pch = 21+unique(groupvar2)[order(unique(groupvar2))], col="black",pt.bg = "white",  bty = "n")
                }
                  mtext(side = 1, line = -3, adj = 0.1, text = paste(group,": ",round(100 * ado[1, 5]), "%  of variation explained", ", p = ", ado[1, 6], sep = ""))
               if(length(group2) != 0) {
                 mtext(side = 1, line = -2, adj = 0.1, text = paste(group2,": ",round(100 * ado[2, 5]), "%  of variation explained", ", p = ", ado[2, 6], sep = ""))
               } } else {
                points(spcoa$sites[, components[2]] ~ spcoa$sites[,components[1]], pch = 21, bg = "black",col="white", cex = pointsize)
            }
          axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
      
    }
    
   if(keep.result)     return(list(pctest,pc))
} 
