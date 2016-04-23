PCoA <- function(taxonomic.table, meta, readcount.cutoff = 0, group = NULL, 
    components = c(1, 2), background.variable = NULL, legendplace = "topright", 
    select.by = NULL, select = NULL, relative = F, pdf = F, quartz = T) {
    
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
        meta[, group] <- meta[, group][drop = T]
        groupvar <- meta[, group]
    }
    
    op <- par(mar = c(3, 3, 1, 1), xpd = T, cex.lab = 1.5, cex.axis = 1.5, 
        mgp = c(1.5, 0.3, 0), tck = -0.01)
    if (length(background.variable) != 0) {
        xlimit <- seq(1.5 * min(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[1]]), 1.5 * max(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[1]]), by = 0.01)
        ylimit <- seq(1.5 * min(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[2]]), 1.5 * max(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[2]]), by = 0.01)
        grd <- expand.grid(x = xlimit, y = ylimit)
        sp::coordinates(grd) <- ~x + y
        sp::gridded(grd) <- TRUE
        temp <- data.frame(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components], meta)
        names(temp)[c(1, 2)] <- c("MDS1", "MDS2")
        IDW <- gstat::krige(as.formula(paste("log(", background.variable, "+1)~1")), 
            locations = ~MDS1 + MDS2, data = na.omit(temp[, c(background.variable, 
                "MDS1", "MDS2")]), newdata = grd)
        if (pdf) {
            pdf("PCoA.pdf")
            plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[1]], type = "n", ylim = (range(ylimit) * 0.94), 
                xlim = (range(xlimit) * 0.94), xlab = paste("Component ", components[1], 
                  " (", 100 * round(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$cont$importance[2, 
                    components[1]], digits = 2), "%)", sep = ""), ylab = paste("Component ", 
                  components[2], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                    1, distance = "bray"))$cont$importance[2, components[2]], 
                    digits = 2), "%)", sep = ""))
            image(IDW, col = heat.colors(50), use.raster = T, add = T)
            if (length(group) != 0) {
                points(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[1]], pch = 21, bg = as.factor(groupvar))
                legend(legendplace, legend = levels(as.factor(groupvar)), title = group, 
                  pch = 21, pt.bg = c(1:length(levels(as.factor(groupvar)))), 
                  bty = "n")
                mtext(side = 1, line = -2, adj = 0.1, text = paste(round(100 * 
                  vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab[1, 
                    5]), "%  of variation explained", ", p = ", vegan::adonis(taxa ~ 
                  groupvar, method = "bray")$aov.tab[1, 6], sep = ""))
            } else {
                points(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[1]], pch = 21, bg = "black")
            }
            par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0), 
                mgp = c(0.5, 0, 0), xpd = F, tck = 0.01)
            plot(seq(1, 50, 1), rep(1, 50), col = heat.colors(50), pch = 15, 
                cex = 10, axes = F, ylab = "", xlab = background.variable)
            axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
            dev.off()
        }
        if (quartz) quartz()
        plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
            components[1]], type = "n", ylim = (range(ylimit) * 0.94), xlim = (range(xlimit) * 
            0.94), xlab = paste("Component ", components[1], " (", 100 * round(summary(vegan::capscale(taxa ~ 
            1, distance = "bray"))$cont$importance[2, components[1]], digits = 2), 
            "%)", sep = ""), ylab = paste("Component ", components[2], " (", 
            100 * round(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$cont$importance[2, 
                components[2]], digits = 2), "%)", sep = ""))
        image(IDW, col = heat.colors(50), use.raster = T, add = T)
        if (length(group) != 0) {
            points(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[1]], pch = 21, bg = as.factor(groupvar))
            legend(legendplace, legend = levels(as.factor(groupvar)), title = group, 
                pch = 21, pt.bg = c(1:length(levels(as.factor(groupvar)))), 
                bty = "n")
            mtext(side = 1, line = -2, adj = 0.1, text = paste(round(100 * 
                vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab[1, 
                  5]), "%  of variation explained", ", p = ", vegan::adonis(taxa ~ 
                groupvar, method = "bray")$aov.tab[1, 6], sep = ""))
        } else {
            points(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[1]], pch = 21, bg = "black")
        }
        par(fig = c(0.1, 0.9, 0.88, 0.99), new = T, mar = c(2, 0, 0, 0), mgp = c(0.5, 
            0, 0), xpd = F, tck = 0.01)
        plot(seq(1, 50, 1), rep(1, 50), col = heat.colors(50), pch = 15, cex = 10, 
            axes = F, ylab = "", xlab = background.variable)
        axis(side = 1, at = c(3, 47), labels = c("Low", "High"))
        
    } else {
        if (length(group) != 0) {
            if (pdf) {
                pdf("PCoA.pdf")
                plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[1]], pch = 21, bg = as.factor(groupvar), xlab = paste("Component ", 
                  components[1], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                    1, distance = "bray"))$cont$importance[2, components[1]], 
                    digits = 2), "%)", sep = ""), ylab = paste("Component ", 
                  components[2], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                    1, distance = "bray"))$cont$importance[2, components[2]], 
                    digits = 2), "%)", sep = ""))
                legend(legendplace, legend = levels(as.factor(groupvar)), title = group, 
                  pch = 21, pt.bg = c(1:length(levels(as.factor(groupvar)))), 
                  bty = "n")
                mtext(side = 1, line = -2, adj = 0.1, text = paste(round(100 * 
                  vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab[1, 
                    5]), "%  of variation explained", ", p = ", vegan::adonis(taxa ~ 
                  groupvar, method = "bray")$aov.tab[1, 6], sep = ""))
                dev.off()
            }
            if (quartz) quartz()
            plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[1]], pch = 21, bg = as.factor(groupvar), xlab = paste("Component ", 
                components[1], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                  1, distance = "bray"))$cont$importance[2, components[1]], 
                  digits = 2), "%)", sep = ""), ylab = paste("Component ", 
                components[2], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                  1, distance = "bray"))$cont$importance[2, components[2]], 
                  digits = 2), "%)", sep = ""))
            legend(legendplace, legend = levels(as.factor(groupvar)), title = group, 
                pch = 21, pt.bg = c(1:length(levels(as.factor(groupvar)))), 
                bty = "n")
            mtext(side = 1, line = -2, adj = 0.1, text = paste(round(100 * 
                vegan::adonis(taxa ~ groupvar, method = "bray")$aov.tab[1, 
                  5]), "%  of variation explained", ", p = ", vegan::adonis(taxa ~ 
                groupvar, method = "bray")$aov.tab[1, 6], sep = ""))
            
        } else {
            if (pdf) {
                pdf("PCoA.pdf")
                plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                  components[1]], pch = 21, bg = "black", xlab = paste("Component ", 
                  components[1], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                    1, distance = "bray"))$cont$importance[2, components[1]], 
                    digits = 2), "%)", sep = ""), ylab = paste("Component ", 
                  components[2], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                    1, distance = "bray"))$cont$importance[2, components[2]], 
                    digits = 2), "%)", sep = ""))
                dev.off()
            }
            if (quartz) quartz()
            plot(summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[2]] ~ summary(vegan::capscale(taxa ~ 1, distance = "bray"))$sites[, 
                components[1]], pch = 21, bg = "black", xlab = paste("Component ", 
                components[1], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                  1, distance = "bray"))$cont$importance[2, components[1]], 
                  digits = 2), "%)", sep = ""), ylab = paste("Component ", 
                components[2], " (", 100 * round(summary(vegan::capscale(taxa ~ 
                  1, distance = "bray"))$cont$importance[2, components[2]], 
                  digits = 2), "%)", sep = ""))
            
        }
    }
} 
