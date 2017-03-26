Organise <- function(meta, sample.names, otutable = list.files(pattern = "_otutable.txt$")[1], 
    taxonomic.tables = list.files(pattern = "_table.txt$"), subject.ID = NULL, 
    time = NULL) {
  
    metadata <- read.delim(meta)
    rownames(metadata) <- metadata[, sample.names]
    readtable <- data.frame(t(read.table(list.files(pattern = "annotated_read_table.txt$")[1], header = T, row.names = 1, check.names = F)))
    readtable <- readtable[rownames(metadata), ]
    genera <- sapply(colnames(readtable), function(x) strsplit(x, split="_")[[1]][5])
    generalist <- list()
     for(i in na.omit(unique(genera))) {
       generalist[[i]] <- cbind(readtable[,names(genera[genera==i&!is.na(genera)])],rep(0,nrow(readtable)))
       metadata[,paste(i,"richness",sep="_")] <- vegan::specnumber(generalist[[i]])
       }
   
   for (i in seq_along(taxonomic.tables)) {
        write.table(data.frame(t(read.table(taxonomic.tables[i], header = T, 
            row.names = "taxon", check.names = F)))[rownames(metadata), ], 
            paste("organised", taxonomic.tables[i], sep = "_"), quote = F, 
            row.names = T, sep = "\t")
   }
   
    otu <- data.frame(t(read.table(otutable, header = T, row.names = "OTUId", check.names = F)))
    otu <- otu[rownames(metadata), ]
    metadata$Richness <- vegan::specnumber(otu)
    metadata$Diversity <- vegan::diversity(otu, "inv")
    
    reads <- read.table(list.files(pattern="readnumbers.txt"), header = T, row.names = "sample", 
        check.names = F)
    reads <- reads[rownames(metadata), ]
    metadata$ReadCount <- reads$processed_reads
    metadata$ReadCountClass <- Hmisc::cut2(metadata$ReadCount, g = 4, digits = 1)
    
    pdf("RichnessReadcount.pdf")
    plot(metadata$Richness ~ metadata$ReadCount, ylab = "OTU Richness", xlab = "Number of reads", 
        pch = 21, bg="skyblue")
    dev.off()
    
    write.table(otu, "organised_otutable.txt", quote = F, row.names = T, sep = "\t")
    
    
    if (length(time) != 0 & length(subject.ID) != 0) {
        time <- as.numeric(metadata[, time])
        reltaxa <- readtable/metadata$ReadCount
        for (i in names(table(metadata[, subject.ID])[table(metadata[, subject.ID]) > 1])) {
            time[metadata[, subject.ID] == i] <- order(time[metadata[, subject.ID] == i])
            for (j in time[metadata[, subject.ID] == i][order(time[metadata[, subject.ID] == i])][-1]) {
                metadata$IntraindividualPearson[metadata[, subject.ID] == i &  time == j] <- cor((log(t(reltaxa[metadata[, subject.ID] == 
                  i & time == j, ] + 1e-11))), (log(t(reltaxa[metadata[, subject.ID] ==  i & time == (j - 1), ] + 1e-11))))
                metadata$IntraindividualBrayCurtis[metadata[, subject.ID] == 
                  i & time == j] <- vegan::vegdist(rbind(reltaxa[metadata[, 
                  subject.ID] == i & time == j, ], reltaxa[metadata[, subject.ID] == 
                  i & time == (j - 1), ]))
            }
        }
    }
    write.table(metadata, meta, quote = F, row.names = F, sep = "\t")
} 
