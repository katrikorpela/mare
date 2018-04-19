ProcessReads <- function(forward.reads = NULL, 
                         reverse.reads = NULL, 
                         name.separator = "_", 
    usearch.path, forward.primer = NULL, reverse.primer = NULL,
    trim.beginning = 0, merge = F, 
    min.merged.length = NULL, min.overlap = NULL, max.mismatches = NULL, 
    pool = F, truncate = T, 
    filter = T, trunc.quality = 2, max.expected.error = 1, 
    readlength = NULL, min.read.abundance = 0.00001, folder.name = "") {
  
  name.separator2 <- name.separator
  if(name.separator==".") name.separator <- "[.]"
  
  if (length(list.files(pattern = paste(folder.name, "ProcessedReads", sep = ""))) == 0) {
       system(paste("mkdir ", folder.name, "ProcessedReads", sep = ""))
  } else {
    system(paste("rm -r ", folder.name, "ProcessedReads", sep = "")) 
    system(paste("mkdir ", folder.name, "ProcessedReads", sep = ""))
   }
    
    wd <- paste(getwd(), "/", folder.name, "ProcessedReads/", sep = "")
    
    readnumbers <- data.frame(sample = sapply(as.character(forward.reads), 
                                              function(x) strsplit(x, split = name.separator,fixed=F)[[1]][1]), 
                              raw_reads = sapply(forward.reads, function(x) Biostrings::fastq.geometry(x)[1])) 
   
    if(pool) readnumbers$raw_reads <-  readnumbers$raw_reads * 2
    
    if (length(forward.primer) != 0) {
        for (i in forward.reads) {
            ShortRead::writeFastq(object = Biostrings::trimLRPatterns(Lpattern = forward.primer, 
                subject = ShortRead::readFastq(i), max.Lmismatch = 3, with.Lindels = T, Lfixed = F), 
                file = paste(wd, sapply(i, function(x) strsplit(x, name.separator, 
                  fixed = F)[[1]][1]), name.separator2 ,"fwd.fastq", sep = ""), mode = "w", 
                full = F, compress = F, qualityType = "Auto")
        }
    } else {
      for (i in seq_along(forward.reads)) {
         if(grepl(x=forward.reads[i], pattern=".gz$")){
         system(paste("cp ", forward.reads[i], " ", wd, forward.reads[i], sep = ""))
         system(paste("gunzip ",wd, forward.reads[i], sep=""))  
         system(paste("mv ", wd, strsplit(forward.reads[i],split=".gz")[[1]][1], " ", wd, strsplit(forward.reads[i],split=name.separator, fixed = F)[[1]][1],  
                    name.separator2 ,"fwd.fastq", sep = ""))
      } else{
        system(paste("cp ", forward.reads[i], " ", wd, forward.reads[i], sep = ""))
         system(paste("mv ", wd, forward.reads[i], " ", wd, strsplit(forward.reads[i],split=name.separator, fixed = F)[[1]][1],  
                    name.separator2 ,"fwd.fastq", sep = ""))
      }
      }
      }
  
    trimmed_fwd <- paste(wd, list.files(path = wd, pattern = "fwd.fastq"),  sep = "")
    
      for (i in seq_along(trimmed_fwd)) {
     system(paste(usearch.path, " -fastx_truncate ", trimmed_fwd[i], 
                " -stripleft ", trim.beginning, " -fastqout ",trimmed_fwd[i],"2", sep = ""))
      system(paste("mv ", trimmed_fwd[i],"2 ",trimmed_fwd[i],sep=""))
      #system(paste("rm ", trimmed_fwd[i],"2",sep=""))
      }
    
    if (length(reverse.reads) != 0) {
        if (length(reverse.primer) != 0) {
            for (i in reverse.reads) {
                ShortRead::writeFastq(Biostrings::trimLRPatterns(Lpattern = reverse.primer, 
                  subject = ShortRead::readFastq(i), max.Lmismatch = 3, with.Lindels = T), 
                  file = paste(wd, sapply(i, function(x) strsplit(x, name.separator, 
                    fixed = F)[[1]][1]), name.separator2 ,"rev.fastq", sep = ""), mode = "w", 
                  full = F, compress = F, qualityType = "Auto")
            }
        } else {
      for (i in seq_along(reverse.reads)) {
         if(grepl(x=reverse.reads[i], pattern=".gz$")){
         system(paste("cp ", reverse.reads[i], " ", wd, reverse.reads[i], sep = ""))
         system(paste("gunzip ",wd, reverse.reads[i], sep=""))  
         system(paste("mv ", wd, strsplit(reverse.reads[i],split=".gz")[[1]][1], " ", wd, strsplit(reverse.reads[i],split=name.separator, fixed = F)[[1]][1],  
                    name.separator2 ,"rev.fastq", sep = ""))
      } else{
        system(paste("cp ", reverse.reads[i], " ", wd, reverse.reads[i], sep = ""))
         system(paste("mv ", wd, reverse.reads[i], " ", wd, strsplit(reverse.reads[i],split=name.separator, fixed = F)[[1]][1],  
                    name.separator2 ,"rev.fastq", sep = ""))
      }
      }
      }
      trimmed_rev <- paste(wd, list.files(path = wd, pattern = "rev.fastq"),  sep = "")
      
      for (i in seq_along(trimmed_rev)) {
     system(paste(usearch.path, " -fastx_truncate ", trimmed_rev[i], 
                " -stripleft ", trim.beginning, " -fastqout ",trimmed_rev[i],"2", sep = ""))
      system(paste("mv ", trimmed_rev[i],"2 ",trimmed_rev[i],sep=""))
      system(paste("rm ", trimmed_rev[i],"2",sep=""))
      }
    }
    
    reads <- paste(wd, list.files(path = wd, pattern = "fwd.fastq"), sep = "")
    
    if (merge) {
        for (i in seq_along(trimmed_fwd)) {
            system(paste(usearch.path, " -fastq_mergepairs ", trimmed_fwd[i], 
                " -reverse ", trimmed_rev[i], " -fastq_minmergelen ", min.merged.length, 
                " -fastq_minovlen ", min.overlap, " -fastq_maxdiffs ", max.mismatches, 
                " -fastqout ", strsplit(trimmed_fwd[i], name.separator , fixed = F)[[1]][1], 
                name.separator2,"merged.fastq", sep = ""))
          system(paste("rm ", trimmed_fwd[i], sep = ""))
          system(paste("rm ", trimmed_rev[i], sep = ""))
        }
        reads <- paste(wd, list.files(path = wd, pattern = "merged.fastq"), sep = "")
    
    } else if (pool) {
       for (i in seq_along(trimmed_fwd)) {
        system(paste("cat ", trimmed_fwd[i], " ", trimmed_rev[i], " > ", strsplit(trimmed_fwd[i], name.separator , fixed = F)[[1]][1], 
                name.separator2,"pooled.fastq", sep = "")) 
         system(paste("rm ", trimmed_fwd[i], sep = ""))
        system(paste("rm ", trimmed_rev[i], sep = ""))
       } 
      reads <- paste(wd, list.files(path = wd, pattern = "pooled.fastq"), sep = "")
    }
      
     
      if (truncate) {
        for (i in seq_along(reads)) {
            system(paste(usearch.path, " -fastx_truncate ", reads[i], 
                " -trunclen ", readlength, " -fastqout ", strsplit(reads[i], 
                  name.separator, fixed = T)[[1]][1], name.separator2 ,"truncated.fastq", sep = ""))
           system(paste("rm ", reads[i], sep = ""))
        }
        reads <- paste(wd, list.files(path = wd, pattern = "truncated.fastq"),  sep = "")
     }
    
    if (filter) {
        for (i in seq_along(reads)) {
            system(paste(usearch.path, " -fastq_filter ", reads[i], " -fastq_truncqual ", 
                trunc.quality, " -fastq_maxee ", max.expected.error, " -fastq_minlen ", 
                readlength, 
                " -fastaout ", strsplit(reads[i], split = ".fastq",fixed=T)[[1]][1], 
                name.separator2, "filtered.fasta", sep = ""))
          system(paste("rm ", reads[i], sep = ""))
        }
    
    } else for (i in seq_along(reads)) {
        ShortRead::writeFasta(ShortRead::readFastq(reads[i]), file = paste(strsplit(reads[i], 
            split = "[.]", fixed = F)[[1]][1], "fasta", sep = "."))
     system(paste("rm ", reads[i], sep = ""))
    }
    reads <- paste(wd, list.files(path = wd, pattern = "*.fasta$"), sep = "")

    for (i in seq_along(reads)) {
        system(paste(usearch.path, " -fastx_truncate ", reads[i], " -trunclen ", 
            readlength, " -fastaout ", paste(strsplit(reads[i], split = ".fasta", fixed = F), 
                "equallength.fasta", sep = name.separator2), sep = ""))
      system(paste("rm ", reads[i], sep = ""))
    }
  
    truncatedreads <- paste(wd, list.files(path = wd, pattern = "equallength.fasta$"), sep = "")
    
    for (i in seq_along(truncatedreads)) {
        system(paste("sed '-es/^>\\(.*\\)/>\\1;barcodelabel=", strsplit(truncatedreads[i], 
            split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2], 
            ";/' <", truncatedreads[i], " > ", paste(strsplit(truncatedreads[i], 
                split = ".fasta", fixed = F), "labelled.fasta", sep = name.separator2), sep = ""))
     system(paste("rm ", truncatedreads[i], sep = "")) 
    }
   
    labelledreads <- paste(wd, list.files(path = wd, pattern = "*labelled.fasta$"), sep = "")
    fl.uparse.string <- ""
    for (i in seq_along(labelledreads)) {
        fl.uparse.string <- paste(fl.uparse.string, labelledreads[i], sep = " ")
    }

    system(paste("cat ", fl.uparse.string, " > ", wd, "all.fasta", sep = ""))
    
    dnaSet = Biostrings::readDNAStringSet(paste(wd, "all.fasta", sep = ""), format = "fasta")
    if (!is(dnaSet, "DNAStringSet")) {
        dnaSet <- Biostrings::DNAStringSet(dnaSet)
    }
    if (is.null(names(dnaSet))) {
        message("No names attribute found in dnaSet object...", "using artifically generated names")
        names(dnaSet) <- paste("read", 1:length(dnaSet), sep = "-")
    }
    dnaSet <- dnaSet[order(dnaSet)]#tämä tarvitaan??!
    counts <- BiocGenerics::table(dnaSet)
    dnaSet <- unique(dnaSet)
    names(dnaSet) <- paste0(paste("read", 1:length(names(dnaSet)), sep = "_"), 
        ";size=", as.integer(counts)) 
    dnaSet <- dnaSet[rev(order(counts))]
    dnaSet <- dnaSet[as.vector(counts[rev(order(counts))] > (min.read.abundance * sum(counts)))]
    Biostrings::writeXStringSet(dnaSet, paste(wd, "dereplicated.fasta", sep = ""), format = "fasta")
 
    system(paste(usearch.path, " -uchime_denovo ", wd,  "dereplicated.fasta", " -nonchimeras ", wd,
            "nonchimeric.fasta", sep = ""))
    system(paste(usearch.path, " -cluster_otus ", wd, "nonchimeric.fasta", 
        " -otus ", wd, "otu.fasta", " -otu_radius_pct 3", 
        " -relabel OTU_", " -sizein ", " -sizeout -uparseout ", wd, "OTUmapping.txt", sep = ""))
  
    system(paste("echo ", shQuote('V1\tV2\tV3\tV4\tV5\tV6'), " > ", wd, "otumappingheader.txt", sep=""))
    system(paste("cat ", wd, "OTUmapping.txt >> ", wd, "otumappingheader.txt",sep=""))
    otumapping <- read.delim(paste(wd,"otumappingheader.txt",sep=""),header=T)
    otumapping$V6 <- as.character(otumapping$V6)
    otumapping$V5 <- as.character(otumapping$V5)
    otumapping$V6[otumapping$V6==""] <-  otumapping$V5[otumapping$V6==""] 
    chimera <- otumapping$V1[otumapping$V2=="chimera"]
   
   derep <- Biostrings::readDNAStringSet(paste(wd,"nonchimeric.fasta",sep=""))
   result <- data.frame(OTUId = setdiff(names(derep),chimera))
   Biostrings::writeXStringSet(derep[setdiff(names(derep),chimera)], paste(wd, "finalreads.fasta", sep = ""), format = "fasta")
   
   for(j in result$OTUId) result[result$OTUId==j,"OTU"] <- otumapping$V6[otumapping$V1==j]
   dict <- Biostrings::PDict(derep[setdiff(names(derep),chimera)], tb.start = (readlength-1))
    for (i in seq_along(labelledreads)) {
        seq <- Biostrings::readDNAStringSet(labelledreads[i])
        result[, strsplit(labelledreads[i], split =  paste("ProcessedReads/|",name.separator,sep=""),fixed=F)[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict,seq))
    }
    
   result2 <- result
   rownames(result2)<-result2$OTUId  
   result2 <- result2[,-c(1,2)]
   readinfo <- data.frame(read=rownames(result2), readtotals=rowSums(result2)/sum(result2),OTU=result$OTU)
   sampleinfo <- result2/colSums(result2)
   abundant <- rownames(readinfo)[readinfo$readtotals>0.01]
   result2[abundant,][sampleinfo[abundant,]<0.001] <- 0
   result2 <- cbind(OTUId = rownames(result2),result2)
   
   if (length(list.files(pattern = paste(folder.name, "TaxonomicTables", sep = ""))) == 0) {
        system(paste("mkdir ", folder.name, "TaxonomicTables", sep = ""))
   }
   
   write.table(result2, paste(getwd(), "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t") 

    rownames(readnumbers) <- readnumbers$sample
    sample = sapply(as.character(labelledreads),  function(x) gsub(pattern=paste(wd,folder.name,"ProcessedReads/",sep=""),
                                                               strsplit(x, split = name.separator,fixed=F)[[1]][1],
                                              replacement="")) 
    for(i in readnumbers$sample) readnumbers[readnumbers$sample==i,"processed_reads"] <- sum(result2[,i])
    write.table(readnumbers, paste(folder.name, "readnumbers.txt",sep="_"), quote = F, row.names = F, sep = "\t")
    readnumbers$raw_reads <- readnumbers$raw_reads - readnumbers$processed_reads
    
    wi <- nrow(readnumbers)/4
    if(wi < 4) wi <- 4
    
    pdf(paste(folder.name,"Readnumbers.pdf",sep="_"), width = wi)
    par(xpd = T, mar = c(max(nchar(rownames(readnumbers)))/2, 6, 4, 4), mgp = c(3, 0.5, 0),tcl = 0.1)
  barplot(as.matrix(t(readnumbers[, c(3, 2)])), ylab = "Number of reads", 
        legend.text = c("Kept", "Discarded"), args.legend = list(ncol = 2, 
            y = max(readnumbers$raw_reads + readnumbers$processed_reads) * 
                1.1, bty = "n"), col = c("royalblue", "aliceblue"), las = 2)
    dev.off()
    
    result_otu <- aggregate(result2[,-c(1)],by=list(OTUId=result$OTU),sum)
    
    write.table(result_otu, paste(getwd(), "/", folder.name, "TaxonomicTables/", folder.name, 
        "_otutable.txt", sep = ""), quote = F, row.names = F, sep = "\t")

    cat(paste("merge = ", merge, ",",
              "pool = ", pool, ",",
              "truncate = ", truncate, ",",
              "filter = ", filter, ",",
              "readlength = ", readlength, ",",
              "min.read.abundance = ", min.read.abundance),
        file = paste(folder.name,"log.txt",sep="_"))
    
    } 
