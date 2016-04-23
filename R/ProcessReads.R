ProcessReads <- function(forward.reads = list.files(pattern = "_R1_"), 
                         reverse.reads = list.files(pattern = "_R2_"), 
                         name.separator = "_", 
    usearch.path, forward.primer = NULL, reverse.primer = NULL, merge = F, 
    min.merged.length = NULL, min.overlap = NULL, max.mismatches = NULL, 
    truncate = F, trucation.length = NULL, 
    filter = F, trunc.quality = 2, max.expected.error = 1, min.readlength = NULL, 
    readlength = NULL, min.read.prevalence = 0.0001, folder.name = "") {
    
    system(paste("mkdir ", folder.name, "ProcessedReads", sep = ""))
    wd <- paste(getwd(), "/", folder.name, "ProcessedReads/", sep = "")
    
    readnumbers <- data.frame(sample = sapply(as.character(forward.reads), 
        function(x) strsplit(x, split = "_")[[1]][1]), raw_reads = sapply(forward.reads, 
        function(x) length(Biostrings::readDNAStringSet(x, format = "fastq"))))
    
    if (length(forward.primer) != 0) {
        for (i in forward.reads) {
            ShortRead::writeFastq(object = Biostrings::trimLRPatterns(Lpattern = forward.primer, 
                subject = ShortRead::readFastq(i), max.Lmismatch = 3, with.Lindels = T), 
                file = paste(wd, sapply(i, function(x) strsplit(x, name.separator, 
                  fixed = T)[[1]][1]), "_fwd.fastq", sep = ""), mode = "w", 
                full = F, compress = F, qualityType = "Auto")
        }
    } else {
        for (i in seq_along(forward.reads)) {
            system(paste("cp ", forward.reads[i], " ", wd, sapply(forward.reads[i], 
                function(x) strsplit(x, name.separator, fixed = T)[[1]][1]), 
                "_fwd.fastq.gz", sep = ""))
        }
        if (length(list.files(path = wd, pattern = "*.gz$")) > 0) {
            system(paste("gunzip", " ", wd, "*_fwd.fastq.gz", sep = ""))
        }
    }
    trimmed_fwd <- paste(wd, list.files(path = wd, pattern = "_fwd.fastq"), 
        sep = "")
    
    if (length(reverse.reads) != 0) {
        if (length(reverse.primer) != 0) {
            for (i in reverse.reads) {
                ShortRead::writeFastq(Biostrings::trimLRPatterns(Lpattern = reverse.primer, 
                  subject = ShortRead::readFastq(i), max.Lmismatch = 3, with.Lindels = T), 
                  file = paste(wd, sapply(i, function(x) strsplit(x, name.separator, 
                    fixed = T)[[1]][1]), "_rev.fastq", sep = ""), mode = "w", 
                  full = F, compress = F, qualityType = "Auto")
            }
        } else {
            for (i in seq_along(reverse.reads)) {
                system(paste("cp ", reverse.reads[i], " ", wd, sapply(reverse.reads[i], 
                  function(x) strsplit(x, name.separator, fixed = T)[[1]][1]), 
                  "_rev.fastq.gz", sep = ""))
            }
            system(paste("gunzip", " ", wd, "*_rev.fastq.gz", sep = ""))
        }
        trimmed_rev <- paste(wd, list.files(path = wd, pattern = "_rev.fastq"), 
            sep = "")
    }
    
    reads <- paste(wd, list.files(path = wd, pattern = "_fwd.fastq"), sep = "")
    
    if (merge) {
        for (i in seq_along(trimmed_fwd)) {
            system(paste(usearch.path, " -fastq_mergepairs ", trimmed_fwd[i], 
                " -reverse ", trimmed_rev[i], " -fastq_minmergelen ", min.merged.length, 
                " -fastq_minovlen ", min.overlap, " -fastq_maxdiffs ", max.mismatches, 
                " -fastqout ", strsplit(trimmed_fwd[i], "_", fixed = T)[[1]][1], 
                "_merged.fastq", sep = ""))
        }
        reads <- paste(wd, list.files(path = wd, pattern = "_merged.fastq"), 
            sep = "")

    } else if (truncate) {
        for (i in seq_along(trimmed_fwd)) {
            system(paste(usearch.path, " -fastx_truncate ", trimmed_fwd[i], 
                " -trunclen ", trucation.length, " -fastqout ", strsplit(trimmed_fwd[i], 
                  "_", fixed = T)[[1]][1], "_truncated.fastq", sep = ""))
        }
        reads <- paste(wd, list.files(path = wd, pattern = "_truncated.fastq"), 
            sep = "")
    
    }
    
    if (filter) {
        for (i in seq_along(reads)) {
            system(paste(usearch.path, " -fastq_filter ", reads[i], " -fastq_truncqual ", 
                trunc.quality, " -fastq_maxee ", max.expected.error, " -fastq_minlen ", 
                min.readlength, " -fastqout ", strsplit(reads[i], split = ".fastq")[[1]][1], 
                "_filtered.fastq", " -fastaout ", strsplit(reads[i], split = ".fastq")[[1]][1], 
                "_filtered.fasta", sep = ""))
        }
    } else for (i in seq_along(reads)) {
        ShortRead::writeFasta(ShortRead::readFastq(reads[i]), file = paste(strsplit(reads[i], 
            split = "[.]")[[1]][1], "fasta", sep = "."))
    }
    reads <- paste(wd, list.files(path = wd, pattern = "*.fasta"), sep = "")
    
    system(paste("rm ", wd, "*.fastq", sep = ""))  
    
    readnumbers$processed_reads <- sapply(list.files(path = wd, pattern = "*.fasta$"), 
        function(x) length(Biostrings::readDNAStringSet(paste(wd, x, sep = ""), 
            format = "fasta")))
    
    #totalreads <- sum(readnumbers$processed_reads)
    #if (length(min.read.prevalence) != 0) 
    #    minprev <- min.read.prevalence else minprev <- 0.0001
    
    for (i in seq_along(reads)) {
        system(paste(usearch.path, " -fastx_truncate ", reads[i], " -trunclen ", 
            readlength, " -fastaout ", paste(strsplit(reads[i], split = ".fasta"), 
                "equallength.fasta", sep = "_"), sep = ""))
    }
    
    truncatedreads <- paste(wd, list.files(path = wd, pattern = "equallength.fasta$"), sep = "")
    
    if (length(list.files(path = wd, "fwd.fasta")) != 0) {
        system(paste("rm ", wd, "*fwd.fasta", sep = ""))  
    }
    if (length(list.files(path = wd, "rev.fasta")) != 0) {
        system(paste("rm ", wd, "*rev.fasta", sep = "")) 
    }
    if (length(list.files(path = wd, "filtered.fasta")) != 0) {
        system(paste("rm ", wd, "*filtered.fasta", sep = "")) 
    }
    if (length(list.files(path = wd, "merged.fasta")) != 0) {
        system(paste("rm ", wd, "*merged.fasta", sep = ""))  
    }
    if (length(list.files(path = wd, "truncated.fasta")) != 0) {
        system(paste("rm ", wd, "*truncated.fasta", sep = ""))  
    }
    
    for (i in seq_along(truncatedreads)) {
        system(paste("sed '-es/^>\\(.*\\)/>\\1;barcodelabel=", strsplit(truncatedreads[i], 
            split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2], 
            ";/' <", truncatedreads[i], " > ", paste(strsplit(truncatedreads[i], 
                split = ".fasta"), "labelled.fasta", sep = "_"), sep = ""))
    }
    system(paste("rm ", wd, "*equallength.fasta", sep = ""))  
    labelledreads <- paste(wd, list.files(path = wd, pattern = "*labelled.fasta$"), sep = "")
    fl.uparse.string <- ""
    for (i in seq_along(labelledreads)) {
        fl.uparse.string <- paste(fl.uparse.string, labelledreads[i], sep = " ")
    }

    system(paste("cat ", fl.uparse.string, " > ", wd, "all.fasta", sep = ""))
    
    dnaSet = Biostrings::readDNAStringSet(paste(wd, "all.fasta", sep = ""), 
        format = "fasta")
    if (!is(dnaSet, "DNAStringSet")) {
        dnaSet <- Biostrings::DNAStringSet(dnaSet)
    }
    if (is.null(names(dnaSet))) {
        message("No names attribute found in dnaSet object...", "using artifically generated names")
        names(dnaSet) <- paste("read", 1:length(dnaSet), sep = "-")
    }
    dnaSet <- suppressWarnings(dnaSet[order(dnaSet)]) 
    counts <- BiocGenerics::table(dnaSet)
    dnaSet <- unique(dnaSet)
    names(dnaSet) <- paste0(paste("read", 1:length(names(dnaSet)), sep = "_"), 
        ";size=", as.integer(counts)) 
    dnaSet <- dnaSet[rev(order(counts))]
    dnaSet <- dnaSet[as.vector(counts[rev(order(counts))] > (min.read.prevalence * sum(counts)))]
    Biostrings::writeXStringSet(dnaSet, paste(wd, "dereplicated.fasta", sep = ""), 
        format = "fasta")
 
    system(paste(usearch.path, " -uchime_denovo ", wd,  "dereplicated.fasta", " -nonchimeras ", wd,
            "nonchimeric.fasta", sep = ""))
    system(paste(usearch.path, " -cluster_otus ", wd, "nonchimeric.fasta", 
        " -otus ", wd, "otu.fasta", " -otu_radius_pct 3", 
        " -relabel OTU_", " -sizein ", " -sizeout ", sep = ""))

    derep <- Biostrings::readDNAStringSet(paste(wd, "otu.fasta", sep = ""))
    result <- data.frame(OTUId = names(derep))
    dict <- Biostrings::PDict(derep, max.mismatch = NA, tb.end = 2)
    for (i in seq_along(labelledreads)) {
        seq <- Biostrings::readDNAStringSet(labelledreads[i])
        result[, strsplit(labelledreads[i], split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict, 
            seq))
    }
    
if (length(list.files(pattern = paste(folder.name, "TaxonomicTables", sep = ""))) == 0) {
        system(paste("mkdir ", folder.name, "TaxonomicTables", sep = ""))
    }
    wd <- getwd()

    write.table(result, paste(wd, "/", folder.name, "TaxonomicTables/", folder.name, 
        "_otutable.txt", sep = ""), quote = F, row.names = F, sep = "\t")

   derep <- Biostrings::readDNAStringSet(paste(wd, "/", folder.name, "ProcessedReads/nonchimeric.fasta", sep = ""))
        result <- data.frame(OTUId = names(derep))
        dict <- Biostrings::PDict(derep, max.mismatch = NA, tb.end = 2)
        for (i in seq_along(labelledreads)) {
            seq <- Biostrings::readDNAStringSet(labelledreads[i])
            result[, strsplit(labelledreads[i], split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict, 
                seq))
        }
      write.table(result, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t") 

    rownames(readnumbers) <- readnumbers$sample
    for(i in readnumbers$sample) readnumbers[readnumbers$sample==i,"processed_reads"] <- sum(result[,i])#new
    write.table(readnumbers, "readnumbers.txt", quote = F, row.names = F, sep = "\t")
    readnumbers$raw_reads <- readnumbers$raw_reads - readnumbers$processed_reads
    
    
    quartz();par(xpd = T, mar = c(9, 7, 4, 4), mgp = c(5, 1, 0))
    barplot(as.matrix(t(readnumbers[, c(3, 2)])), ylab = "Number of reads", 
        legend.text = c("Kept", "Discarded"), args.legend = list(ncol = 2, 
            y = max(readnumbers$raw_reads + readnumbers$processed_reads) * 
                1.1, bty = "n"), col = c("royalblue", "aliceblue"), las = 2)
    
    cat(paste("merge = ", merge, ",",
              "truncate = ", truncate, ",",
              "filter = ", filter, ",",
              "readlength = ", readlength, ",",
              "min.read.prevalence = ", min.read.prevalence),
        file = paste(folder.name,"log.txt",sep="_"))
    
    } 
