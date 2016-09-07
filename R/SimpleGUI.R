SimpleGUI <- function(forward.primer = NULL, name.separator = "_", 
    folder.name = "", usearch.path) {
  
  all.reads <- list.files(pattern=".fastq")
  rev.reads <- list.files(pattern="R2")
  forward.reads <- setdiff(all.reads,rev.reads)
  refDB <- system.file("extdata/silva_full.fasta",package="mare")
  
        
    system(paste("mkdir ", folder.name, "ProcessedReads", sep = ""))
    wd <- paste(getwd(), "/", folder.name, "ProcessedReads/", sep = "")
    
    readnumbers <- data.frame(sample = sapply(as.character(forward.reads), 
        function(x) strsplit(x, split = name.separator)[[1]][1]), raw_reads = sapply(forward.reads, 
        function(x) length(Biostrings::readDNAStringSet(x, format = "fastq"))))
    
    if (length(forward.primer) != 0) {
        for (i in forward.reads) {
            ShortRead::writeFastq(object = Biostrings::trimLRPatterns(Lpattern = forward.primer, 
                subject = ShortRead::readFastq(i), max.Lmismatch = 3, with.Lindels = T, Lfixed = F), 
                file = paste(wd, sapply(i, function(x) strsplit(x, name.separator, 
                  fixed = F)[[1]][1]), name.separator ,"fwd.fastq", sep = ""), mode = "w", 
                full = F, compress = F, qualityType = "Auto")
        }
    } else {
        for (i in seq_along(forward.reads)) {
            system(paste("cp ", forward.reads[i], " ", wd, forward.reads[i], sep = ""))
        }
        if (length(list.files(path = wd, pattern = "*.gz$")) > 0) {
            system(paste("gunzip", " ", wd, "*fastq.gz", sep = ""))
        }

          for (i in grep(list.files(path=wd), pattern='R2', inv=T, value=T)) {
      system(paste("mv ", wd, i, " " , wd, sapply(i, function(x) strsplit(x, name.separator, fixed = F)[[1]][1]), 
                name.separator ,"fwd.fastq", sep = ""))
       }
    }
    trimmed_fwd <- paste(wd, list.files(path = wd, pattern = "fwd.fastq"),  sep = "")
    
        for (i in seq_along(trimmed_fwd)) {
            system(paste(usearch.path, " -fastx_truncate ", trimmed_fwd[i], 
                " -trunclen 150", " -fastqout ", strsplit(trimmed_fwd[i], 
                  name.separator, fixed = T)[[1]][1], name.separator ,"truncated.fastq", sep = ""))
        }
        reads <- paste(wd, list.files(path = wd, pattern = "truncated.fastq"), sep = "")
    
for (i in seq_along(reads)) {
        ShortRead::writeFasta(ShortRead::readFastq(reads[i]), file = paste(strsplit(reads[i], 
            split = "[.]")[[1]][1], "fasta", sep = "."))
    }
    reads <- paste(wd, list.files(path = wd, pattern = "*.fasta"), sep = "")
    
    system(paste("rm ", wd, "*.fastq", sep = ""))  
    
    for (i in seq_along(reads)) {
        system(paste(usearch.path, " -fastx_truncate ", reads[i], " -trunclen ", 
            "150", " -fastaout ", paste(strsplit(reads[i], split = ".fasta"), 
                "equallength.fasta", sep = name.separator), sep = ""))
    }
    
    truncatedreads <- paste(wd, list.files(path = wd, pattern = "equallength.fasta$"), sep = "")
    
    if (length(list.files(path = wd, "fwd.fasta")) != 0) {
        system(paste("rm ", wd, "*fwd.fasta", sep = ""))  
    }
   
    if (length(list.files(path = wd, "truncated.fasta")) != 0) {
        system(paste("rm ", wd, "*truncated.fasta", sep = ""))  
    }
    
    for (i in seq_along(truncatedreads)) {
        system(paste("sed '-es/^>\\(.*\\)/>\\1;barcodelabel=", strsplit(truncatedreads[i], 
            split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2], 
            ";/' <", truncatedreads[i], " > ", paste(strsplit(truncatedreads[i], 
                split = ".fasta"), "labelled.fasta", sep = name.separator), sep = ""))
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
    dnaSet <- dnaSet[as.vector(counts[rev(order(counts))] > (0.00001 * sum(counts)))]
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
        result[, strsplit(labelledreads[i], split =  paste("ProcessedReads/|",name.separator,sep=""))[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict, 
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
            result[, strsplit(labelledreads[i], split =  paste("ProcessedReads/|",name.separator,sep=""))[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict, 
                seq))
        }
      write.table(result, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t") 

    rownames(readnumbers) <- readnumbers$sample
    for(i in readnumbers$sample) readnumbers[readnumbers$sample==i,"processed_reads"] <- sum(result[,i])#new
    write.table(readnumbers, paste(folder.name, "readnumbers.txt",sep="_"), quote = F, row.names = F, sep = "\t")
    readnumbers$raw_reads <- readnumbers$raw_reads - readnumbers$processed_reads
    
    
    quartz();par(xpd = T, mar = c(9, 7, 4, 4), mgp = c(5, 1, 0))
    barplot(as.matrix(t(readnumbers[, c(3, 2)])), ylab = "Number of reads", 
        legend.text = c("Kept", "Discarded"), args.legend = list(ncol = 2, 
            y = max(readnumbers$raw_reads + readnumbers$processed_reads) * 
                1.1, bty = "n"), col = c("royalblue", "aliceblue"), las = 2)
    
    cat(paste("merge = F", ",",
              "truncate T= ", ",",
              "filter = F", ",",
              "readlength = 150", ",",
              "min.read.prevalence = 0.0001"),
        file = paste(folder.name,"log.txt",sep="_"))
    
        wd <- getwd()
    

        tobeannotated <- paste(wd, "/", folder.name, "ProcessedReads/nonchimeric.fasta", sep = "")
        readfile <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), check.names = F)
        OTU <- ""
    

    system(paste(usearch.path, " -utax ", tobeannotated, " -db ", refDB, " -utax_cutoff ", 
        "0", " -strand both -utaxout ", wd, "/", folder.name, 
        "TaxonomicTables/", folder.name, OTU, "taxonomy.txt", sep = ""))
        
        taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "taxonomy.txt", sep = ""), header = F)[, -c(3:4)]
        names(taxonomy) <- c("OTUId", "taxonomy")
        taxonomy$species <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$species[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "s:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$species <- gsub("_","",taxonomy$species)
        taxonomy$genus <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$genus[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "g:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$family <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$family[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "f:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$order <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$order[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "o:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$class <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$class[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "c:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$phylum <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$phylum[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "p:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$species <- paste(taxonomy$genus, taxonomy$species, sep = ".")
        
        taxonomy$species <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, taxonomy$species, sep = "_")
        taxonomy$genus <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, sep = "_")
        taxonomy$family <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, sep = "_")
        taxonomy$order <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            sep = "_")
        taxonomy$class <- paste(taxonomy$phylum, taxonomy$class, sep = "_")
        
        annotated_readtable <- merge(taxonomy, readfile, by = "OTUId", all = T)
        species_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$species), 
            sum)
        genus_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$genus), 
            sum)
        family_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$family), 
            sum)
        order_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$order), 
            sum)
        class_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$class), 
            sum)
        phylum_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$phylum), 
            sum)
        write.table(species_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_species_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(genus_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_genus_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(family_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_family_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(order_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_order_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(class_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_class_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(phylum_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_phylum_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        
} 

    
    

    
