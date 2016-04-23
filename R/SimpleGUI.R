SimpleGUI <- function(forward.reads = NULL, forward.primer = NULL, name.separator = "_", 
    folder.name = "", usearch.path, refDB, gut.specific = F, 
    taxman = F, named.species = F) {
    
  refDB <- system.file("extdata/silva_full.fasta",package="mare")
  
      if (taxman) {
        
        system(paste("cp ", refDB, " refDB.fasta", sep = ""))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/;/,/g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/ //g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/subsp./subsp/g"), 
            "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/Bacteria/;tax=d:Bacteria/g"), 
            "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/Bacteria,/Bacteria,p:/g"), 
            "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/[(]//g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/[)]//g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/[|]//g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/[/]/_/g"), "refDB.fasta"))
        system(paste("sed -i", shQuote(""), "-E", shQuote("s/[-]//g"), "refDB.fasta"))
        
        database <- seqinr::read.fasta("refDB.fasta", as.string = T, forceDNAtolower = F)
        for (i in seq_along(names(database))) {
            names(database)[i] <- paste(strsplit(names(database)[i], split = ",")[[1]][1], 
                ",", strsplit(names(database)[i], split = ",")[[1]][2], ",c:", 
                strsplit(names(database)[i], split = ",")[[1]][3], ",o:", strsplit(names(database)[i], 
                  split = ",")[[1]][4], ",f:", strsplit(names(database)[i], 
                  split = ",")[[1]][5], ",g:", strsplit(names(database)[i], 
                  split = ",")[[1]][6], ",s:", strsplit(names(database)[i], 
                  split = ",")[[1]][7], sep = "")
        }
        
        database <- database[names(database)[grepl(":NA", names(database), fixed = T) == F]]
        database <- database[names(database)[grepl(":,", names(database), fixed = T) == F]]
        
        if (named.species) {
            database <- database[names(database)[grepl("uncultured", names(database)) == F]]
            database <- database[names(database)[grepl("unidentified", names(database)) == F]]
            database <- database[names(database)[grepl("s:bacterium", names(database)) == F]]
            database <- database[names(database)[grepl("metagenome", names(database)) == F]]
            database <- database[names(database)[grepl("sp.", names(database), fixed = T) == F]]
        }
        
        if (gut.specific) {
            databasetaxonomy <- data.frame(p = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][2]), c = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][3]), o = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][4]), f = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][5]), g = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][6]), s = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][7]))
            
            databaseselected <- database[databasetaxonomy$g == "g:Millisia" | 
                databasetaxonomy$g == "g:Rhodococcus" | databasetaxonomy$g == 
                "g:Actinomyces" | databasetaxonomy$g == "g:Bifidobacterium" | 
                databasetaxonomy$g == "g:Corynebacterium" | databasetaxonomy$g == 
                "g:Micrococcus" | databasetaxonomy$g == "g:Rothia" | databasetaxonomy$g == 
                "g:Mycobacterium" | databasetaxonomy$g == "g:Propionibacterium" | 
                databasetaxonomy$g == "g:Streptomyces" | databasetaxonomy$f == 
                "f:Coriobacteriaceae" | databasetaxonomy$g == "g:Alcaligenes" | 
                databasetaxonomy$g == "g:Bordetella" | databasetaxonomy$g == 
                "g:Parasutterella" | databasetaxonomy$g == "g:Sutterella" | 
                databasetaxonomy$g == "g:Burkholderia" | databasetaxonomy$g == 
                "g:Ralstonia" | databasetaxonomy$g == "g:Neisseria" | databasetaxonomy$g == 
                "g:Oxalobacter" | databasetaxonomy$g == "g:Aeromonas" | databasetaxonomy$g == 
                "g:Granulicatella" | databasetaxonomy$g == "g:Citrobacter" | 
                databasetaxonomy$g == "g:Enterobacter" | databasetaxonomy$g == 
                "g:Escherichia-Shigella" | databasetaxonomy$g == "g:Klebsiella" | 
                databasetaxonomy$g == "g:Leminorella" | databasetaxonomy$g == 
                "g:Morganella" | databasetaxonomy$g == "g:Proteus" | databasetaxonomy$g == 
                "g:Salmonella" | databasetaxonomy$g == "g:Serratia" | databasetaxonomy$g == 
                "g:Yersinia" | databasetaxonomy$g == "g:Pasteurella" | databasetaxonomy$g == 
                "g:Haemophilus" | databasetaxonomy$g == "g:Pseudomonas" | databasetaxonomy$g == 
                "g:Bilophila" | databasetaxonomy$g == "g:Desulfovibrio" | databasetaxonomy$g == 
                "g:Campylobacter" | databasetaxonomy$g == "g:Helicobacter" | 
                databasetaxonomy$g == "g:Fusobacterium" | databasetaxonomy$g == 
                "g:Bacteroides" | databasetaxonomy$f == "f:Porphyromonadaceae" | 
                databasetaxonomy$f == "f:Prevotellaceae" | databasetaxonomy$g == 
                "g:Alistipes" | databasetaxonomy$g == "g:Aerococcus" | databasetaxonomy$g == 
                "g:Bacillus" | databasetaxonomy$g == "g:Enterococcus" | databasetaxonomy$f == 
                "f:Lactobacillaceae" | databasetaxonomy$g == "g:Listeria" | 
                databasetaxonomy$g == "g:Akkermansia" | databasetaxonomy$g == 
                "g:Staphylococcus" | databasetaxonomy$f == "f:Streptococcaceae" | 
                databasetaxonomy$g == "g:Leuconostoc" | databasetaxonomy$g == 
                "g:Weissella" | databasetaxonomy$f == "f:Christensenellaceae" | 
                databasetaxonomy$g == "g:Clostridium" | databasetaxonomy$g == 
                "g:Anaerofustis" | databasetaxonomy$g == "g:Eubacterium" | 
                databasetaxonomy$g == "g:Anaerococcus" | databasetaxonomy$g == 
                "g:Finegoldia" | databasetaxonomy$g == "g:Gemella" | databasetaxonomy$g == 
                "g:IncertaeSedis" | 
                databasetaxonomy$g == "g:Parvimonas" | databasetaxonomy$g == 
                "g:Peptoniphilus" | databasetaxonomy$g == "g:Tissierella" | 
                databasetaxonomy$g == "g:Anaerosporobacter" | databasetaxonomy$g == 
                "g:Anaerostipes" | databasetaxonomy$g == "g:Blautia" | databasetaxonomy$g == 
                "g:Butyrivibrio" | databasetaxonomy$g == "g:Coprococcus" | 
                databasetaxonomy$g == "g:Dorea" | databasetaxonomy$g == 
                "g:Lachnospira" |
                databasetaxonomy$g == "g:Roseburia" | databasetaxonomy$g == 
                "g:Peptostreptococcus" | databasetaxonomy$g == "g:Sporacetigenium" | 
                databasetaxonomy$f == "f:Erysipelotrichaceae" | databasetaxonomy$f == 
                "f:Acidaminococcaceae" | databasetaxonomy$g == "g:Veillonella" | 
                databasetaxonomy$g == "g:Allisonella" | databasetaxonomy$g == 
                "g:Dialister" | databasetaxonomy$g == "g:Megamonas" | databasetaxonomy$g == 
                "g:Megasphaera" | databasetaxonomy$g == "g:Mitsuokella" | databasetaxonomy$g == 
                "g:Anaerotruncus" | databasetaxonomy$g == "g:Faecalibacterium" | 
                databasetaxonomy$g == "g:Flavonifractor" | databasetaxonomy$g == 
                "g:Oscillibacter" | databasetaxonomy$g == "g:Oscillospira" | 
                databasetaxonomy$g == "g:Ruminococcus" | databasetaxonomy$g == 
                "g:Pseudoflavonifractor" | databasetaxonomy$g == "g:Sporobacter" | 
                databasetaxonomy$g == "g:Subdoligranulum" | databasetaxonomy$g == 
                "g:Butyricicoccus" | databasetaxonomy$g == "g:Anaerofilum"]
            
            seqinr::write.fasta(sequences = databaseselected, file.out = "refDB.fasta", 
                names = paste(names(databaseselected), ";", sep = ""), as.string = T)
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/tax=d:Bacteria;//g"), 
                "refDB.fasta"))
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/:,/:NA,/g"), 
                "refDB.fasta"))

        } else {
            seqinr::write.fasta(sequences = database, file.out = "refDB.fasta", 
                names = paste(names(database), ";", sep = ""), as.string = T)
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/tax=d:Bacteria;//g"), 
                "refDB.fasta"))
        }
    } else {
        if (gut.specific) {
            database <- seqinr::read.fasta(refDB, as.string = T, forceDNAtolower = F)
            databasetaxonomy <- data.frame(p = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][2]), c = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][3]), o = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][4]), f = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][5]), g = sapply(names(database), function(x) strsplit(x, 
                split = ",")[[1]][6]))
            databaseselected <- database[databasetaxonomy$g == "g:Millisia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Rhodococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Actinomyces;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Bifidobacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Corynebacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Micrococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Rothia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Mycobacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Propionibacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Streptomyces;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Coriobacteriaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Alcaligenes;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Bordetella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Parasutterella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Sutterella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Burkholderia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Ralstonia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Neisseria;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Oxalobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Aeromonas;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Granulicatella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Citrobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Enterobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Escherichia-Shigella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Klebsiella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Leminorella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Morganella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Proteus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Salmonella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Serratia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Yersinia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Pasteurella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Haemophilus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Pseudomonas;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Bilophila;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Desulfovibrio;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Campylobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Helicobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Fusobacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Bacteroides;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Porphyromonadaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$f == "f:Prevotellaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Alistipes;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Aerococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Bacillus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Enterococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Lactobacillaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Listeria;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Akkermansia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Staphylococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Streptococcaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Leuconostoc;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Weissella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Christensenellaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Clostridium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerofustis;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Eubacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Finegoldia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Gemella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:IncertaeSedis;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Parvimonas;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Peptoniphilus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Tissierella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerosporobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerostipes;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Blautia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Butyrivibrio;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Coprococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Dorea;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Lachnospira;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Roseburia;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Peptostreptococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Sporacetigenium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$f == "f:Erysipelotrichaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$f == "f:Acidaminococcaceae" & 
                !is.na(databasetaxonomy$f) | databasetaxonomy$g == "g:Veillonella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Allisonella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Dialister;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Megamonas;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Megasphaera;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Mitsuokella;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerotruncus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Faecalibacterium;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Flavonifractor;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Oscillibacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Oscillospira;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Ruminococcus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Pseudoflavonifractor;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Sporobacter;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Subdoligranulum;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Butyricicoccus;" & 
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Anaerofilum;" & 
                !is.na(databasetaxonomy$g)]
            seqinr::write.fasta(sequences = na.omit(databaseselected), file.out = "refDB.fasta", 
                names = names(databaseselected), as.string = T)
        } else{
      system(paste("cp ", refDB, "refDB.fasta"))
    }
} 

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
    trimmed_fwd <- paste(wd, list.files(path = wd, pattern = "_fwd.fastq"), sep = "")
    for (i in seq_along(trimmed_fwd)) {
            system(paste(usearch.path, " -fastx_truncate ", trimmed_fwd[i], 
                " -trunclen ", "150", " -fastqout ", strsplit(trimmed_fwd[i], 
                  "_", fixed = T)[[1]][1], "_truncated.fastq", sep = ""))
        }
    reads <- paste(wd, list.files(path = wd, pattern = "_truncated.fastq"), sep = "")
    for (i in seq_along(reads)) {
        ShortRead::writeFasta(ShortRead::readFastq(reads[i]), file = paste(strsplit(reads[i], 
            split = "[.]")[[1]][1], "fasta", sep = "."))
    }
    reads <- paste(wd, list.files(path = wd, pattern = "*.fasta"), sep = "")
    
    system(paste("rm ", wd, "*.fastq", sep = ""))  
    
    rownames(readnumbers) <- readnumbers$sample
    
    fl.uparse.string <- ""
    for (i in seq_along(reads)) {
        fl.uparse.string <- paste(fl.uparse.string, reads[i], sep = " ")
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
    dnaSet <- dnaSet[as.vector(counts[rev(order(counts))] > (0.0001 * sum(counts)))]
    Biostrings::writeXStringSet(dnaSet, paste(wd, "dereplicated.fasta", sep = ""), 
        format = "fasta")
 system(paste(usearch.path, " -uchime_denovo ", wd,  "dereplicated.fasta", " -nonchimeras ", wd,
            "nonchimeric.fasta", sep = ""))
    
if (length(list.files(pattern = paste(folder.name, "TaxonomicTables", sep = ""))) == 0) {
        system(paste("mkdir ", folder.name, "TaxonomicTables", sep = ""))
    }
    wd <- getwd()
      
    derep <- Biostrings::readDNAStringSet(paste(wd, "/", folder.name,
            "ProcessedReads/nonchimeric.fasta", sep = ""))
    result <- data.frame(OTUId = names(derep))
    dict <- Biostrings::PDict(derep, max.mismatch = NA, tb.end = 2)
    for (i in seq_along(reads)) {
            seq <- Biostrings::readDNAStringSet(reads[i])
            result[, strsplit(reads[i], split = "\\_\\ProcessedReads/|\\_|\\ProcessedReads/")[[1]][2]] <- rowSums(Biostrings::vcountPDict(dict, 
                seq))
        }
       
    for(i in readnumbers$sample) readnumbers[readnumbers$sample==i,"processed_reads"] <- sum(result[,i])
    write.table(readnumbers, "readnumbers.txt", quote = F, row.names = F, sep = "\t")
    readnumbers$raw_reads <- readnumbers$raw_reads - readnumbers$processed_reads 
  
   quartz();par(xpd = T, mar = c(9, 7, 4, 4), mgp = c(5, 1, 0)) 
   barplot(as.matrix(t(readnumbers[, c(3, 2)])), ylab = "Number of reads", 
        legend.text = c("Kept", "Discarded"), args.legend = list(ncol = 2, 
            y = max(readnumbers$raw_reads + readnumbers$processed_reads) * 
                1.1, bty = "n"), col = c("royalblue", "aliceblue"), 
        las = 2)
       
  write.table(result, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
  tobeannotated <- paste(wd, "/", folder.name, "ProcessedReads/nonchimeric.fasta", sep = "")
  readfile <- result
    
  system(paste(usearch.path, " -utax ", tobeannotated, " -db refDB.fasta", " -utax_cutoff 0", 
        " -strand both -utaxout ", wd, "/", folder.name, 
        "TaxonomicTables/", folder.name, "taxonomy.txt", sep = ""))
    
  taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "taxonomy.txt", sep = ""), header = F)[, -c(3:4)]
        names(taxonomy) <- c("OTUId", "taxonomy")
        taxonomy$species <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$species[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "s:")[[1]][2], split = "[(]")[[1]][1]
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
            folder.name, "_species_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(genus_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_genus_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(family_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_family_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(order_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_order_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(class_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_class_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
        write.table(phylum_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_phylum_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
    } 
