FormatRefDB <- function(refDB, usearch.path, confidence.file, gut.specific = F, 
    taxman = F, named.species = F) {
    
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
                                           databasetaxonomy$g == "g:Succinivibrio" |
                                           databasetaxonomy$g == "g:Parabacteroides" |
                                           databasetaxonomy$g == "g:Paraprevotella" |
                                           databasetaxonomy$g == "g:Treponema" |
                                           databasetaxonomy$g == "g:Fibrobacter" |
                                           databasetaxonomy$g == "g:Methanobrevibacter" |
                                           databasetaxonomy$g == "g:Butyricoccus" |
                                           databasetaxonomy$g == "g:Phascolarctobacterium" |
                                           databasetaxonomy$g == "g:Butyricimonas" |
                                            databasetaxonomy$g == "g:Turicibacter" |
                                           databasetaxonomy$g == "g:Odoribacter" |
                databasetaxonomy$g == "g:Rhodococcus" | databasetaxonomy$g == 
                "g:Actinomyces" | databasetaxonomy$g == "g:Bifidobacterium" | 
                databasetaxonomy$g == "g:Corynebacterium" | databasetaxonomy$g == 
                "g:Micrococcus" | databasetaxonomy$g == "g:Rothia" | databasetaxonomy$g == 
                "g:Mycobacterium" | databasetaxonomy$g == "g:Propionibacterium" | 
                databasetaxonomy$g == "g:Streptomyces" | databasetaxonomy$f == 
                "f:Coriobacteriaceae" | databasetaxonomy$g == "g:Alcaligenes" |
                databasetaxonomy$g == "g:Gemmiger" | 
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
                databasetaxonomy$f == "f:Erysipelotrichaceae" | 
                   databasetaxonomy$f == "f:Lachnospiraceae" | 
                  databasetaxonomy$f == "f:Ruminococcaceae" | 
                  databasetaxonomy$f == "f:Veillonellaceae" | 
                  databasetaxonomy$f == "f:Acidaminococcaceae" | databasetaxonomy$g == "g:Veillonella" | 
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
            
            seqinr::write.fasta(sequences = databaseselected, file.out = "GutDB.fasta", 
                names = paste(names(databaseselected), ";", sep = ""), as.string = T)
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/tax=d:Bacteria;//g"), 
                "GutDB.fasta"))
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/:,/:NA,/g"), 
                "GutDB.fasta"))
            
            system(paste(usearch.path, " -makeudb_utax ", "GutDB.fasta", " -output ", 
                "GutDB.udb", " -taxconfsin ", confidence.file, sep = ""))
            
            
        } else {
            seqinr::write.fasta(sequences = database, file.out = "refDB.fasta", 
                names = paste(names(database), ";", sep = ""), as.string = T)
            system(paste("sed -i", shQuote(""), "-E", shQuote("s/tax=d:Bacteria;//g"), 
                "refDB.fasta"))
            
            system(paste(usearch.path, " -makeudb_utax ", "refDB.fasta", " -output ", 
                "refDB.udb", " -taxconfsin ", confidence.file, sep = ""))
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
            databaseselected <- database[databasetaxonomy$g == "g:Gemmiger;" & 
                                          !is.na(databasetaxonomy$g) |  databasetaxonomy$g == "g:Succinivibrio" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Parabacteroides" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Paraprevotella" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Treponema" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Fibrobacter" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Methanobrevibacter" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Butyricoccus" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Phascolarctobacterium" &
                                          !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Butyricimonas" &
                                          !is.na(databasetaxonomy$g) |  databasetaxonomy$g == "g:Turicibacter" &
                                         !is.na(databasetaxonomy$g) |  databasetaxonomy$g == "g:Odoribacter" &
                !is.na(databasetaxonomy$g) | databasetaxonomy$g == "g:Millisia;" & 
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
            seqinr::write.fasta(sequences = na.omit(databaseselected), file.out = "GutDB.fasta", 
                names = names(databaseselected), as.string = T)
            
            system(paste(usearch.path, " -makeudb_utax ", "GutDB.fasta", " -output ", 
                "GutDB.udb", " -taxconfsin ", confidence.file, sep = ""))
            
        } else {
            system(paste(usearch.path, " -makeudb_utax ", refDB, " -output ", 
                "refDB.udb", " -taxconfsin ", confidence.file, sep = ""))
        }
    }
} 
