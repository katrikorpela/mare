TaxonomicTable <- function(usearch.path, refDB, annotate.reads = T, 
    folder.name = "", confidence.cutoff = 0.5) {
    
    wd <- getwd()
    
    if (annotate.reads) {
        tobeannotated <- paste(wd, "/", folder.name, "ProcessedReads/nonchimeric.fasta", sep = "")
        readfile <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), check.names = F)
        OTU <- ""
    } else {
        tobeannotated <- paste(wd, "/", folder.name, "ProcessedReads/otu.fasta", sep = "")
        readfile <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_otutable.txt", sep = ""), check.names = F)
        OTU <- "OTU"
    }
    

    system(paste(usearch.path, " -utax ", tobeannotated, " -db ", refDB, " -utax_cutoff ", 
        confidence.cutoff, " -strand both -utaxout ", wd, "/", folder.name, 
        "TaxonomicTables/", folder.name, OTU, "taxonomy.txt", sep = ""))
    
    
    if (confidence.cutoff == 0) {
        
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
        
    } else {
        taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "taxonomy.txt", sep = ""), header = F)[, c(1, 
            3)]
        names(taxonomy) <- c("OTUId", "taxonomy")
        taxonomy$genus <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$genus[j] <- strsplit(as.character(taxonomy$taxonomy[j]), 
            split = ",g:")[[1]][2]
        taxonomy$family <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$family[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = ",f:")[[1]][2], split = ",")[[1]][1]
        taxonomy$order <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$order[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = ",o:")[[1]][2], split = "[,]")[[1]][1]
        taxonomy$class <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$class[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = ",c:")[[1]][2], split = "[,]")[[1]][1]
        taxonomy$phylum <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$phylum[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = ",p:")[[1]][2], split = "[,]")[[1]][1]
        taxonomy$genus <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, sep = "_")
        taxonomy$family <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, sep = "_")
        taxonomy$order <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            sep = "_")
        taxonomy$class <- paste(taxonomy$phylum, taxonomy$class, sep = "_")
        annotated_readtable <- merge(taxonomy, readfile, by = "OTUId", all = T)
        genus_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$genus), 
            sum)
        family_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$family), 
            sum)
        order_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$order), 
            sum)
        class_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$class), 
            sum)
        phylum_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$phylum), 
            sum)
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
} 
