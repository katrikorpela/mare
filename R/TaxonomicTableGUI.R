TaxonomicTableGUI <- function(usearch.path, refDB, annotate.reads = T, 
    folder.name = "", confidence.cutoff = 0) {
    
   refDB <- system.file(paste("extdata/",refDB,sep=""),package="mare")
 wd <- getwd()
    
        tobeannotated <- paste(wd, "/", folder.name, "ProcessedReads/finalreads.fasta", sep = "")
        readfile <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "_readtable.txt", sep = ""), check.names = F)
        OTU <- ""

    system(paste(usearch.path, " -utax ", tobeannotated, " -db ", refDB, " -utax_cutoff ", 
        confidence.cutoff, " -strand both -utaxout ", wd, "/", folder.name, 
        "TaxonomicTables/", folder.name, OTU, "taxonomy.txt", sep = ""))

    
    if (confidence.cutoff == 0) {    
      
        taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "taxonomy.txt", sep = ""), header = F)[, -c(3:4)]
        names(taxonomy) <- c("OTUId", "taxonomy")
        
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
       
        taxonomy$genus <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, sep = "_")
        taxonomy$family <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, sep = "_")
        taxonomy$order <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            sep = "_")
        taxonomy$class <- paste(taxonomy$phylum, taxonomy$class, sep = "_")
        
        annotated_readtable <- merge(taxonomy, readfile, by = "OTUId", all = T)
       
        genus_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$genus), sum)
        family_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$family), sum)
        order_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$order), sum)
        class_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$class), sum)
        phylum_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$phylum), sum)
        
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
        #readtable <- annotated_readtable[,-c(2:7)]
    
      taxonomy$species <- "NA"
        for (j in 1:nrow(taxonomy)) taxonomy$species[j] <- strsplit(strsplit(as.character(taxonomy$taxonomy[j]), 
            split = "s:")[[1]][2], split = "[(]")[[1]][1]
        taxonomy$species <- gsub("_","",taxonomy$species)
        taxonomy$species <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, taxonomy$species, sep = "_")
         annotated_readtable <- merge(taxonomy, readfile, by = "OTUId", all = T)
         species_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$species), sum)
         write.table(species_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_species_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
         readtable <- annotated_readtable[,-c(2:8)]
         
    } else {
         if(grepl(pattern=".udb",x=refDB)){   
        taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "taxonomy.txt", sep = ""), header = F)[, -c(2,4)]
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
       
        genus_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$genus), sum)
        family_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$family), sum)
        order_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$order), sum)
        class_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$class), sum)
        phylum_table <- aggregate(annotated_readtable[, -c(1:7)], by = list(taxon = annotated_readtable$phylum), sum)
        
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
        readtable <- annotated_readtable[,-c(2:7)]
  
         } else print("RefDB must be in .udb format if confidence.cutoff is not 0!")
    }
        
    
    readtable$OTUId<-as.character(readtable$OTUId)
    names(readtable)[1]<-"taxon"
    for(i in na.omit(unique(annotated_readtable$genus))) readtable$taxon[annotated_readtable$genus==i&!is.na(annotated_readtable$genus)] <- paste(i,"read",1:nrow(annotated_readtable[annotated_readtable$genus==i&!is.na(annotated_readtable$genus),]),sep="_")
    readtable <- readtable[order(readtable$taxon),] 
    
     write.table(readtable, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, OTU, "_annotated_read_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
     
library(metacoder)     
     
gen <- genus_table
rownames(gen) <- gen$taxon
gen <- gen[,-1]

seqs <- list()
for(i in rownames(gen)){
  for(h in colnames(gen)){
     if(gen[i,h]>0) seqs[[paste(h,i)]]<- round(gen[i,h]/1000)
  }
}
newseqs <- rep(names(seqs), seqs)

tmp <-  extract_taxonomy(input=newseqs,regex = "^(.*)\\ (.*)",
                         key=c(id = "obs_info","class"),class_sep = "_")

heat_tree(tmp, node_size = n_obs*1000,
                    node_label = ifelse(name == "NA", NA, name), 
                    node_color = n_obs*1000,
                    node_color_axis_label = "Number of reads per taxon",
                    node_color_range=c("gray90","skyblue","pink","hotpink","red","firebrick"),
                    initial_layout = "reingold-tilford", layout = "davidson-harel",
                    node_size_range=c(0.01,0.04),
                     overlap_avoidance = 0.5,
                    node_label_max=200,
                    node_label_size_range=c(0.01,0.015),
                    title="Abundance of taxa in the dataset",title_size=0.03,
          output_file=paste(folder.name,"HeatTree.pdf",sep="_"))

detach(package:metacoder)
} 
