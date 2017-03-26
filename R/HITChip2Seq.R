HITChip2Seq <- function(HITChip_species_table){

phylogenypath <-  system.file("extdata/HITChip_phylogeny.txt",package="mare")
phylogeny <- data.frame(read.delim(phylogenypath))
phylogeny$taxonomy <- as.character(phylogeny$taxonomy)
L3 <- read.delim(HITChip_species_table, check.names = F)
names(L3)[1] <- "X"
L3$X <- sapply(L3$X, function(x) gsub(x, pattern = " ", replacement = "."))
L3taxonomy <- data.frame(X=L3$X)

for(i in L3taxonomy$X) L3taxonomy$taxonomy[L3taxonomy$X==i] <- phylogeny$taxonomy[phylogeny$Phylotype==i]

L3taxonomy$genus <- "NA"; for (j in 1:nrow(L3taxonomy)) L3taxonomy$genus[j] <- strsplit(as.character(L3taxonomy$taxonomy[j]), split = "_")[[1]][5]
L3taxonomy$family <- "NA"; for (j in 1:nrow(L3taxonomy)) L3taxonomy$family[j] <- strsplit(as.character(L3taxonomy$taxonomy[j]), split = "_")[[1]][4]
L3taxonomy$order <- "NA"; for (j in 1:nrow(L3taxonomy)) L3taxonomy$order[j] <- strsplit(as.character(L3taxonomy$taxonomy[j]), split = "_")[[1]][3]
L3taxonomy$class <- "NA"; for (j in 1:nrow(L3taxonomy)) L3taxonomy$class[j] <- strsplit(as.character(L3taxonomy$taxonomy[j]), split = "_")[[1]][2]
L3taxonomy$phylum <- "NA"; for (j in 1:nrow(L3taxonomy)) L3taxonomy$phylum[j] <- strsplit(as.character(L3taxonomy$taxonomy[j]), split = "_")[[1]][1]

L3taxonomy$genus <- paste(L3taxonomy$phylum, L3taxonomy$class, L3taxonomy$order,  L3taxonomy$family, L3taxonomy$genus, sep = "_")
L3taxonomy$family <- paste(L3taxonomy$phylum, L3taxonomy$class, L3taxonomy$order, L3taxonomy$family, sep = "_")
L3taxonomy$order <- paste(L3taxonomy$phylum, L3taxonomy$class, L3taxonomy$order, sep = "_")
L3taxonomy$class <- paste(L3taxonomy$phylum, L3taxonomy$class, sep = "_")

genus_table <- aggregate(L3[, -c(1)], by = list(taxon = L3taxonomy$genus),  sum)
family_table <- aggregate(L3[, -c(1)], by = list(taxon = L3taxonomy$family),  sum)
order_table <- aggregate(L3[, -c(1)], by = list(taxon = L3taxonomy$order),  sum)
class_table <- aggregate(L3[, -c(1)], by = list(taxon = L3taxonomy$class),  sum)
phylum_table <- aggregate(L3[, -c(1)], by = list(taxon = L3taxonomy$phylum),  sum)
rownames(phylum_table)<-phylum_table$taxon

species_table <- L3
names(species_table)[1]<-"taxon"
species_table$taxon <- paste(L3taxonomy$phylum, L3taxonomy$class, L3taxonomy$order,  L3taxonomy$family, L3taxonomy$genus, L3$X,sep = "_")

otutable <- L3
otutable[,-1][otutable[,-1]<mean(as.numeric(unlist(otutable[,-1])))]<-0
names(otutable)[1]<-"OTUId"

readnumbers <- data.frame(sample=colnames(L3)[-c(1)],raw_reads=colSums(L3[,-1]),processed_reads=colSums(L3[,-1]))

 write.table(genus_table, "HITChip_genus_table.txt", quote = F, row.names = F, 
            sep = "\t")
        write.table(family_table,  "HITChip_family_table.txt", quote = F, row.names = F, 
            sep = "\t")
        write.table(order_table, "HITChip_order_table.txt", quote = F, row.names = F, 
            sep = "\t")
        write.table(class_table,  "HITChip_class_table.txt",  quote = F, row.names = F, 
            sep = "\t")
        write.table(phylum_table, "HITChip_phylum_table.txt",  quote = F, row.names = F, 
            sep = "\t")
        write.table(species_table, "HITChip_species_table.txt",  quote = F, row.names = F, 
            sep = "\t")
          write.table(species_table, "HITChip_annotated_read_table.txt",  quote = F, row.names = F, 
            sep = "\t")
        write.table(otutable, "HITChip_otutable.txt",  quote = F, row.names = F, 
            sep = "\t")
         write.table(readnumbers, "HITChip_readnumbers.txt",  quote = F, row.names = F, 
            sep = "\t")
}
