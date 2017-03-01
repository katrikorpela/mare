Blast <- function(finalreads, read.table, NCBItaxonomy, folder.name){

wd <- getwd()
    
system(paste("blastn -db rRNA_typestrains/prokaryotic_16S_ribosomal_RNA -query ",finalreads, 
" -out ", wd, "/", folder.name, 
"TaxonomicTables/blasttaxonomy.txt -culling_limit 1 -outfmt ", shQuote("6 qseqid staxid pident evalue"), " -remote", 
sep=""))

taxonomy <- read.delim(paste(wd, "/", folder.name, "TaxonomicTables/blasttaxonomy.txt", sep = ""), header = F)[, c(1:2)]
names(taxonomy) <- c("OTUId", "id")
taxonomy[,c("subspecies","species","genus","family","order","class","phylum")]<-NA

library(CHNOSZ)
nod <- getnodes(taxdir=NCBItaxonomy) 
nam <- getnames(taxdir=NCBItaxonomy) 

ranks <- list()
tax <- list()
for(i in taxonomy$id){ 
  ranks[[paste("rank",i,sep="")]] <- getrank(allparents(id=i, taxdir=NCBItaxonomy, nodes=nod), taxdir=NCBItaxonomy, nodes=nod)
  tax[[paste("tax",i,sep="")]] <- sciname(allparents(id=i, taxdir=NCBItaxonomy, nodes=nod), taxdir=NCBItaxonomy,names=nam)
  }
for(i in taxonomy$id){
  for(j in intersect(ranks[[paste("rank",i,sep="")]],names(taxonomy)[3:9])){
taxonomy[taxonomy$id==i,j] <- tax[[paste("tax",i,sep="")]][ranks[[paste("rank",i,sep="")]]==j]
}}
taxonomy$subspecies <- sapply(taxonomy$subspecies, function(x) gsub(" ", "", x, fixed = TRUE))
taxonomy$species[!is.na(taxonomy$subspecies)] <- taxonomy$subspecies[!is.na(taxonomy$subspecies)]
taxonomy$species <- sapply(taxonomy$species, function(x) gsub(" ", ".", x, fixed = TRUE))
taxonomy$species <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, taxonomy$species, sep = "_")
taxonomy$genus <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, taxonomy$genus, sep = "_")
taxonomy$family <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            taxonomy$family, sep = "_")
taxonomy$order <- paste(taxonomy$phylum, taxonomy$class, taxonomy$order, 
            sep = "_")
taxonomy$class <- paste(taxonomy$phylum, taxonomy$class, sep = "_")
   
readfile <- read.delim(read.table, check.names = F)     
annotated_readtable <- merge(taxonomy[,-2], readfile, by = "OTUId", all = T)
species_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$species), sum)
genus_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$genus),  sum)
family_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$family), sum)
order_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$order), sum)
class_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$class), sum)
phylum_table <- aggregate(annotated_readtable[, -c(1:8)], by = list(taxon = annotated_readtable$phylum), sum)

write.table(species_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_species_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
write.table(genus_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_genus_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
write.table(family_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_family_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
write.table(order_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_order_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
write.table(class_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_class_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
write.table(phylum_table, paste(wd, "/", folder.name, "TaxonomicTables/", 
            folder.name, "blast_phylum_table.txt", sep = ""), quote = F, row.names = F, 
            sep = "\t")
}