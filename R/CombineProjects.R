CombineProjects <- function(project1folder,project2folder, meta1, meta2, batch.effect = F, all.taxa = T){

project1meta <- read.delim(paste(project1folder,"/",meta1,sep=""))
project1meta$project <- "project1"
project2meta <- read.delim(paste(project2folder,"/",meta2,sep=""))
project2meta$project <- "project2"
for(i in setdiff(names(project1meta),names(project2meta))) project2meta[,i] <- NA
for(i in setdiff(names(project2meta),names(project1meta))) project1meta[,i] <- NA
project1meta <- as.data.frame(apply(project1meta,FUN = as.character,MARGIN = 2))
project2meta <- as.data.frame(apply(project2meta,FUN = as.character,MARGIN = 2))
project1meta$ReadCount <- as.numeric(as.character(project1meta$ReadCount))
project2meta$ReadCount <-as.numeric(as.character(project2meta$ReadCount))
combinedmeta <- rbind(project1meta,project2meta)
write.table(combinedmeta, "combined_meta.txt", quote = F, row.names = F, sep = "\t")
  
if(length(list.files(path=project1folder,pattern="species_table", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="species_table", full.names=T))])>0) { 
  if(length(list.files(path=project2folder,pattern="species_table", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="species_table", full.names=T))])>0) { 
project1species <- read.delim(list.files(path=project1folder,pattern="_species_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="species_table", full.names=T))])
project2species <- read.delim(list.files(path=project2folder,pattern="_species_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="species_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1species),names(project2species))) project2species[,i] <- 0
for(i in setdiff(names(project2species),names(project1species))) project1species[,i] <- 0
} else {
  project1species <- project1species[,intersect(names(project1species),names(project2species))]
  project2species <- project2species[,intersect(names(project1species),names(project2species))]
}
if(ncol(project1species)>0){
  if(ncol(project2species)>0){
if(batch.effect){
speciesdiff <- list()
for(i in colnames(project2species)) speciesdiff[[i]] <- mean((project2species[,i]+1)/project2meta$ReadCount) / mean((project1species[,i]+1)/project1meta$ReadCount)
speciesStandard <- project2species
for(i in colnames(project2species)) {speciesStandard[,i]<- (project2species[,i]+1)/speciesdiff[[i]]}
for(i in rownames(speciesStandard)) speciesStandard[i,] <- speciesStandard[i,]/(rowSums(speciesStandard)[i]/rowSums(project2species)[i])
project2species <- round(speciesStandard)
}
 combinedspecies <- rbind(project1species,project2species)
write.table(combinedspecies, "combined_species_table.txt", quote = F, row.names = T, sep = "\t")
 }
}
}}

project1genus <- read.delim(list.files(path=project1folder,pattern="_genus_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="genus_table", full.names=T))])
project2genus <- read.delim(list.files(path=project2folder,pattern="_genus_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="genus_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1genus),names(project2genus))) project2genus[,i] <- 0
for(i in setdiff(names(project2genus),names(project1genus))) project1genus[,i] <- 0
} else {
  project1genus <- project1genus[,intersect(names(project1genus),names(project2genus))]
  project2genus <- project2genus[,intersect(names(project1genus),names(project2genus))]
}
if(batch.effect){
genusdiff <- list()
for(i in colnames(project2genus)) genusdiff[[i]] <- mean((project2genus[,i]+1)/project2meta$ReadCount) / mean((project1genus[,i]+1)/project1meta$ReadCount)
genusStandard <- project2genus
for(i in colnames(project2genus)) {genusStandard[,i]<- (project2genus[,i]+1)/genusdiff[[i]]}
for(i in rownames(genusStandard)) genusStandard[i,] <- genusStandard[i,]/(rowSums(genusStandard[i,])/rowSums(project2genus[i,]))
project2genus <- round(genusStandard)

}
combinedgenus <- rbind(project1genus,project2genus)
write.table(combinedgenus, "combined_genus_table.txt", quote = F, row.names = T, sep = "\t")

project1family <- read.delim(list.files(path=project1folder,pattern="_family_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="family_table", full.names=T))])
project2family <- read.delim(list.files(path=project2folder,pattern="_family_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="family_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1family),names(project2family))) project2family[,i] <- 0
for(i in setdiff(names(project2family),names(project1family))) project1family[,i] <- 0
} else {
  project1family <- project1family[,intersect(names(project1family),names(project2family))]
  project2family <- project2family[,intersect(names(project1family),names(project2family))]
}
if(batch.effect){
familydiff <- list()
for(i in colnames(project2family)) familydiff[[i]] <- mean((project2family[,i]+1)/project2meta$ReadCount) / mean((project1family[,i]+1)/project1meta$ReadCount)
familyStandard <- project2family
for(i in colnames(project2family)) {familyStandard[,i]<- (project2family[,i]+1)/familydiff[[i]]}
for(i in rownames(familyStandard)) familyStandard[i,] <- familyStandard[i,]/(rowSums(familyStandard[i,])/rowSums(project2family[i,]))
project2family <- round(familyStandard)
}
combinedfamily <- rbind(project1family,project2family)
write.table(combinedfamily, "combined_family_table.txt", quote = F, row.names = T, sep = "\t")

project1class <- read.delim(list.files(path=project1folder,pattern="_class_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="class_table", full.names=T))])
project2class <- read.delim(list.files(path=project2folder,pattern="_class_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="class_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1class),names(project2class))) project2class[,i] <- 0
for(i in setdiff(names(project2class),names(project1class))) project1class[,i] <- 0
} else {
  project1class <- project1class[,intersect(names(project1class),names(project2class))]
  project2class <- project2class[,intersect(names(project1class),names(project2class))]
}
if(batch.effect){
classdiff <- list()
for(i in colnames(project2class)) classdiff[[i]] <- mean((project2class[,i]+1)/project2meta$ReadCount) / mean((project1class[,i]+1)/project1meta$ReadCount)
classStandard <- project2class
for(i in colnames(project2class)) {classStandard[,i]<- (project2class[,i]+1)/classdiff[[i]]}
for(i in rownames(classStandard)) classStandard[i,] <- classStandard[i,]/(rowSums(classStandard[i,])/rowSums(project2class[i,]))
project2class <- round(classStandard)
}
combinedclass <- rbind(project1class,project2class)
write.table(combinedclass, "combined_class_table.txt", quote = F, row.names = T, sep = "\t")

project1order <- read.delim(list.files(path=project1folder,pattern="_order_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="order_table", full.names=T))])
project2order <- read.delim(list.files(path=project2folder,pattern="_order_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="order_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1order),names(project2order))) project2order[,i] <- 0
for(i in setdiff(names(project2order),names(project1order))) project1order[,i] <- 0
} else {
  project1order <- project1order[,intersect(names(project1order),names(project2order))]
  project2order <- project2order[,intersect(names(project1order),names(project2order))]
}
if(batch.effect){
orderdiff <- list()
for(i in colnames(project2order)) orderdiff[[i]] <- mean((project2order[,i]+1)/project2meta$ReadCount) / mean((project1order[,i]+1)/project1meta$ReadCount)
orderStandard <- project2order
for(i in colnames(project2order)) {orderStandard[,i]<- (project2order[,i]+1)/orderdiff[[i]]}
for(i in rownames(orderStandard)) orderStandard[i,] <- orderStandard[i,]/(rowSums(orderStandard[i,])/rowSums(project2order[i,]))
project2order <- round(orderStandard)
}
combinedorder <- rbind(project1order,project2order)
write.table(combinedorder, "combined_order_table.txt", quote = F, row.names = T, sep = "\t")

project1phylum <- read.delim(list.files(path=project1folder,pattern="_phylum_table.txt", full.names=T)[grepl("organised",x=list.files(path=project1folder,pattern="phylum_table", full.names=T))])
project2phylum <- read.delim(list.files(path=project2folder,pattern="_phylum_table.txt", full.names=T)[grepl("organised",x=list.files(path=project2folder,pattern="phylum_table", full.names=T))])
if(all.taxa){
for(i in setdiff(names(project1phylum),names(project2phylum))) project2phylum[,i] <- 0
for(i in setdiff(names(project2phylum),names(project1phylum))) project1phylum[,i] <- 0
} else {
  project1phylum <- project1phylum[,intersect(names(project1phylum),names(project2phylum))]
  project2phylum <- project2phylum[,intersect(names(project1phylum),names(project2phylum))]
}
if(batch.effect){
phylumdiff <- list()
for(i in colnames(project2phylum)) phylumdiff[[i]] <- mean((project2phylum[,i]+1)/project2meta$ReadCount) / mean((project1phylum[,i]+1)/project1meta$ReadCount)
phylumStandard <- project2phylum
for(i in colnames(project2phylum)) {phylumStandard[,i]<-(project2phylum[,i]+1)/phylumdiff[[i]]}
for(i in rownames(phylumStandard)) phylumStandard[i,] <- phylumStandard[i,]/(rowSums(phylumStandard[i,])/rowSums(project2phylum[i,]))
project2phylum <- round(phylumStandard)
}
combinedphylum <- rbind(project1phylum,project2phylum)
write.table(combinedphylum, "combined_phylum_table.txt", quote = F, row.names = T, sep = "\t")

}
