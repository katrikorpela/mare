CAZy <- function(species=NULL){

system(paste("grep bacteria ", species, " > species.txt"))
system("awk '{OFS=RS;$1=$1}1' species.txt > species2.txt")
system("rm species.txt")
system("mv species2.txt species.txt")
system("sed -n -e 's/^.*_//p' species.txt > species2.txt")
system("rm species.txt")
system("mv species2.txt species.txt")

sp <- read.delim(species)
sp <- data.frame(t(sp))
sp$species <- sapply(rownames(sp),function(x) strsplit(x, split="[.]")[[1]][2])
sp <- aggregate(sp[,-length(names(sp))],by=list(species=sp$species),sum)
rownames(sp)<- sp$species
sp <- sp[,-1]
sp <- sp[rownames(sp)!="",]

dbpath <- system.file("extdata/CE.txt",package="mare")
ce <- data.frame(read.delim(dbpath,header=F,sep="|")[,c(2,5)])
ce$V5 <- gsub(pattern="[.]", replacement="",x=ce$V5)
cetable <- array(dim=c(length(rownames(sp)),length(unique(ce$V2))))
rownames(cetable)<-rownames(sp)
colnames(cetable)<-unique(ce$V2)
for(i in rownames(cetable)){
  for(j in unique(ce$V2[grepl(pattern = i, x = ce$V5)])){
    cetable[i,j]<-1
  }
}
cetable[is.na(cetable)]<-0
cetable<-cetable[,colSums(cetable)>0]
cesample <- array(dim=c(length(names(sp)),length(colnames(cetable))))
rownames(cesample)<-names(sp)
colnames(cesample)<-colnames(cetable)
tmp <- list()
for(i in rownames(cesample)){
    for(k in colnames(cesample)){
      tmp[[paste(i,k)]] <- data.frame(sp=cetable[rownames(sp)[sp[,i]>0],k],
                        count=sp[sp[,i]>0,i])
      tmp[[paste(i,k)]]$summed <- tmp[[paste(i,k)]]$sp*tmp[[paste(i,k)]]$count
      cesample[i,k] <- sum(tmp[[paste(i,k)]]$summed)
     }
}
rownames(cesample)<-gsub(rownames(cesample),pattern = "[.]", replacement = "-")
write.table(cesample,"CAZy_table_CarbohydrateEsterases.txt",quote=F,sep="\t",row.names=T)

dbpath <- system.file("extdata/AA.txt",package="mare")
aa <- data.frame(read.delim(dbpath,header=F,sep="|")[,c(2,5)])
aa$V5 <- gsub(pattern="[.]", replacement="",x=aa$V5)
aatable <- array(dim=c(length(rownames(sp)),length(unique(aa$V2))))
rownames(aatable)<-rownames(sp)
colnames(aatable)<-unique(aa$V2)
for(i in rownames(aatable)){
  for(j in unique(aa$V2[grepl(pattern = i, x = aa$V5)])){
    aatable[i,j]<-1
  }
}
aatable[is.na(aatable)]<-0
aatable<-aatable[,colSums(aatable)>0]
aasample <- array(dim=c(length(names(sp)),length(colnames(aatable))))
rownames(aasample)<-names(sp)
colnames(aasample)<-colnames(aatable)
tmp <- list()
for(i in rownames(aasample)){
    for(k in colnames(aasample)){
      tmp[[paste(i,k)]] <- data.frame(sp=aatable[rownames(sp)[sp[,i]>0],k],
                        count=sp[sp[,i]>0,i])
      tmp[[paste(i,k)]]$summed <- tmp[[paste(i,k)]]$sp*tmp[[paste(i,k)]]$count
      aasample[i,k] <- sum(tmp[[paste(i,k)]]$summed)
     }
}
rownames(aasample)<-gsub(rownames(aasample),pattern = "[.]", replacement = "-")
write.table(aasample,"CAZy_table_AuxiliaryActivities.txt",quote=F,sep="\t",row.names=T)

dbpath <- system.file("extdata/GH.txt",package="mare")
gh <- data.frame(read.delim(dbpath,header=F,sep="|")[,c(2,5)])
gh$V5 <- gsub(pattern="[.]", replacement="",x=gh$V5)
ghtable <- array(dim=c(length(rownames(sp)),length(unique(gh$V2))))
rownames(ghtable)<-rownames(sp)
colnames(ghtable)<-unique(gh$V2)
for(i in rownames(ghtable)){
  for(j in unique(gh$V2[grepl(pattern = i, x = gh$V5)])){
    ghtable[i,j]<-1
  }
}
ghtable[is.na(ghtable)]<-0
ghtable<-ghtable[,colSums(ghtable)>0]
ghsample <- array(dim=c(length(names(sp)),length(colnames(ghtable))))
rownames(ghsample)<-names(sp)
colnames(ghsample)<-colnames(ghtable)
tmp <- list()
for(i in rownames(ghsample)){
    for(k in colnames(ghsample)){
      tmp[[paste(i,k)]] <- data.frame(sp=ghtable[rownames(sp)[sp[,i]>0],k],
                        count=sp[sp[,i]>0,i])
      tmp[[paste(i,k)]]$summed <- tmp[[paste(i,k)]]$sp*tmp[[paste(i,k)]]$count
      ghsample[i,k] <- sum(tmp[[paste(i,k)]]$summed)
     }
}
rownames(ghsample)<-gsub(rownames(ghsample),pattern = "[.]", replacement = "-")
write.table(ghsample,"CAZy_table_GlycosideHydrolases.txt",quote=F,sep="\t",row.names=T)

dbpath <- system.file("extdata/GT.txt",package="mare")
gt <- data.frame(read.delim(dbpath,header=F,sep="|")[,c(2,5)])
gt$V5 <- gsub(pattern="[.]", replacement="",x=gt$V5)
gttable <- array(dim=c(length(rownames(sp)),length(unique(gt$V2))))
rownames(gttable)<-rownames(sp)
colnames(gttable)<-unique(gt$V2)
for(i in rownames(gttable)){
  for(j in unique(gt$V2[grepl(pattern = i, x = gt$V5)])){
    gttable[i,j]<-1
  }
}
gttable[is.na(gttable)]<-0
gttable<-gttable[,colSums(gttable)>0]
gtsample <- array(dim=c(length(names(sp)),length(colnames(gttable))))
rownames(gtsample)<-names(sp)
colnames(gtsample)<-colnames(gttable)
tmp <- list()
for(i in rownames(gtsample)){
    for(k in colnames(gtsample)){
      tmp[[paste(i,k)]] <- data.frame(sp=gttable[rownames(sp)[sp[,i]>0],k],
                        count=sp[sp[,i]>0,i])
      tmp[[paste(i,k)]]$summed <- tmp[[paste(i,k)]]$sp*tmp[[paste(i,k)]]$count
      gtsample[i,k] <- sum(tmp[[paste(i,k)]]$summed)
     }
}
rownames(gtsample)<-gsub(rownames(gtsample),pattern = "[.]", replacement = "-")
write.table(gtsample,"CAZy_table_GlycosylTransferases.txt",quote=F,sep="\t",row.names=T)

#polysaccharide lyases
dbpath <- system.file("extdata/PL.txt",package="mare")
pl <- data.frame(read.delim(dbpath,header=F,sep="|")[,c(2,5)])
pl$V5 <- gsub(pattern="[.]", replacement="",x=pl$V5)
pltable <- array(dim=c(length(rownames(sp)),length(unique(pl$V2))))
rownames(pltable)<-rownames(sp)
colnames(pltable)<-unique(pl$V2)
for(i in rownames(pltable)){
  for(j in unique(pl$V2[grepl(pattern = i, x = pl$V5)])){
    pltable[i,j]<-1
  }
}
pltable[is.na(pltable)]<-0
pltable<-pltable[,colSums(pltable)>0]
plsample <- array(dim=c(length(names(sp)),length(colnames(pltable))))
rownames(plsample)<-names(sp)
colnames(plsample)<-colnames(pltable)
tmp <- list()
for(i in rownames(plsample)){
    for(k in colnames(plsample)){
      tmp[[paste(i,k)]] <- data.frame(sp=pltable[rownames(sp)[sp[,i]>0],k],
                        count=sp[sp[,i]>0,i])
      tmp[[paste(i,k)]]$summed <- tmp[[paste(i,k)]]$sp*tmp[[paste(i,k)]]$count
      plsample[i,k] <- sum(tmp[[paste(i,k)]]$summed)
     }
}
rownames(plsample)<-gsub(rownames(plsample),pattern = "[.]", replacement = "-")
write.table(plsample,"CAZy_table_PolysaccharideLyases.txt",quote=F,sep="\t",row.names=T)

#laita yhteen datat?

} 
