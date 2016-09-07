CopyFiles <- function(folders = "*TaxonomicTables", files = "*table.txt", meta = NULL, folder.name = "",
    copy.to) {
    system(paste("mkdir", copy.to, sep = " "))
    system(paste("cp ", getwd(), "/", folders, "/", files, " ", copy.to, sep = ""))
    system(paste("cp ", getwd(), "/", folder.name , "_readnumbers.txt", " ", copy.to, sep = ""))
    if (length(meta) != 0) {
        system(paste("cp ", getwd(), "/", meta, " ", copy.to, sep = ""))
    }
    setwd(copy.to)
} 
