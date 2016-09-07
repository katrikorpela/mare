mareGUI <- function(){
library(fgui)
  
helpfunction <- function(function.name){
h <- help(function.name)
print(h)
}

guidefunction <- function(guide.name){
browseURL(system.file(paste(guide.name,".pdf",sep=""), package = "mare"))
}

wdfunction <- function(dir){
  setwd(dir)
}

fguiWindow(title="mare", text="Welcome to mare! 
Start by reading the help files and guides.
To analyse your data, first do the pre-processing
(functions ProcessReads and TaxonomicTables or Simple),
copy the files to a new folder for analysis and organise the data
(functions CopyFiles and Organise),
then explore the results by plotting
(functions PCoA, GroupPlot, CovariatePlot, Clusters, CorrelationMap), 
and finally conduct sigficance testing
(functions GroupTest, CovariateTest, ChangeTest),
and carbohydrate utilization prediction
(function CAZy).")
fguiNewMenu(c("Help", "SEPARATOR"))
mgui(helpfunction, title = c("Help","Help files"), modal=F, output=NULL,
    argText = list(function.name="Write the name of the function in quotations"))
mgui(guidefunction, title = c("Help","Guides"), 
     argOption=list(guide.name=c("mareGuide","mareBackground")), modal=F, output=NULL)
fguiNewMenu(c("Preprocess","SEPARATOR"))
mgui(wdfunction, title = c("Preprocess", "Set working directory"), output=NULL, 
     argText = list(dir = "Write the full path to the working directory"))
mgui(SimpleGUI, title = c("Preprocess", "Simple short cut to taxonomic tables"), output=NULL, 
     argFilename=list(usearch.path=NULL))
mgui(ProcessReads, title = c("Preprocess", "ProcessReads"), output=NULL, 
     argFilename=list(usearch.path=NULL))
mgui(TaxonomicTableGUI, title = c("Preprocess", "TaxonomicTable"), output=NULL,
     argFilename=list(usearch.path=NULL),#,refDB=NULL)
     argList=list(refDB=c("RDP.fasta","RDP_Gut.fasta","silva_full.fasta",            
 "silva_full_Gut.fasta", "silva_full_Gut_namedSP.fasta", "silva_full_namedSP.fasta",   
 "silva_v3v4.fasta","silva_v3v4_Gut.fasta","silva_v3v4_Gut_namedSP.fasta","silva_v3v4_namedSP.fasta",
 "RDP.udb","RDP_Gut.udb","silva_full.udb","silva_full_Gut_namedSP.udb",
"silva_full_namedSP.udb","silva_v3v4_Gut.udb","silva_v3v4_Gut_namedSP.udb","silva_v3v4_namedSP.udb")))
mgui(CopyFiles, title = c("Preprocess", "CopyFiles"), output=NULL)
mgui(Organise, title = c("Preprocess", "Organise"), output=NULL)
fguiNewMenu(c("Plot","SEPARATOR"))
mgui(PCoA, title = c("Plot", "PCoA"),argFilename=list(taxonomic.table=NULL, meta=NULL,output="m"))
mgui(Clusters, title = c("Plot", "Clusters"),argFilename=list(taxonomic.table=NULL, meta=NULL,output="m"))
mgui(CorrelationMap, title = c("Plot", "CorrelationMap"),argFilename=list(taxonomic.table=NULL, meta=NULL,output="m"))
mgui(GroupPlot, title = c("Plot", "GroupPlot"),argFilename=list(taxonomic.table=NULL, meta=NULL), output="m")
mgui(CovariatePlot, title = c("Plot", "CovariatePlot"),
     argFilename=list(taxonomic.table=NULL, meta=NULL), output="m")
fguiNewMenu(c("Test","SEPARATOR"))
mgui(GroupTest, title = c("Test", "GroupTest"),argFilename=list(taxonomic.table=NULL, meta=NULL), output="m")
mgui(CovariateTest, title = c("Test", "CovariateTest"),argFilename=list(taxonomic.table=NULL, meta=NULL), output="m")
mgui(ChangeTest, title = c("Test", "ChangeTest"),argFilename=list(taxonomic.table=NULL, meta=NULL), output="m")
fguiNewMenu(c("CAZy prediction","SEPARATOR"))
mgui(CAZy, title = c("CAZy prediction","CAZy"), argFilename=list(species=NULL),output=NULL)
}
