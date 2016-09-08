\name{CopyFiles}
\alias{CopyFiles}
\title{
Copy the taxonomic tables and move the R session to the new directory
}
\description{
Creates a new directory, copies the taxonomic tables, OTU table, metadata file and the readnumber file to the new directory, and moves the R session to the new directory. This leaves the processed reads and taxonomic tables intact where they are and allows the further analysis steps to be done in clean new folder.
}
\usage{
CopyFiles(folders = "*TaxonomicTables", files = "*table.txt", meta = NULL, folder.name = "",
    copy.to)
}
\arguments{
  \item{folders}{
Names of the directories where the taxonomic tables to be copied are. Defaults to all folders in the current working directory with names ending 'TaxonomicTables'.
}
  \item{files}{
List of files to be copied. Defaults to all files in the specified folders with names ending '_table.txt'.
}
  \item{meta}{
Metadata file. Should be the name of a text file.
}
  \item{folder.name}{
Folder name that was given to ProcessReads and Taxonomictable functions.
}
  \item{copy.to}{
Directory that should be created (not an existing folder!).
}
}
\author{
Katri Korpela
}

\examples{
	\dontrun{
#Create a new folder called NewFolder and copy all tables from the folder 
#ExampleTaxonomicTables into the new folder.
CopyFiles(folders = 'ExampleTaxonomicTables', 
          files = '*table.txt', 
          copy.to = '/new/folder/directory/NewFolder')

}
}
