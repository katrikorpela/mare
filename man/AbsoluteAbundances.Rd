\name{AbsoluteAbundances}
\alias{AbsoluteAbundances}


\title{
Translate the relative abundances into absolute abundances using total bacterial qPCR
}
\description{
The function takes the taxonomic tables and a file with qPCR result (copies of 16S per 1g of sample) and calculates the absolute abundances of the 16S copies present in the sample. Optionally, also 16S rRNA gene copy number correction is done resulting in an estimate of absolute abundances of genomes. 
}
\usage{
TaxonomicTable(qPCR.file, qPCR.var, folder.name = "", copynumberDB = NULL)
}

\arguments{
  \item{qPCR.file}{
Name of file in that contains a variable with the qPCR data. The samples should be in the same order as in the "readnumbers" file produced by ProcessReads.
}
  \item{qPCR.var}{
Name of the qPCR variable in the qPCR file.
}

  \item{folder.name}{
Name for the folder where the taxonomic tables are. Should be the same as the "folder.name" given to functions ProcessReads and TaxonomicTable.
}
  \item{copynumberDB}{
Optional. Path to a database with 16S rRNA gene copy number information for the bacterial taxa. Should have at least the variable "name", specifying the name of the organism and the variable "mean" specifying the mean copy number in the organism. See rrnDB for an example.
}
}

\references{
Stoddard S.F, Smith B.J., Hein R., Roller B.R.K. and Schmidt T.M. (2015) rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucleic Acids Research 2014; doi: 10.1093/nar/gku1201 [PMID:25414355]

}
\author{
Katri Korpela
}
