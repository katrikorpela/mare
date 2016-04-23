\name{TaxonomicTable}
\alias{TaxonomicTable}
\alias{TaxonomicTableGUI}

\title{
Create taxonomic tables
}
\description{
The function assigns taxonomies (using USEARCH) for the OTUs or the unclustered reads, and creates taxonomic tables.
}
\usage{
TaxonomicTable(usearch.path, refDB, annotate.reads = T, 
folder.name = "", confidence.cutoff = 0.5)
}

\arguments{
  \item{usearch.path}{
Directory of the USEARCH program.
}
  \item{refDB}{
Reference database in UDB format.
}
  \item{annotate.reads}{
Should the reads be taxonomically annotated instead of the OTUs? TRUE is yes, FALSE if no.
}
  \item{folder.name}{
Name for the new folder where the taxonomic tables will be written. Should be the same as the folder.name given to function ProcessReads in order to analyse the correct reads.
}
  \item{confidence.cutoff}{
Lowest acceptable taxonomic confidence value.
}
}


\examples{
	\dontrun{
#Taxonomically annotate the reads (instead of the OTUs) that are in folder ExampleProcessedReads, 
#using 0.6 confidence cut-off and the gut-specifc reference database created by FormatRefDB. 
#Then create taxonomic tables, which are written into a new folder called ExampleTaxonomicTables.

TaxonomicTable(usearch.path = '/path/to/USEARCH8',
                          reference.database =  'GutDB.udb',
                          annotate.reads = T,
                          folder.name = 'Example',
                          confidence.cutoff = 0.6)

}
}
\references{

See http://drive5.com/usearch/manual/ for details on taxonomic annotation.

}
\author{
Katri Korpela
}
