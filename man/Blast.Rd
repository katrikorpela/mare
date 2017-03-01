\name{Blast}
\alias{Blast}

\title{
Create taxonomic tables based on Blast alignment of the reads
}
\description{
The function BLASTS the reads againts the NCBI database, selects the best hit for each read, 
and creates taxonomic tables based on the best hits. A blast taxonomy file is created with the 
following columns: read name, GI of the best hit, similarity, and evalue. 

}
\usage{
Blast(finalreads, read.table, NCBItaxonomy, folder.name = "")
}

\arguments{
  \item{finalreads}{
Directory of the processed reads, a file called finalreads.fasta produced by ProcessReads.
}
  \item{read.table}{
Directory of the read table, produced by ProcessReads.
}
  \item{NCBItaxonomy}{
Directory of the NCBI taxonomy files.
}
  \item{folder.name}{
Name for the new folder where the taxonomic tables will be written. Should be the same as the folder.name given to function ProcessReads in order to analyse the correct reads.
}
}

  \details{
To run Blast, you need to install the BLAST+ command line application installed 
<ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/>, and 
the NCBI taxonomy files  <ftp://ftp.ncbi.nih.gov/pub/taxonomy/>. You also need to install the R package CHNOSZ.
  }

\examples{
	\dontrun{

Blast(finalreads = "finalreads.fasta", read.table = "read_table.txt", NCBItaxonomy = "directory/to/NCBIdb/", folder.name = "")
}
}

\author{
Katri Korpela
}
