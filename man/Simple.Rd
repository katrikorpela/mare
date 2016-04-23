\name{Simple}
\alias{Simple}
\alias{SimpleGUI}
\title{
Pre-processing of raw 16S amplicon sequences and creation of taxonomic tables in one simple function
}
\description{
The function processes the forward reads. It trims the primer sequences and truncates the reads to 150 bases. The reads are dereplicated and rare reads (representing <0.1 percent of the reads) are removed. OTU clustering and taxonomic annotation of the reads is done using USEACH.
}
\usage{
Simple(forward.reads = NULL, forward.primer = NULL, name.separator = "_", 
    folder.name = "", usearch.path, refDB, gut.specific = F, 
    taxman = F, named.species = F)
}

\arguments{
  \item{forward.reads}{
List of the forward reads.
}
  \item{forward.primer}{
Primer sequence to be removed from the forward reads.
}
  \item{name.separator}{
Character that separates the sample name in the beginning of the read files.
}
  \item{folder.name}{
Name for the new folder where the processed reads are written. 
} 
  \item{refDB}{
The reference database in fasta format. This should either be the reference RDP database 
(http://drive5.com/usearch/manual/utax_downloads.html.) or a reference database from 
Taxman (http://www.ibi.vu.nl/programs/taxmanwww/).
}
  \item{usearch.path}{
Directory of the USEARCH program.
}
  \item{gut.specific}{
TRUE if you want a gut-specific database.
}
  \item{taxman}{
TRUE if your database is from Taxman.
}
  \item{named.species}{
TRUE if you want to restrict the database to named species.
}
}

\references{
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,?? dx.doi.org/10.1038/nmeth.2604].
}
\author{
Katri Korpela
}

\examples{
	\dontrun{

Simple(forward.reads = list.files(pattern = '_R1_'),
             usearch.path = 'path/to/USEARCH8',
             forward.primer = c('AAATCGTACG'),
             folder.name = 'Example')
}
}

