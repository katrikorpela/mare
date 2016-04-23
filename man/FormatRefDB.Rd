\name{FormatRefDB}
\alias{FormatRefDB}
\title{
Create a UDB-format reference database
}
\description{
Creates a reference database in UDB-format, which is compatible with the USEARCH function 
for taxonomic annotation. Optionally, the reference database can be restricted to 
gut-associated taxa, and/or only to named species.
}
\usage{
FormatRefDB(refDB, usearch.path, confidence.file, gut.specific = F, 
    taxman = F, named.species = F)

}
\arguments{
  \item{refDB}{
The reference database in fasta format. This should either be the reference RDP database 
(http://drive5.com/usearch/manual/utax_downloads.html.) or a reference database from 
Taxman (http://www.ibi.vu.nl/programs/taxmanwww/).
}
  \item{usearch.path}{
Directory of the USEARCH program.
}
  \item{confidence.file}{
Taxonomy confidence file, which comes with the RDP reference data files at 
http://drive5.com/usearch/manual/utax_downloads.html. 
Select a suitable confidence file for your readlength.
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
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, 
Nature Methods [Pubmed:23955772,?? dx.doi.org/10.1038/nmeth.2604].

Brandt, B.W., Bonder, M.J., Huse, S.M. and Zaura, E. (2012) TaxMan: a server to trim rRNA 
reference databases and inspect taxonomic coverage. Nucleic Acids Research 40:W82-W87.
}
\author{
Katri Korpela
}

\examples{
	\dontrun{
#Create a gut-specific RDP reference database in UBD-format
FormatRefDB(refDB = '/path/to/the/RDPreferencedatabase/from/UTAXdownloads/at/drive5.com',
            usearch.path = 'path/to/USEARCH8',
            confidence.file = '/path/to/the/taxonomicconfidencefile/from/UTAXdownloads/at/drive5.com',
            gut.specific = T)
}
}