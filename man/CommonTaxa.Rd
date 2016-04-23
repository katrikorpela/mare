\name{CommonTaxa}
\alias{CommonTaxa}

\title{
Identify the most abundant and prevalent taxa
}
\description{
Creates a list with two components, 'Abundant' and 'Prevalent', listing the names of the most abundant and most prevalent taxa. This is the only function in the package mare that saves the output as an R object (and therefore has to be assigned a name).
}
\usage{
CommonTaxa(taxonomic.table, mean.abundance, prevalence)
}

\arguments{
  \item{taxonomic.table}{
Name of the taxonomic table. Should be a text file.
}
  \item{mean.abundance}{
Relative abundance cut-off: taxa with higher relative abundances are shown. Should be a value between 0-1.
}
  \item{prevalence}{
Prevalence cut-off: taxa occurring in a higher proportion of samples will be shown. Should be a value between 0-1.
}
}
\details{
}
\value{
A list of the most abundant and prevalent taxa.
 \item{Abundant}{The most abundant taxa in the samples.}
  \item{Prevalent}{The most prevalent taxa in the samples.}

}

\author{
Katri Korpela
}

\examples{
\dontrun{
#Identify the genus-level taxa whose mean relative abundance is >5 percent and those that 
#occur in >50 percent of all samples.
common <- CommonTaxa(taxonomic.table = 'organised_genus_table.txt',
                     mean.abundance = 0.05,
                     prevalence = 0.5)
abundant <- unlist(common[1])
prevalent <- unlist(common[2])
}}
