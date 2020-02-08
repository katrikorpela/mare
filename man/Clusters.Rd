\name{Clusters}
\alias{Clusters}
\title{
Clustering of bacterial taxa
}
\description{
Hierarchical clustering based on significant correlations between bacterial taxa.
    }
\usage{
Clusters(taxonomic.table, meta, N.taxa = NULL,  readcount.cutoff = 0,
                     minimum.correlation = 0.5, minimum.network = 1,
                      select.by = NULL, select = NULL, keep.result = F, pdf = F, relative = T)
}

\arguments{
  \item{taxonomic.table}{
Name of the taxonomic table. Should be the name of a text file.
}
  \item{meta}{
Name of the metadata file containing the grouping variable. Should be the name of a text file.
}
  \item{N.taxa}{
The number of most abundant taxa to include, if you don't want to include all taxa.
}
  \item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}
  \item{minimum.correlation}{
Minimum acceptable correlation between bacterial taxa. Lower values are considered 0.
}
  \item{minimum.network}{
Minimum network size, i.e. number of correlating taxa. Taxa with smaller networks are ignored.
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected.
}
  \item{keep.result}{
Should the results be returned as an R object? TRUE or FALSE.
}
\item{pdf}{
Should the figure be saved as pdf? TRUE or FALSE.
}
  \item{relative}{
Use relative abundances? Default is TRUE. FALSE will use data as is.
}
}

\value{
Returns a list with two components:
1. Data frame with relative abundances of the bacterial clusters
2. Cluster assignments of the bacterial taxa
}

\author{
Katri Korpela
}

\examples{
\dontrun{
#Select the 50 most abundant bacterial genera and calculate correlations between them.
#Set all correlation <0.5 to 0. 
#Include genera with >1 significantly correlating partners.

Clusters(taxonomic.table = "organised_genus_table.txt", 
          meta="meta.txt", N.taxa = 50, readcount.cutoff = 0,
          minimum.correlation = 0.5,  minimum.network = 1)

	
}
}

