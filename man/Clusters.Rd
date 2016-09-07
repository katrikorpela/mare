\name{Clusters}
\alias{Clusters}
\title{
Clustering of bacterial taxa
}
\description{
Hierarchical clustering based on significant correlations between bacterial taxa.
    }
\usage{
Clusters(taxonomic.table, meta, N.taxa = NULL, readcount.cutoff = 0,
          minimum.correlation = 0.5,  minimum.network = 1,  cluster.similarity = 1, quartz = T, pdf = F)
}

\arguments{
  \item{taxonomic.table}{
Name of the taxonomic table. Should be the name of a text file.
}
  \item{meta}{
Name of the metadata file containing the grouping variable. Should be the name of a text file.
}
  \item{N.taxa}{
The number of most abundant taxa to include.
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
  \item{cluster.similarity}{
Minimum cluster similarity, i.e. maximum distance between taxa to be included in the same cluster.
}
\item{pdf}{
Should the figure be saved as pdf? If yes, specify TRUE.
}
  \item{quartz}{
Plot the figure into a quartz device (default). If your OS does not support quartz, set to FALSE.
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
#Split taxa with ward-distance > 1 into different clusters.

Clusters(taxonomic.table = "organised_genus_table.txt", 
          meta="meta.txt", N.taxa = 50, readcount.cutoff = 0,
          minimum.correlation = 0.5,  minimum.network = 1,  cluster.similarity = 1, quartz = T, pdf = F)

	
}
}

