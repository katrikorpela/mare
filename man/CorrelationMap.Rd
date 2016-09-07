\name{CorrelationMap}
\alias{CorrelationMap}
\title{
Heatmap of correlations between bacterial taxa and host variables
}
\description{
Calculates the correlations between bacterial taxa and selected variables in the metadata file, and plots the correlations as a heatmap.
  }
\usage{
CorrelationMap(taxonomic.table, meta, variables, readcount.cutoff = 0,
              select.by = NULL, selection = NULL,  outlier.cutoff = 3,  
               min.abundance = 0, min.prevalence = 0, pdf = F, quartz = T)
}

\arguments{
  \item{taxonomic.table}{
Name of the taxonomic table. Should be the name of a text file.
}
  \item{meta}{
Name of the metadata file containing the grouping variable. Should be the name of a text file.
}
  \item{variables}{
Names or index numbers of continuous variables in the metadata file to be included in the analysis.
}
  \item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}
 
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected for plotting.
}
  \item{outlier.cutoff}{
Highest acceptable standard deviation from the overall mean abundance for a given taxon. Data points exceeding this value will be replaced by the cut-off value to eliminate disproportional effect of outliers on the results.
}
\item{min.abundance}{
Minimum acceptable relative abundance of bacterial taxa. Used with min.prevalence to include only taxa with abundance > min.abundance in > min.prevalence samples.
}
\item{min.prevalence}{
Minimum acceptable prevalence of taxa (proportion of all samples). Used with min.abundance to include only taxa with abundance > min.abundance in > min.prevalence cases.
}
\item{pdf}{
Should the figure be saved as pdf? If yes, specify TRUE.
}
  \item{quartz}{
Plot the figure into a quartz device (default). If your OS does not support quartz, set to FALSE.
}
}

\value{
Returns a table with correlations of the bacterial taxa with the selected principal coordinates and p-values of the correlations.
}

\author{
Katri Korpela
}

\examples{
\dontrun{
#Plot the correlations between genus-level bacterial taxa and the first 15 variables in the metadata file.
#Select only taxa whose relative abundance is > 0.1 in > half of all samples.
#Save the plot as PDF.

CorrelationMap(taxonomic.table = 'organised_genus_table.txt',
                meta = 'meta.txt', 
                variables = c(1:15), 
                min.abundance = 0.01,
                min.prevalence = 0.5,
                pdf = T)
	
}
}

