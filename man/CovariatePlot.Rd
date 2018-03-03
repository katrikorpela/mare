\name{CovariatePlot}
\alias{CovariatePlot}
\alias{CovariatePlotGUI}
\title{
Plot trend lines of bacterial abundances against a continuous covariate
}
\description{
The function plots the relative abundances of selected bacterial taxa against a continuous covariate.
}
\usage{
CovariatePlot(meta, taxonomic.table, covariate, taxa, smooth.method = "loess", 
    readcount.cutoff = 0, select.by = NULL, select = NULL, pdf = F, logtrans = F, nozeros = F)
}

\arguments{
\item{meta}{
Name of the metadata file containing the covariate. Should be the name of a text file.
} 
 \item{taxonomic.table}{
Name of the taxonomic table. Should be the name of a text file.
}
  \item{covariate}{
The covariate to plot against the bacterial taxa. Should be the name of a column in the metadata file.
}
\item{taxa}{
A vector specifying the bacterial taxa that should be plotted.
}  
\item{smooth.method}{
Smoothing method to use. Options include 'lm' for a linear fit, and 'loess' for locally-weighted polynomial regression fit.
}
  \item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected.
}
   \item{pdf}{
Should the figure be saved as pdf? If yes, specify TRUE.
}
  \item{logtrans}{
Should the bacterial abundances be log-transformed? TRUE or FALSE.
}
  \item{logtrans}{
Should zero values be omitted from the plot? TRUE or FALSE.
}
}
\author{
Katri Korpela
}
\examples{
\dontrun{
#Plot the most abundant and prevalent genera against age, excluding samples with <2000 reads and 
#selecting only cases in treatment group 1.

CovariatePlot(meta = 'meta.txt',
              taxonomic.table = 'organised_genus_table.txt', 
              readcount.cutoff = 10000, 
              covariate = 'Age', 
              taxa = unlist(CommonTaxa(taxonomic.table = 'organised_genus_table.txt', 
              mean.abundance = 0.05, prevalence = 1)), 
              smooth.method = 'loess',
              select.by = 'Treatment', 
              select = '1')
	
}
}

