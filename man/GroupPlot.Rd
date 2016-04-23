\name{GroupPlot}
\alias{GroupPlot}
\alias{GroupPlotGUI}
\title{
Visualise group differences
}
\description{
The function prints five different plots to visualise group differences in selected bacterial taxa.
}
\usage{
GroupPlot(taxa, group = NULL, taxonomic.table, meta, readcount.cutoff = 0, stacked = T, 
bar = T, box = T, bean = T, covariate = NULL, smooth.method = 'loess',select.by = NULL, 
select = NULL, pdf = F, quartz = T, label.direction = 1)
}
\arguments{
  \item{taxa}{
Names of the taxa that should be visualised.
}
  \item{group}{
Grouping variable.
}
  \item{taxonomic.table}{
Taxonomic table file containing the taxa of interest. Should be the name of a text file.
}
  \item{meta}{
Metadata file containing the grouping variable. Should be the name of a text file.
}
  \item{readcount.cutoff}{
Lowest acceptable read count. Samples with fewer reads will be ignored.
}
  \item{stacked}{
Should a stacked barplot be printed? If yes, specify TRUE.
}
  \item{bar}{
Should a barplot of the means and standard errors with be printed? If yes, specify TRUE.
}
  \item{box}{
Should a box-and-whiskers plot be printed? If yes, specify TRUE.
}
  \item{bean}{
Should a bean plot be printed? If yes, specify TRUE.
}
  \item{covariate}{
Continuous variable to show associations with the bacterial taxa in the different groups.
}
  \item{smooth.method}{
Method to calculate the association between the taxa and the covariate. Set to 'lm' for a linear fit line, and 'loess' for a smooth fit line.
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected.
}
  \item{pdf}{
Should the figures be printed as pdf? If yes, specify TRUE.
}
  \item{quartz}{
Plot the figure into a quartz device (default). If your OS does not support quartz, set to FALSE.
}
  \item{label.direction}{
Direction of the x-axis labels: 1 = horizontal,  2 = vertical.
}
}

\author{
Katri Korpela
}

\seealso{
\code{\link{CommonTaxa}}
}
\examples{
\dontrun{	
#Make all five plots, showing the association between the most abundant
#genus-level bacterial taxa and treatment group, including a plot showing how the 
#taxa correlate with OTU richness in the different treatment groups, omitting samples 
#with <2000 reads.

common <- CommonTaxa(taxonomic.table = 'organised_genus_table.txt', mean.abundance = 0.05, prevalence = 10)
GroupPlot(taxa = unlist(common[1]),
          group = 'Treatment',
          taxonomic.table = 'organised_genus_table.txt', 
          meta = 'meta.txt', 
          readcount.cutoff = 2000,  
          covariate = 'Richness', 
          smooth.method = 'lm')	
}
}
