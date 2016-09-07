\name{PCoA}
\alias{PCoA}
\title{
Principal coordinates analysis
}
\description{
The function conducts principal coordinates analysis using Bray-Curtis dissimilarities, interpolates values of a continuous variable into the principal coordinate space, and plots the PCoA with the interpolated continuous variable as the background color. Group differences are shown by different-coloured symbols and the proportion of overall variation explained by group is printed on the bottom of the figure with a p-value from adonis in package vegan.
}
\usage{
PCoA(taxonomic.table, meta, readcount.cutoff = 0, group = NULL, 
    components = c(1, 2), background.variable = NULL, legendplace = "topright", 
    select.by = NULL, select = NULL, relative = F, pdf = F, quartz = T)
}

\arguments{
  \item{taxonomic.table}{
Name of the taxonomic table. Should be the name of a text file.
}
  \item{meta}{
Name of the metadata file containing the grouping variable. Should be the name of a text file.
}
  \item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}
  \item{group}{
Name of the grouping variable.
}
  \item{components}{
Which principal components should be shown? Defaults to the first two components.
}
  \item{background.variable}{
Name of the continuous variable from the metadata file or a bacterial taxon to be shown as the background colour.
}
  \item{legendplace}{
Place of the group legend in the figure. Options are "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected for plotting.
}
  \item{relative}{
Should relative abundances rather than absolute counts be used? If yes, specify TRUE.
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
#Plot the first two principal coordinates in the genus-level data, ignoring samples with 
#fewer than 2000 reads, showing differences between treatment groups, and the OTU 
#richness as the background colour.	Save the plot as PDF.

pctable <- PCoA(taxonomic.table = 'organised_genus_table.txt',
                meta = 'meta.txt', 
                group = 'Treatment', 
                background.variable = 'Richness', 
                readcount.cutoff = 2000,
                pdf = T)
	
}
}

