\name{PCoA}
\alias{PCoA}
\title{
Principal coordinates analysis
}
\description{
The function conducts principal coordinates analysis using Bray-Curtis dissimilarities, interpolates values of a continuous variable into the principal coordinate space, and plots the PCoA with the interpolated continuous variable as the background color. Group differences are shown by different-coloured symbols and the proportion of overall variation explained by group is printed on the bottom of the figure with a p-value from adonis in package vegan.
}
\usage{
PCoA(taxonomic.table, meta, readcount.cutoff = 0, group = NULL, group2 = NULL,
    components = c(1, 2), background.variable = NULL, colour.scheme = "terrain.colors", 
    ellipse = F, hull = F, spider = F, legendplace = "topright", 
    select.by = NULL, select = NULL, relative = F, pdf = F, keep.result = F,
    subjectID = NULL, time = NULL, distance = "bray", constrain = NULL)
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
  \item{group2}{
Name of a second grouping variable.
}
  \item{components}{
Which principal components should be shown? Defaults to the first two components.
}
  \item{background.variable}{
Name of the continuous variable from the metadata file or a bacterial taxon to be shown as the background colour.
}
  \item{colour.scheme}{
Backround colour scheme. Can be a vector of colour codes or names, or pre-specified colour scheme such as "terrain.colors" (default), "heat.colors","topo.colors", or "cm.colors". NOTE: These must be in quotes!
}
  \item{ellipse}{
Draw standard deviation ellipse for each group using function ordiellipse from package vegan.
}
  \item{hull}{
Draw outline for each group using function ordihull from package vegan.
}
  \item{spider}{
Draw lines joining the samples within each group using function ordispider from package vegan.
}
  \item{legendplace}{
Place of the group legend in the figure. Options are "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value of the selection variable will be selected for plotting.
}
  \item{relative}{
Should relative abundances rather than absolute counts be used? If yes, specify TRUE.
}
\item{pdf}{
Should the figure be saved as pdf? If yes, specify TRUE.
}
\item{keep.result}{
Should the PC-component scores be returned? If yes, specify TRUE.
}
\item{subjectID}{
Column name for subject ID in the metadata file, if several samples from the same subjects are present and should be linked in the plot.
}
\item{time}{
Column name for time/order variable in the metadata file, if several samples from the same subjects are present and should be linked in the plot.
}
\item{distance}{
Use Bray-Curtis distance ("bray"), Pearson correlation distance ("pearson"), Euclidean distance ("euc").
}
\item{constrain}{
To constrain the first principal component to be a specified variable, specify a variable name.
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
#richness as the background colour.	Save the plot as PDF. Save the result as R object "pctable".

pctable <- PCoA(taxonomic.table = 'organised_genus_table.txt',
                meta = 'meta.txt', 
                group = 'Treatment', 
                background.variable = 'Richness', 
                readcount.cutoff = 2000,
                pdf = T, keep.results = T)
	
}
}

