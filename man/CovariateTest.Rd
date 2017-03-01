\name{CovariateTest}
\alias{CovariateTest}
\alias{GroupTest}
\alias{ChangeTest}
\title{
Significance Testing
}
\description{
CovariateTest estimates the significance of the association between each bacterial taxon and a continuous 
covariate.

GroupTest finds group differences in bacterial abundance.

ChangeTest finds differences between groups in microbial changes, or correlations between microbial changes and a covariate.

}
\usage{
CovariateTest(taxonomic.table, meta, covariate, readcount.cutoff = 0, 
    confounders = NULL, subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05,
    group = NULL, select.by = NULL, select = NULL, pdf = F, 
    min.prevalence = 0, min.abundance = 0, quartz = T, keep.result = F)

GroupTest(taxonomic.table, meta, group, compare.to = NULL, readcount.cutoff = 0, confounders = NULL, 
    subject.ID = NULL, outlier.cutoff = 3, p.cutoff = 0.05, 
    select.by = NULL, select = NULL, pdf = F,  min.prevalence = 0, min.abundance = 0, 
    label.direction = 1, quartz = T, keep.result = F)
    
ChangeTest(taxonomic.table, meta, group = NULL, compare.to = NULL, 
    covariate = NULL, readcount.cutoff = 0, confounders = NULL, subject.ID,  time,
    outlier.cutoff = 3, p.cutoff = 0.05, select.by = NULL, select = NULL, pdf = F, 
    consecutive = T, min.prevalence = 0, min.abundance = 0, label.direction = 1, quartz = T, keep.result =F)  
}
\arguments{
  \item{taxonomic.table}{
Taxonomic table containing the bacterial taxa to be tested. Should be the organised file returned by function Organise.
}
  \item{meta}{
Metadata file.
}
  \item{covariate}{
Covariate to test for associations with the bacterial taxa.
}
  \item{group}{
In group test, group is the grouping variable to test for associations with the bacterial taxa. In covariate test, group will be used to stratify the data and run the analyses separately in each group.
}
  \item{compare.to}{
Which level of the grouping variable should be used as the reference level?
} 
\item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}
  \item{confounders}{
Vector of maximally five potentially confounding variables to be included in the models. Should all be named in the metadata file.
}
  \item{subject.ID}{
If some samples are paired (e.g. from the same individual or family), specify which variable in the metadata file indicates the sample grouping (e.g. subject ID). In CovariateTest and GroupTest, a mixed model will be used, with the subject.ID as the random factor. In ChangeTest, this is required.
}
  \item{outlier.cutoff}{
Highest acceptable standard deviation from the overall mean abundance for a given taxon. Data points exceeding this value will be replaced by the cut-off value to eliminate disproportional effect of outliers on the results.
}
\item{p.cutoff}{
P-value cut-off for selecting which taxa are displayed in the figures.
}
  \item{min.prevalence}{
Minimum acceptable prevalence of a taxon. Taxa occurring less frequently will be ignored.
}
  \item{min.abundance}{
Minimum acceptable relative abundance of a taxon. Taxa occurring at a lower relative abundance will be considered absent when calculating prevalence.
}
  \item{select.by}{
Name of a variable in the metadata file by which a subset will be selected for plotting.
}
  \item{select}{
Determines which value on the selection variable will be selected.
}
  \item{pdf}{
Should be the result figure be printed as PDF? If yes, set to TRUE.
}
  \item{time}{
Name of the variable representing sampling time points.
}
  \item{consecutive}{
Do you want to calculate the change between consecutive time points (instead of calculating the change from baseline at all time points)? If yes, specify TRUE (default).
}
  \item{label.direction}{
Direction of the x-axis labels in plots: 1 = horizontal,  2 = vertical.
}
  \item{quartz}{
Plot the results in quartz device.
}
  \item{keep.result}{
Return the results to the R session.
}
}

\details{
Functions CovariateTest and GroupTest use generalised linear negative binomial models and fit two models to each taxon: one with all data points and one with only the datapoints with non-zero values. Sequencing data are commonly zero-inflated and the zeros may or may not indicate true absence of the taxon in the sample. Including the zeros may thus distort the observed association between a variable of interest and the bacterial taxa. It is also possible that the presence/absence of a taxon is determined by different factors than the abundance when the taxon is present. These problems are addressed by testing the association with and without zeros.

Violation against model assumptions (importantly the homogeneity of residuals) is commonly ignored and leads to meaningless p-values and potentially false conclusions. CovariateTest, GroupTest, and ChangeTest perform model checking on the negative binomial models, and if violation against model assumptions are observed, the model is discarded. In that case, a linear model on log-transformed relative abundances is fitted and tested. If this model fails to fullfill the assumptions,  no estimates or p-values are given, as these would be false.

It is possible to adjust for potentially confounding variables in the models. In that case, the results are presented as deviance from the expected abundance, where the expected value is based on the confounding variables. In ChangeTest, the models are always adjusted for baseline abundance of the taxon (due to unavoidable density dependence in the bacterial population dynamics) and therefore the results are always presented as deviance from expected change, rather than change. 

In ChangeTest, for bacterial taxa that are not observed in any sample from the same individual (change between time points cannot be reliably calculated), these individuals are excluded from the analysis.



}

\value{
Returns a table with the results if keep.result = T.
}

\author{
Katri Korpela
}

\examples{
\dontrun{
#Find genus-level taxa that are significantly associated with the categorical variable 
#Treatment, comparing the different treatments to the placebo group,
#visualise the association in the taxa whose p-values are <0.001, and save the figures as pdf, 
#and save the results into óbject called result.

result <- GroupTest(taxonomic.table = 'organised_genus_table.txt', 
          meta = 'meta.txt', 
          group = 'Treatment',
          compare.to = 'Placebo',
          p.cutoff = 0.001, 
          pdf = T,
          keep.result = T)

#Find genus-level taxa that are significantly associated with OTU richness, 
#and visualise the association in the taxa whose p-values are <0.05.

CovariateTest(taxonomic.table = 'organised_genus_table.txt', 
              meta = 'meta.txt', 
              covariate = 'Richness')

#Find genus-level taxa that are significantly associated with OTU richness, replacing values >3 
#standard deviations higher than the mean with a value = 3 standard deviations higher than the 
#mean, using only samples with >2000 reads, ignoring taxa that occur in <30 percent of the samples,
#and controlling for repeated measurements from the same individuals by specifying subject ID,
#and visualise the association in the taxa whose p-values are <0.05.

CovariateTest(taxonomic.table = 'organised_genus_table.txt', 
              meta = 'meta.txt', 
              covariate = 'Richness',
              readcount.cutoff = 2000,
              outlier.cutoff = 3, 
              min.prevalence = 0.3,
              subject.ID = 'ID')

#Find class-level taxa that are significantly associated with OTU richness in any treatment group,
#selecting only samples from time point 1.

CovariateTest(taxonomic.table = 'organised_class_table.txt', 
              meta = 'meta.txt', 
              covariate = 'Richness', 
              group = 'Treatment',
              select.by = 'Timepoint'
              select = '1')

#Find taxa whose abundances change differently in the different treatment groups compared to the placebo group

ChangeTest(taxonomic.table = 'organised_class_table.txt', 
              meta = 'meta.txt', 
              covariate = 'Richness', 
              group = 'Treatment',
              compare.to = 'Placebo',
              time = 'Timepoint') 
}
}	






