\name{PathModel}
\alias{PathModel}

\title{
Test how bacterial taxa and host variables together influence the response variable
}
\description{
Combine information on bacterial taxa and a set of host variables into one model for the response variable. PathModel finds the bacterial taxa and host variables that together are significantly associated with the response variable, and the host variables that are associated with the response-variable-associated bacterial taxa.
}
\usage{
PathModel(taxonomic.table, meta, response = NULL, 
                      binomial.response = F, count.response = F, variables, readcount.cutoff = 0, 
                      outlier.cutoff = 3, 
                      select.by = NULL, select = NULL, pdf = F, 
                      min.prevalence = 0, min.abundance = 0, keep.result = F)
}

\arguments{
  \item{taxonomic.table}{
Taxonomic table containing the bacterial taxa to be tested. Should be the organised file returned by function Organise.
}
  \item{meta}{
Metadata file.
}
  \item{response}{
The response variable. This has to be a numerical variable: binary (0s and 1s), counts of observations, or continuous. If the response variable is omitted, the model will use all bacterial taxa as the response variables, and finds the best total model for each taxon.
}
  \item{binomial.response}{
Set to TRUE if the response variable contains only 0s and 1s.
}
  \item{count.response}{
Set to TRUE if the response variable is measured as counts of observations.
}
  \item{variables}{
Indices or names of the metadata file that you want to use as additional explanatory variables in the model. These must be numerical, categorical variables cannot be used.
} 
\item{readcount.cutoff}{
Lowest acceptable read count per sample. Samples with fewer reads are ignored.
}

  \item{outlier.cutoff}{
Highest acceptable standard deviation from the overall mean abundance for a given taxon. Data points exceeding this value will be replaced by the cut-off value to eliminate disproportional effect of outliers on the results.
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

  \item{quartz}{
Plot the results in quartz device.
}
  \item{keep.result}{
Store the results as a list. The first item on the list is the response model, and the second item is a list of all the bacterial models. The resulting final best models can then be accessed and inspected by using the summary command. 
}
}
\details{
PathModel first fits a linear model (make sure your variable is normally distributed) or a binomial model (for 0/1 response variables) or a negative binomial model (for count variables) separately to all bacterial taxa and host variables, then selects the best of those into a combined model and performs model reduction using AIC to simplify the model. It then performs the same model selection for each of the bacterial taxa retained in response variable model, now fitting negative binomial models. Then the best overall model is plotted, including the bacterial taxa and host variables associated with the response variable, and the host variables associated with these bacterial taxa.

}

\author{
Katri Korpela
}

