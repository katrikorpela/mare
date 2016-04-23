\name{mare-package}
\alias{mare-package}
\alias{mare}
\docType{package}
\title{
Microbiota Analysis in R Easily
~~ mare ~~
}
\description{
This package processes 16S rRNA gene amplicon reads produced by Illumina MiSeq, visualises the data and performs some basic analyses. The package utilises the USEARCH software to process the reads, so for full functionality, USEARCH version 8 should be installed (http://drive5.com/). 
}
\details{
\tabular{ll}{
Package: \tab mare\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-04-13\cr
License: \tab What license is it under?\cr
}
The package offers an accessible, easy-to-use pipeline for bacterial community analysis based on sequenced 16S amplicons. The package has functions to run the analysis from the raw 16S reads to p-value tables and figures. The package formats a reference database into the UTAX-compatible UDB-format and to optionally limits the reference database to gut-associated taxa, performs various optional preprocessing steps, OTU clustering, and taxonomic annotation of the OTUs or the reads themselves. Taxonomic tables at different phylogenetic levels from phylum to genus are produced. The results are visualised and bacterial taxa significantly associated with variables of interest are identified. All functions read data in and save the results as text files, not as R objects, so large files are not stored in the R session. Detailed instructions and examples can be found in the help files of the functions and in the user guide.
}
\author{
Katri Korpela

Maintainer: Katri Korpela <katri.korpela@helsinki.fi>
}
\references{
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,  dx.doi.org/10.1038/nmeth.2604].

See http://drive5.com/usearch/manual/ for details on reference database construnction, taxonomic annotation, paired-end read merging, quality filtering, and OTU clustering.
}

