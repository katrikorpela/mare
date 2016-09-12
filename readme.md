# mare

The mare R package is an easy-to-use pipeline for microbiota analysis based on 16S-amplicon reads. It takes the raw reads, creates taxonomic tables, visualises the results, and finally identifies organisms significantly associated with variables of interest. For OTU clustering and taxonomic annotation the package relies on USEARCH (at least version 8.1.1756_i86osx32, possibly also later versions), which you need to obtain first.

For more information see the files mareGuide and mareBackground. To open the files, when you have mare installed, you can run: 

browseURL(system.file("mareGuide.pdf", package = "mare"))

browseURL(system.file("mareBackground.pdf", package = "mare"))

# Citation

Katri Korpela (2016). mare: Microbiota Analysis in R Easily. R package version 1.0. https://github.com/katrikorpela/mare
[![DOI](https://zenodo.org/badge/21622/katrikorpela/mare.svg)](https://zenodo.org/badge/latestdoi/21622/katrikorpela/mare)


# Installation

First install the required packages:

install.packages(c("devtools", "R2admb", "R2admb", "vegan", "sp", "gstat", "Hmisc", "beanplot", "stringr", "MASS", "seqinr", "ggplot2", "reshape2", "qgraph", "gplots"))

If you want to use the graphical user interface, install also the package fgui:
install.packages("fgui")

If you want to use the Blast option, install also the package CHNOZS:
install.packages("CHNOZS")

If some of the packages cannot be installed, install them manually (see R help pages on how to install packages). 

Then continue to run this code to install packages from Bioconductor and R-forge:

devtools::source_url("https://bioconductor.org/biocLite.R")

biocLite(c("Biostrings", "ShortRead","BiocGenerics"))

install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")), type="source")

When all packages are installed, you can install mare:

devtools::install_github("katrikorpela/mare")



# Functions

Blast: Creates a BLAST-based taxonomic table

CAZy: Does carbohydrate enzyme abundance predictions based on the species table

Clusters: Performs clustering of the bacterial taxa and plots a correlation network

ChangeTest: Tests for differences between groups or associations with covariates in the change of bacterial abundances from one time point to another

CommonTaxa: Identifies the most abundant and common taxa in the dataset

CopyFiles: Copies the taxonomic tables and metadata to a new folder for analysis

CorrelationMap: Plots a heatmap of correlations between bacterial taxa and selected variables in the metadata file

CovariatePlot: Plots the selected bacterial taxa against a covariate

CovariateTest: Tests for associations between the bacterial taxa and a covariate

FormatRefDB: Formats a reference database from fasta to UDB-format

GroupPlot: Plots group comparisons

GroupTest: Tests for differences between groups in bacterial abundances

mareGUI: Graphical user interface for mare

Organise: Organises the taxonomic tables based on the metadata file

PCoA: Pricipal Coordinates Analysis

ProcessReads: Processes the sequencing reads

Simple: A short-cut function to format the database, process the reads and create taxonomic tables

TaxonomicTable: Creates taxonomic tables

# Databases

Some reference databases (Silva, and CAZy for carbohydrate utilisation prediction) come with the package to to make getting started easier. To find location of the database, type:
filepath <- system.file("extdata", "NameOfTheDatabase", package="mare")

List of databases:

silva_full.fasta (full-length Silva-database)

silva_full.udb (udb-formatted full-length Silva-database)

AA.txt ("Auxiliary Activities" database from CAZy)

CE.txt ("Carbohydrate Esterases" database from CAZy)

GH.txt ("Glycoside Hydrolases" database from CAZy)

GT.txt ("Glycosyl Transferases" database from CAZy)

PL.txt ("Polysaccharide Lyases" database from CAZy)
