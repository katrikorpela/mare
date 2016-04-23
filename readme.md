# mare

https://zenodo.org/badge/21622/katrikorpela/mare.svg

The mare R package is an easy-to-use pipeline for microbiota analysis based on 16S-amplicon reads. It takes the raw reads, creates taxonomic tables, visualises the results, and finally identifies organisms significantly associated with variables of interest. For OTU clustering and taxonomic annotation the package relies on USEARCH (at least version 8.1.1756_i86osx32, possibly also later versions), which you need to obtain first.

For more information see the files mareGuide and mareBackground. To open the files, when you have mare installed, you can run: 

browseURL(system.file("mareGuide.pdf", package = "mare"))

browseURL(system.file("mareBackground.pdf", package = "mare"))

# Installation

First install the required packages:

install.packages("devtools", "R2admb", "R2admb", "vegan", "sp", "gstat", "Hmisc", "beanplot", "stringr", "MASS", "seqinr", "ggplot2", "reshape2")

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

# Citation

Katri Korpela (2016). mare: Microbiota Analysis in R Easily. R package version 1.0. https://github.com/katrikorpela/mare

# Databases

Some reference databases (Silva, RDP, CAZy for carbohydrate utilisation prediction) come with the package to to make getting started easier. To find location of the database, type:
filepath <- system.file("extdata", "NameOfTheDatabase", package="mare")

List of databases:


silva_full.fasta (full-length Silva-database)

silva_full.udb (udb-formatted full-length Silva-database)



AA.txt ("Auxiliary Activities" database from CAZy)

CE.txt ("Carbohydrate Esterases" database from CAZy)

GH.txt ("Glycoside Hydrolases" database from CAZy)

GT.txt ("Glycosyl Transferases" database from CAZy)

PL.txt ("Polysaccharide Lyases" database from CAZy)
