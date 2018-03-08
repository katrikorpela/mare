\name{Organise}
\alias{Organise}
\title{
Organise the taxonomic tables to match the metadata file
}
\description{
Rearranges the rows in the taxonomic table files to match those of the metadata file. This is necessary for all analyses which utilise both taxonomic data and metadata. Also calculates OTU richness and diversity (inverse Simpson, from package vegan), and richness of unique reads per genus, saves them into the metadata file, and adds the read count information into the metadata file. The metadata file should have one column with sample names matching exactly the sample names in the taxonomic tables. If you have several samples from the same individual or location, Organise will also calculate the Bray-Curtis dissimilarity and the Pearson correlation (similarity) in community composition between the samples.
}
\usage{
Organise(meta, sample.names, 
                     readtable = list.files(pattern = "annotated_read_table.txt$")[1],
                     otutable = list.files(pattern = "otutable.txt$")[1], 
                     taxonomic.tables = list.files(pattern = "_table.txt$")[-1], 
                     subject.ID = NULL, time = NULL)
}

\arguments{
  \item{meta}{
Name of a text file containing one column for the sample names, which must match exactly with the names of the samples in the taxonomic tables. For each sample, only one row should be specified. Otherwise the file can contain any additional information. Should be a text file.
}
  \item{sample.names}{
Name of the column in the metadata file which contains the sample names.
}
  \item{readtable}{
Name of the annotated read table. Should be a text file.
}
  \item{otutable}{
Name of the OTU table. Should be a text file.
}
  \item{taxonomic.tables}{
Names of the taxonomic tables. Should be text files.
}
 \item{subject.ID}{
If you have several samples from the same subject, specifying subject.ID and time will cause the 
function to calculate within-subject similarity/dissimilarity (Pearson correlation and Bray-Curtis distance based on the readtable) between the different time points.
}
 \item{time}{
If you have several samples from the same subject, specifying subject.ID and time will cause the 
function to calculate within-subject similarity/dissimilarity (Pearson correlation and Bray-Curtis distance based on the otutable) between the different time points.
}
}
\author{
Katri Korpela
}

\references{
Diversity index and dissimilarity between samples from the same individual is calculated using vegdist from package vegan.

}

\examples{
\dontrun{
Organise(meta = 'meta.txt', sample.names = 'sample_name', otutable = 'otutable.txt')
}
}

