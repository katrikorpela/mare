\name{Organise}
\alias{Organise}
\title{
Organise the taxonomic tables to match the metadata file
}
\description{
Rearranges the rows in the taxonomic table files to match those of the metadata file. This is necessary for all analyses which utilise both taxonomic data and metadata. Also calculates OTU richness and diversity and saves them into the metadata file, and adds the read count information into the metadata file. The metadata file should have one column with sample names matching exactly the sample names in the taxonomic tables.
}
\usage{
Organise(meta, sample.names, otutable = list.files(pattern = "_otutable.txt$")[1], 
    taxonomic.tables = list.files(pattern = "_table.txt$"), subject.ID = NULL, 
    time = NULL)
}

\arguments{
  \item{meta}{
Name of a text file containing one column for the sample names, which must match exactly with the names of the samples in the taxonomic tables. For each sample, only one row should be specified. Otherwise the file can contain any additional information. Should be a text file.
}
  \item{sample.names}{
Name of the column in the metadata file which contains the sample names.
}
  \item{otutable}{
Name of the OTU table. Should be a text file.
}
  \item{taxonomic.tables}{
Names of the taxonomic tables. Should be text files.
}
 \item{subject.ID}{
If you have several samples from the same subject, specifying subject.ID and time will cause the 
function to calculate within-subject similarity/dissimilarity (Pearson correlation and Bray-Curtis distance based on the otutable) between the different time points.
}
 \item{time}{
If you have several samples from the same subject, specifying subject.ID and time will cause the 
function to calculate within-subject similarity/dissimilarity (Pearson correlation and Bray-Curtis distance based on the otutable) between the different time points.
}
}
\author{
Katri Korpela
}

\examples{
\dontrun{
Organise(meta = 'meta.txt', sample.names = 'sample_name', otutable = 'otutable.txt')
}
}

