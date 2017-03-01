\name{HITChip2Seq}
\alias{HITChip2Seq}

\title{
Transform HITChip species-level data into taxonomic tables
}
\description{
Takes the species-level HITChip data table and transforms it into taxonomic tables in the same format as the output of TaxonomicTable.
}
\usage{
HITChip2Seq(HITChip_species_table)
}

\arguments{
  \item{HITChip_species_table}{
Name of the species-level HITChip data set.
}
}
\details{
The resulting taxonomic tables can be combined with data from sequencing projects using the function CombineProjects.
}

\author{
Katri Korpela
}
