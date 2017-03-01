\name{CombineProjects}
\alias{CombineProjects}

\title{
Combine taxonomic tables and metadata of two projects
}
\description{
Combines the organised taxonomic tables and metadata files of two projects that were processed separately, and (optionally) removes batch effects (mean-level differences between projects).
}
\usage{
CombineProjects(project1folder,project2folder, meta1, meta2, batch.effect = F, all.taxa = T)
}

\arguments{
  \item{project1folder}{
Path to the folder where data files of project 1 are located.
}
  \item{project2folder}{
Path to the folder where data files of project 2 are located.
}
  \item{meta1}{
Name of the metadatafile for project 1.
}
  \item{project1folder}{
Name of the metadatafile for project 2.
}
  \item{batch.effect}{
Set to TRUE, if you want to remove a batch effect, i.e. to standardise the mean relative abundances of the taxa.
}
  \item{all.taxa}{
Set to TRUE, if you want to include all taxa that appear at least in one of the projects. The default (FALSE) includes only taxa that are observed in both projects.
}
}
\details{
This function combines the organised taxonomic tables and metadata files. The function assumes that the metadata files  are in the same project folders (project1folder and project2folder). 
}

\author{
Katri Korpela
}

\examples{
\dontrun{
CombineProjects(project1folder="/path/to/project1",
                project2folder="/path/to/project2",
                meta1="meta1.txt",meta2="meta2.txt")

}}
