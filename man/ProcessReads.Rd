\name{ProcessReads}
\alias{ProcessReads}
\title{
Pre-processing of raw 16S amplicon sequences
}
\description{
The function trims primer sequences, and using USEARCH merges paired-end reads, quality-filters the reads, and truncates the reads to a uniform lenght. The reads are dereplicated and chimera-filtered, rare reads are removed, and OTU clustering is performed.  Most of the steps are optional.
}
\usage{
ProcessReads(forward.reads = list.files(pattern = "_R1_"), 
                         reverse.reads = list.files(pattern = "_R2_"), 
                         name.separator = "_", 
    usearch.path, forward.primer = NULL, reverse.primer = NULL, merge = F, 
    min.merged.length = NULL, min.overlap = NULL, max.mismatches = NULL, 
    truncate = F, trucation.length = NULL, 
    filter = F, trunc.quality = 2, max.expected.error = 1, min.readlength = NULL, 
    readlength = NULL, min.read.prevalence = 0.0001, folder.name = "")
}

\arguments{
  \item{forward.reads}{
List of the forward reads.
}
  \item{reverse.reads}{
List of the reverse reads (optional).
}
  \item{name.separator}{
Character that separates the sample name in the beginning of the read files.
}
  \item{usearch.path}{
Directory of the USEARCH program.
}
  \item{forward.primer}{
Primer sequence to be removed from the forward reads.
}
  \item{reverse.primer}{
Primer sequence to be removed from the reverse reads.
}
  \item{merge}{
Should the paired-end reads be merged? TRUE if yes, FALSE if no.
}
  \item{min.merged.length}{
(Merging) Minimum length of the merged paired-end read.
}
  \item{min.overlap}{
(Merging) Minimum overlap between the forward and reverse reads when merging.
}
  \item{max.mismatches}{
(Merging) Maximum number of mismatches in the overlapping region when merging.
}
  \item{truncate}{
Should the forward reads be truncated to remove low-quality end regions? TRUE if yes, FALSE if no.
}
  \item{trucation.length}{
(Truncation) Length of the truncated reads.
}
  \item{filter}{
Should the reads be quality filtered? TRUE if yes, FALSE if no.
}
  \item{trunc.quality}{
(Quality-filtering) Minimum quality score before the read is truncated.
}
  \item{max.expected.error}{
(Quality-filtering) Maximum expected error.
}
  \item{min.readlength}{
(Quality-filtering) Minimum acceptable read length.
}
  \item{readlength}{
The length to which all reads should be trimmed.
}
  \item{min.read.prevalence}{
Minimum acceptable prevalence of a sequence in the whole data set. Sequenes that occur less frequently are removed. 
}
  \item{folder.name}{
Name for the new folder where the processed reads are written. 
}
}

\references{
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,?? dx.doi.org/10.1038/nmeth.2604].

See http://drive5.com/usearch/manual/ for details on the paired-end merging and quality filtering.
}
\author{
Katri Korpela
}

\examples{
	\dontrun{
#Trim primer sequences from the reads and then merge paired-end reads requiring that the 
#minimum length of the merged reads is 400 bases, the minimum overlap between the forward 
#and reverse reads is 50 bases and the maximum number of mismathces in the overlap region 
#is 10. Then quality-filter the merged reads, requiring a minimum quality score of 8 (the 
#reads are truncated when quality scores are lower) and maximum expected error for the 
#reads to be 1, and remove all reads that are shorter than 400 bases. Then truncate all 
#reads to 400 bases, dereplicate the reads, and remove all rare reads (by default reads that 
#occur represent less than 0.01 percent of all reads). The processed reads are saved into a new folder called 
#MergedFilteredProcessedReads.	

ProcessReads(forward.reads = list.files(pattern = '_R1_'),
             reverse.reads= list.files(pattern = '_R2_'),
             usearch.path = 'path/to/USEARCH8',
             forward.primer = c('AAATCGTACG'),
             reverse.primer = c('GTCGGGATTTCCG'),
             min.merged.length = 400,
             min.overlap = 50,
             max.mismatches = 10,
             trunc.quality = 8,
             max.expected.error = 1, 
             min.readlength = 400,
             readlength = 400, 
             folder.name = 'MergedFiltered')

#Use only forward reads. Trim primer sequences from the reads and then truncate the reads to 
#150 bases to remove the low-quality end region.	Dereplicate the reads, and remove reads 
#that occur at less than 1 percent frequency. The processed reads are saved into a new 
#folder called ShortNotfilteredProcessedReads.	

ProcessReads(forward.reads = list.files(pattern='_R1_'),
             usearch.path = 'path/to/USEARCH8',
             forward.primer = c('AAATCGTACG'),
             truncate = T,
             trucation.length = 150,
             merge = F, 
             filter = F,
             readlength = 150, 
             min.read.prevalence = 0.01,
             folder.name = 'ShortNotfiltered')
}
}

