\name{ProcessReads}
\alias{ProcessReads}
\title{
Pre-processing of raw 16S amplicon sequences
}
\description{
The function performs preprocessing of 16S sequences using USEARCH version 8: trimming of primer sequences, paired-end read merging or pooling, read truncation to set length, quality filtering, chimera filtering, dereplication, and OTU clustering.
}
\usage{
ProcessReads(forward.reads = list.files(pattern = "_R1_"), 
             reverse.reads = list.files(pattern = "_R2_"), 
             name.separator = "_", 
             usearch.path, 
             forward.primer = NULL, 
             reverse.primer = NULL,
             trim.beginning = 0, 
             merge = F, 
             min.merged.length = NULL, 
             min.overlap = NULL, 
             max.mismatches = NULL, 
             pool = F, 
             truncate = T, 
             filter = T, 
             trunc.quality = 2, 
             max.expected.error = 1, 
             readlength = NULL, 
             min.read.abundance = 0.00001, 
             folder.name = "")
}

\arguments{
  \item{forward.reads}{
List of the forward reads.
}
  \item{reverse.reads}{
List of the reverse reads (optional).
}
  \item{name.separator}{
Character that separates the sample name in the beginning of the sequence data files.
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
  \item{trim.beginning}{
Number of bases to be removed from the beginning of each read, after the primer has been trimmed (if provided). 
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
  \item{pool}{
Should the forward and reverse reads be pooled and treated as forward reads? TRUE if yes, FALSE if no.
}
  \item{truncate}{
Should the reads be truncated to remove low-quality end regions? TRUE if yes, FALSE if no.
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
The length to which all reads should be trimmed. Reads shorter than this will be discarded.
}
  \item{min.read.abundance}{
Minimum acceptable relative abundance of a unique read in the whole data set. Reads that occur less frequently are discarded. 
}
  \item{folder.name}{
Name for the new folder where the processed reads are written. 
}
}

\details{
ProcessReads has been tested with Illumina HiSeq and MiSeq reads, as well as 454 and IonTorrent reads. It requires as input .fastq files (or .fastq.gz), one file per sample. The quality scores are assumed to be Illumina quality scores, but if the quality scores are something else, the quality filtering step can be skipped in most cases, as filtering-by-readcount takes care of erroneous reads.

Sample names: The function splits the file names at "name.separator" and uses the beginning of the file name as the sample name. Each resulting sample name should be unique (so select "name.separator" carefully) and should not contain "_" or "/"!

There should not be "_" or " " anywhere in the path to the folder where the samples are.

Trimming the beginning of the read is done after primer trimming, so if your primer is 15nt and you want to remove the first 30nt(due to low quality), set trim.beginning to 15.
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
             readlength = 400, 
             folder.name = 'MergedFiltered')

ProcessReads(forward.reads = list.files(pattern='_R1_'),
             usearch.path = 'path/to/USEARCH8',
             forward.primer = c('AAATCGTACG'),
             truncate = T,
             trucation.length = 150,
             merge = F, 
             filter = F,
             readlength = 150, 
             min.read.abundance = 0.01,
             folder.name = 'ShortNotfiltered')
}
}

