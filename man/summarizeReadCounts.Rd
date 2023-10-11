\name{summarizeReadCounts}
\alias{summarizeReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes normalized or raw aligned read counts from the BAM file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The summarizeReadCounts() function calculates normazlied or raw aligned read counts for each variant in the BAM file. Normalized aligned read counts are computed by dividing the raw read counts by the total number of reads in the sample and then multiplying by the scaling factor.
}
\usage{
summarizeReadCounts(input, scale=1e+7, chrom=c(1:22,'X','Y'), output=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{x} here~~
  A single character string of the name of the BAM file.
  }
  \item{scale}{
  A numeric value of the scaling factor. If `scale` is set to 0, the raw aligned read counts are computed. If `scale` is negative, the normalized aligned read counts are computed using the default scaling factor.
}
  \item{chrom}{
  A vector of chromosome names specifying which chromosome to include. If `chrom` is set to NULL, all chromosomes are include.
}
  \item{output}{
  A single character string of the output file name. By default, no output file is written.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.table in the RLE format with columns of *Chr*, *Length*, and *NormRC* or *RC*.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Jingqi Duan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## BAM file from the ENCODE eCLIP-seq experiment targeting HNRNPK in K562 cell line.
experiment = 'https://www.encodeproject.org/files/ENCFF525XIM/@@download/ENCFF525XIM.bam'

## Normalized aligned read counts
normRC = summarizeeReadCounts(experiment, chrom=c(13,14))

## Raw aligned read counts
RC = summarizeeReadCounts(experiment, scale=0, chrom=c(13,14))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
