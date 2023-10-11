\name{summarizePeakSignals}
\alias{summarizePeakSignals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarizes peak signal values from the narrowPeak BED file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The summarizeReadCounts() function summarizes peak signal values in the narrowPeak BED file. Only regions with a p-value less than 0.05 under FDR control are included. If multiple files are provided, peak signal values are combined at the region level.
}
\usage{
summarizePeakSignals(input, suffix=NULL, chrom=c(1:22,'X','Y'), output=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{x} here~~
  A vector of the names of the narrowPeak BED files.
  }
  \item{suffix}{
  A vector of additional information about each file. By default, no additional information is provided. If `input` contains multiple files but `suffix` is not provided, `suffix` is set to the index of the files.
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
A data.table containing columns of *Chr*, *Start*, and *End* by default, and columns of *signalValue*, with one column for each file in the `input`.
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
## narrowPeak BED files from the ENCODE eCLIP-seq experiment targeting HNRNPK in K562 cell line.
experiment = c('https://www.encodeproject.org/files/ENCFF079MYZ/@@download/ENCFF079MYZ.bed.gz',
               'https://www.encodeproject.org/files/ENCFF769QNG/@@download/ENCFF769QNG.bed.gz')
suffix = c('Rep1','Rep2')

## Peak signal values from multiple files
signal1 = summarizePeakSignals(experiment, suffix, chrom=c(13,14))

## Peak signal values from one file
signal2 = summarizePeakSignals(experiment[1], chrom=c(13,14))
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
