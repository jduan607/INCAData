\name{standardizeChr}
\alias{standardizeChr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Standardizes chromosome names for internal use
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The standardizeChr() function is used in the summarizePeakSignals() and summarizeReadCounts() functions. The function transforms the given chromosome names to the standard format beginning with "chr".
}
\usage{
standardizeChr(chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{
%%     ~~Describe \code{x} here~~
  A vector of chromosome names.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
A vector of chromosome names beginning with "chr".
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
