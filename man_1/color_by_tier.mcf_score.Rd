\name{color_by_tier.mcf_score}
\alias{color_by_tier.mcf_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
color_by_tier.mcf_score(mcf, tier = Inf, palette = rainbow)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcf}{
%%     ~~Describe \code{mcf} here~~
}
  \item{tier}{
%%     ~~Describe \code{tier} here~~
}
  \item{palette}{
%%     ~~Describe \code{palette} here~~
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mcf, tier = Inf, palette = rainbow) 
{
    res <- do10_1(mcf[, 2]/sum(mcf[, 2], na.rm = TRUE))
    r <- res$x49_2
    r[r > tier] <- tier + 1
    r <- data.frame(mcf, color = palette(max(r))[r])
    class(r) <- c("mcf_score", "data.frame")
    ref(r) <- ref(mcf)
    r
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
