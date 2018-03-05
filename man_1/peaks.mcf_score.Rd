\name{peaks.mcf_score}
\alias{peaks.mcf_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
peaks.mcf_score(mcf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcf}{
%%     ~~Describe \code{mcf} here~~
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
function (mcf) 
with(ref(mcf), {
    options(stringsAsFactors = FALSE)
    x2 <- data.frame(mcf, ref(mcf)$gene[mcf$idx, ])
    x2 <- data.frame(x2[order(x2$chr, x2$begin), ], 1:nrow(x2))
    res <- do10_2(x2[, 2], x2$chr)
    x2 <- cbind(x2, tier = res$x49_2)
    x2 <- x2[x2[, 8] > 0, ]
    x2 <- cbind(x2, cumsum(c(TRUE, tail(x2$chr, -1) != head(x2$chr, 
        -1)) | c(FALSE, head(x2[, 7], -1) + 1 != tail(x2[, 7], 
        -1))))
    x2 <- x2[order(x2[, 9], -x2[, 2]), ]
    x3 <- x2[!duplicated(x2[, 9]), 2]
    x2 <- data.frame(x2, x3[x2[, 9]])
    x2 <- x2[x2[, 2] == x2[, 10], ]
    x5 <- table(x2[, 9])
    x5 <- ifelse(x5 > 1, paste("(+", x5 - 1, ")", sep = ""), 
        "")
    x6 <- sapply(split(x2$gene, x2[, 9]), function(x) paste(x, 
        collapse = ","))
    x7 <- paste(paste(x2$gene, x2$chr, sep = "@")[!duplicated(x2[, 
        9])], x5, sep = "")
    x8_1 <- tapply(x2[, 5], x2[, 9], min)
    x8_2 <- tapply(x2[, 6], x2[, 9], max)
    x4 <- data.frame(x2[!duplicated(x2[, 9]), c(1:2, 8, 4)], 
        begin = x8_1, end = x8_2, name = x7, genes = x6)
    x4 <- x4[order(x4[, 1]), ]
    x4
  })
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
