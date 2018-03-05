\name{trfr.mcf_score}
\alias{trfr.mcf_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trfr.mcf_score(mcf, perms = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcf}{
%%     ~~Describe \code{mcf} here~~
}
  \item{perms}{
%%     ~~Describe \code{perms} here~~
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
function (mcf, perms = NULL) 
with(ref(mcf), {
    x2 <- data.frame(gene[mcf[, 1], 2:3], mcf[, 2], 1:nrow(mcf), 
        mcf[, 1])
    x2 <- x2[order(x2[, 1], x2[, 2]), ]
    res <- do10_2(x2[, 3], x2[, 1])
    x2 <- cbind(x2, res$x49_2)
    x5 <- unique(x2[, 3])
    x5 <- x5[order(-x5)]
    x5 <- cbind(x5, matrix(NA, length(x5), 3))
    x6 <- tail(x2[, 1], -1) == head(x2[, 1], -1)
    for (i1 in 1:nrow(x5)) {
        x7 <- x2[, 3] >= x5[i1, 1]
        x7_1 <- x7[-c(1, length(x7))]
        x7_h <- c(head(x7, 1), x7_1)
        x7_t <- c(x7_1, tail(x7, 1))
        if (!is.null(perms)) 
            x8 <- perms[, colnames(mcf)[2]] >= x5[i1, 1]
        x9 <- c(x7[1], x6 & !x7_h & x7_t | !x6 & x7_t)
        x5[i1, 2:4] <- c(sum(x7), sum(x9), ifelse(!is.null(perms), 
            sum(x8), NA))
    }
    x5 <- cbind(x5, x5[, 2]/nrow(x2), x5[, 4]/nrow(mcf))
    x5 <- cbind(x5, ifelse(x5[, 6] <= x5[, 5], x5[, 6]/x5[, 5], 
        1))
    x10 <- match(x2[, 3], x5[, 1])
    x2 <- cbind(x2, x5[x10, 2], x5[x10, 7], x5[x10, 3])
    x2 <- x2[order(x2[, 4]), c(5, 3, 6:9)]
    colnames(x2) <- c(colnames(mcf), paste(colnames(mcf)[2], 
        c("tier", "rank", "fdr", "regions"), sep = "_"))
    class(x2) <- c("mcf_score", class(x2))
    ref(x2) <- ref(mcf)
    x2
  })
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
