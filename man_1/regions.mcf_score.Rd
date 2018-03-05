\name{regions.mcf_score}
\alias{regions.mcf_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regions.mcf_score(mcf, num_regs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcf}{
%%     ~~Describe \code{mcf} here~~
}
  \item{num_regs}{
%%     ~~Describe \code{num_regs} here~~
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
function (mcf, num_regs) 
with(ref(mcf), {
    x42 <- data.frame(gene[mcf[, 1], 2:3], mcf[, 2], rep(1, nrow(mcf)))
    x42 <- x42[order(x42[, 1], x42[, 2]), ]
    x42_9 <- unique(x42[x42[, 3] > 0, 3])
    x42_9 <- x42_9[order(-x42_9)]
    x42_9 <- cbind(x42_9, matrix(0, length(x42_9), 3))
    for (i1 in 1:nrow(x42_9)) {
        x42_5 <- x42[, 3] >= x42_9[i1, 1]
        x42_6 <- c(x42_5[1], tail(x42[, 1], -1) == head(x42[, 
            1], -1) & !head(x42_5, -1) & tail(x42_5, -1) | tail(x42[, 
            1], -1) != head(x42[, 1], -1) & tail(x42_5, -1))
        x42_9[i1, 2:4] <- c(sum(x42_5), sum(x42_6), sum(x42_6 & 
            x42[, 4] > 0))
    }
    x90_6 <- x42_9[which(x42_9[, 4] >= min(num_regs, tail(x42_9[, 
        4], 1)))[1], 2]
    x60_5 <- mcf[, 2]/sum(mcf[, 2], na.rm = TRUE)
    x49_13 <- order(-x60_5)
    x49_13[x49_13] <- 1:length(x49_13)
    x60_1 <- x49_13 <= x90_6
    x49_2_1 <- cbind(x60_1, x60_5)
    x49_4 <- cbind(gene[mcf[, 1], 2:4], mcf[, 1], x49_2_1)
    x49_4 <- x49_4[order(x49_4[, 1], x49_4[, 2], x49_4[, 3]), 
        ]
    x49_4[, 7] <- NA
    x49_4[, 8] <- c(x49_4[1, 5] > 0, tail(x49_4[, 1], -1) == 
        head(x49_4[, 1], -1) & head(x49_4[, 5], -1) == 0 & tail(x49_4[, 
        5], -1) > 0 | tail(x49_4[, 1], -1) != head(x49_4[, 1], 
        -1) & tail(x49_4[, 5], -1) > 0)
    x49_4 <- x49_4[x49_4[, 5] > 0, ]
    x49_4[, 8] <- cumsum(x49_4[, 8])
    x49_4 <- x49_4[order(x49_4[, 8], -x49_4[, 6]), ]
    x49_5 <- which(!duplicated(x49_4[, 8]))
    x49_6 <- x49_4[x49_5, 4:6]
    x49_6 <- cbind(x49_6, tapply(rep(1, nrow(x49_4)), x49_4[, 
        8], sum), 1:nrow(x49_6))
    x49_6 <- x49_6[order(-x49_6[, 3]), ]
    x50 <- x49_6[, 1:2]
    x50_12 <- data.frame(idx = x50[, 1], leader = NA, score = x49_6[, 
        3], chr = gene[x50[, 1], 2], begin = NA, end = NA, genes = NA)
    for (i1 in 1:nrow(x50)) {
        if (x49_6[i1, 4] > 1) {
            x50_12[i1, "leader"] <- sprintf("\%s@\%s(+\%.0f)", gene[x50[i1, 
                1], 1], gene[x50[i1, 1], 2], x49_6[i1, 4] - 1)
        }
        else {
            x50_12[i1, "leader"] <- sprintf("\%s@\%s", gene[x50[i1, 
                1], 1], gene[x50[i1, 1], 2])
        }
        x49_7 <- x49_4[x49_4[, 8] == x49_6[i1, 5], c(4, 6)]
        x49_8 <- c(Inf, -Inf)
        x49_9 <- ""
        for (i2 in 1:nrow(x49_7)) {
            x49_8 <- c(min(x49_8[1], gene[x49_7[i2, 1], 3]), 
                max(x49_8[2], gene[x49_7[i2, 1], 4]))
            if (x49_9 == "") {
                x49_9 <- gene[x49_7[i2, 1], 1]
            }
            else {
                x49_9 <- sprintf("\%s\t\%s", x49_9, gene[x49_7[i2, 
                  1], 1])
            }
            if (i2 > 1 && x49_7[i2, 2] == x49_7[1, 2]) 
                x49_9 <- sprintf("\%s*", x49_9)
        }
        x50_12[i1, "begin"] <- x49_8[1]
        x50_12[i1, "end"] <- x49_8[2]
        x50_12[i1, "genes"] <- x49_9
    }
    x50_12
  })
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
