\name{plot.mcf_score}
\alias{plot.mcf_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.mcf_score(mcf, genes = NULL, shrink = 0.7, pch = 21, point_cex = 0.5, text_cex = 0.5, lb = 1/60, lb1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcf}{
%%     ~~Describe \code{mcf} here~~
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{shrink}{
%%     ~~Describe \code{shrink} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{point_cex}{
%%     ~~Describe \code{point_cex} here~~
}
  \item{text_cex}{
%%     ~~Describe \code{text_cex} here~~
}
  \item{lb}{
%%     ~~Describe \code{lb} here~~
}
  \item{lb1}{
%%     ~~Describe \code{lb1} here~~
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
function (mcf, genes = NULL, shrink = 0.7, pch = 21, point_cex = 0.5, 
    text_cex = 0.5, lb = 1/60, lb1 = 0) 
with(ref(mcf), {
    x90_4 <- FALSE
    x60_5 <- mcf[, 2]
    x60_4 <- mcf[, 3]
    x48 <- data.frame(mcf[, 1], x60_5, x60_4)
    x50 <- genes[, c(1, 3)]
    x50_12 <- genes[, 2]
    x190 <- data.frame(gene[x48[, 1], 2:4], matrix(NA, nrow(x48), 
        3))
    x190[, 6] <- (x190[, 2] + x190[, 3])/2
    x48 <- data.frame(x190, x48)
    x48 <- x48[order(x48[, 1], x48[, 6]), ]
    x190 <- x48[, 1:6]
    x48 <- x48[, -(1:6)]
    x200 <- x48[, 2]
    x300 <- x48[, 3]
    m0 <- max(x200)
    n0 <- min(x200[x200 != Inf])
    x53 <- match(sort(unique(x190[, 1])), chr[, 4])
    x54 <- match(chr[, 4], sort(unique(x190[, 1])))
    x23 <- cumsum(0 + chr[x53, 1])
    x25 <- c(0, head(x23, -1))
    x24 <- 1 + (x25 + x23)/2
    x52 <- chr[x53, 5:6]
    x52 <- x25 + (x52[, 1] + x52[, 2])/2
    if (length(x50) > 0) {
        x50_1 <- x25[x54[as.numeric(gene[x50[, 1], 2])]] + rowSums(gene[x50[, 
            1], 3:4], 2)/2
        if (!all(is.na(x50_1))) {
            x50 <- x50[!is.na(x50_1), ]
            x50_16 <- x50_12[!is.na(x50_1)]
            x50_1 <- x50_1[!is.na(x50_1)]
            x70_2 <- gene[x50[, 1], 2]
            x50_2 <- rep(0, nrow(gene))
            x50_2[x48[, 1]] <- x200
            x50_1 <- cbind(x50[, 1:2], x50_1, x50_2[x50[, 1]])
            x50_14 <- order(x50_1[, 3])
            x50_1 <- x50_1[x50_14, ]
            x70_2 <- x70_2[x50_14]
            x50_16 <- x50_16[x50_14]
        }
        else {
            x50_1 <- c()
        }
    }
    else {
        x50_1 <- c()
    }
    x210 <- x25[x54[as.numeric(x190[, 1])]]
    x210 <- x210 + x190[, 6]
    if (length(x50_1) > 0) {
        x50_13_1 <- par()$plt
        if (is.null(shrink)) {
            x50_16_1 <- par()$ps
            x50_16_2 <- max(sapply(x50_16, function(x) nchar(x) * 
                x50_16_1 * 0.5 * (1/72)))
            x50_16_3 <- par()$pin
            shrink <- ((x50_16_3[1] - x50_16_2)/1.05)/x50_16_3[1]
        }
        x50_13 <- c(x50_13_1[1], x50_13_1[1] + shrink * (x50_13_1[2] - 
            x50_13_1[1]), x50_13_1[3:4])
        par(plt = x50_13)
    }
    plot(NA, NA, xlim = c(n0, m0), ylim = c(tail(x23, 1), 1), 
        yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
    for (j1 in 1:length(x53)) {
        j <- x53[j1]
        if (!is.na(x52[j1])) {
            polygon(c(n0, n0, m0, m0), c(x52[j1], x23[j1], x23[j1], 
                x52[j1]), col = rgb(0.95, 0.95, 0.95), border = NA)
            lines(x200[x190[, 1] == j & x210 < x52[j1]], x210[x190[, 
                1] == j & x210 < x52[j1]], col = "black", lty = 3)
            lines(x200[x190[, 1] == j & x210 >= x52[j1]], x210[x190[, 
                1] == j & x210 >= x52[j1]], col = "black", lty = 3)
            lines(c(n0, m0), c(x23[j1], x23[j1]), col = rgb(0.5, 
                0.5, 0.5), lty = 1)
        }
        else {
            lines(x200[x190[, 1] == j], x210[x190[, 1] == j], 
                col = "black", lty = 3)
        }
    }
    axis(side = 2, at = ifelse(is.na(x52), x24, x52), labels = chr[x53, 
        4])
    for (j in unique(x300)) points(x200[x300 == j], x210[x300 == 
        j], col = j, bg = j, pch = pch, cex = point_cex)
    if (length(x50_1) > 0) {
        x50_13 <- c(n0, 1, m0 - n0, tail(x23, 1) - 1)
        x50_6 <- c(1, tail(x23, 1), n0, m0)
        x70_1 <- do18_6((x50_1[, 3] - x50_6[1])/(x50_6[2] - x50_6[1]), 
            x70_2, lb, lb1)
        par(xpd = TRUE)
        for (i1 in 1:nrow(x50_1)) {
            text(x50_13[1] + x50_13[3] * 1.05, x50_13[2] + x50_13[4] * 
                x70_1[i1], x50_16[i1], col = x50_1[i1, 2], adj = c(0, 
                0.5), cex = text_cex)
            lines(c(x50_13[1] + x50_13[3] * 1.05, x50_13[1] + 
                x50_13[3], x50_1[i1, 4]), c(x50_13[2] + x50_13[4] * 
                x70_1[i1], x50_1[i1, 3], x50_1[i1, 3]))
        }
        par(plt = x50_13_1)
    }
  })
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
