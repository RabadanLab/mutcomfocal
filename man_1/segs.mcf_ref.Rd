\name{segs.mcf_ref}
\alias{segs.mcf_ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segs.mcf_ref(r, cnv, mut)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{cnv}{
%%     ~~Describe \code{cnv} here~~
}
  \item{mut}{
%%     ~~Describe \code{mut} here~~
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
function (r, cnv, mut) 
{
    options(stringsAsFactors = FALSE)
    require(hash)
    require(data.table)
    x2 <- sort(unique(c(cnv[, 1], mut[, 1])))
    x2 <- data.frame(x2, 1:length(x2))
    x3 <- data.frame(cnv[, -1], cnv[, 1])
    x3[, 1] <- gsub("24", "Y", gsub("23", "X", x3[, 1]))
    x3 <- data.frame(x3[order(x3[, 1], x3[, 2], x3[, 3]), ], 
        1:nrow(x3))
    chr <- list()
    for (c in sort(unique(x3[, 1]))) {
        chr[[c]] <- list()
        chr[[c]]$x1 <- data.frame(0, x3[x3[, 1] == c, -1])
    }
    g.x9.2 <- data.frame(mut[, -1], mut[, 1], 1:nrow(mut))
    x9.3 <- sort(unique(g.x9.2[, 1]))
    for (c in names(chr)) chr[[c]] <- within(chr[[c]], {
        printf("\%s\n", c)
        x9.7 <- r$chr_list[[c]]
        x9.1 <- data.frame(paste(x9.7[, 4], c, sep = "@"), x9.7[, 
            4], x9.7[, 1], x9.7[, 2])
        x9.2 <- sort(unique(x9.1[, 1]))
        x9.3 <- g.x9.2
        x9.3 <- x9.3[x9.3[, 1] \%in\% x9.2, ]
        x9.4 <- data.frame(x9.3, x9.1[match(x9.3[, 1], x9.1[, 
            1]), -1])
        x9.4 <- data.frame(1, x9.4[, 4], x9.4, 1, x9.4[, 3])
        x9.4 <- x9.4[, c(1, 2, 4, 7:ncol(x9.4))]
        x9.5 <- data.frame(x9.1[, 1], c, x9.1[, 2:3])
        x2 <- data.frame(x1[, 2:3], 1, x1[, 7], 0 + (x1[, 5] < 
            0), abs(x1[, 5]), x1[, 6])
        x2.2.1 <- data.frame(paste(c, 0, x2[, 5], sep = "_"), 
            c, x2[, 2:3], x2[6:7])
        colnames(x2.2.1) <- 1:6
        x2.2.2 <- data.frame(paste(c, x9.4[, 1], x9.4[, 7], sep = "_"), 
            c, x9.4[, 4:5], 1 + x9.4[, 1], x9.4[, 6])
        colnames(x2.2.2) <- 1:6
        x2.2 <- rbind(x2.2.1, x2.2.2)
        x3.1 <- segsmerge(list(x2, x9.7))
        f <<- hash()
        x2_9 <- as.matrix(x2[, c(1:2, 5)])
        x2_10 <- as.matrix(x9.7[, 1:2])
        x3 <- lapply(x3.1, function(x) {
            x2_1 <- x[[2]][x[[2]][, "i"] == 1, "n"]
            x2_2 <- x[[2]][x[[2]][, "i"] == 2, "n"]
            if (length(x2_1) > 0 && length(x2_2) > 0) {
                x3_3 <- as.matrix(expand.grid(x2_1, x2_2))
                x3_3 <- unique(x3_3)
                x3_4 <- paste(x3_3[, 1], x3_3[, 2])
                x3_5 <- !has.key(x3_4, f)
                x3_4 <- x3_4[x3_5]
                x3_3 <- x3_3[x3_5, , drop = FALSE]
                x3_3 <- cbind(x3_3, x2_9[x3_3[, 1], , drop = FALSE], 
                  x2_10[x3_3[, 2], , drop = FALSE])
                x3_5 <- x3_3[, 5] == 1 | x3_3[, 5] != 1 & x3_3[, 
                  3] <= x3_3[, 6] & x3_3[, 7] <= x3_3[, 4]
                x3_4 <- x3_4[x3_5]
                if (length(x3_4) > 0) 
                  f[x3_4] <- 1
                as.data.frame(x3_3[x3_5, 1:2, drop = FALSE])
            }
        })
        x3 <- as.matrix(rbindlist(x3))
        x3 <- data.frame(x9.7[x3[, 2], c(4, 1, 2)], paste(c, 
            0, x2[x3[, 1], 4], sep = "_"), x2[x3[, 1], 5:7])
        setnames(x3, as.character(1:ncol(x3)))
        x3.2 <- data.frame(x9.4[, c(2, 4, 5)], paste(c, x9.4[, 
            1], x9.4[, 7], sep = "_"), x9.4[, 1] + 1, x9.4[, 
            c(6, 3)])
        colnames(x3.2) <- 1:ncol(x3.2)
        x7.0 <- rbind(x3, x3.2)
        x7.0_1 <- paste(x7.0[, 1], x7.0[, 4], sep = "@")
        x7.0_1 <- do.call(rbind, strsplit(x7.0_1, "_", fixed = TRUE))
        x7.0 <- data.frame(x7.0_1[, 1], x7.0)
        x7.0 <- data.frame(x7.0[, 8], x7.0)
        x7.0 <- data.frame(x7.0[, c(1, 2, 6, 7, 8)])
    })
    x7.0 <- do.call(rbind, lapply(chr, function(x) x$x7.0))
    x7.2 <- data.frame(paste(r$gene[, 1], r$gene[, 2], sep = "@"), 
        1:nrow(r$gene))
    x7.3 <- sort(unique(x7.0[, 3]))
    x7.3 <- data.frame(x7.3, 1:length(x7.3))
    x7.1 <- data.frame(x7.0, x2[match(x7.0[, 1], x2[, 1]), -1])[, 
        -1]
    x7.4 <- data.frame(x7.1, x7.2[match(x7.1[, 1], x7.2[, 1]), 
        2])[, -1]
    x7 <- data.frame(x7.4, x7.3[match(x7.4[, 1], x7.3[, 1]), 
        -1])[, -1]
    colnames(x7) <- 1:ncol(x7)
    rownames(x7) <- 1:nrow(x7)
    x7.5 <- do.call(rbind, lapply(chr, function(x) x$x2.2))
    x8.1 <- x2[, 1]
    x7[, 3] <- factor(x2[x7[, 3], 1], levels = x2[, 1])
    colnames(x7) <- c("type", "d", "sample", "gene", "seg")
    class(x7) <- c("mcf_data", class(x7))
    ref(x7) <- r
    x7
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
