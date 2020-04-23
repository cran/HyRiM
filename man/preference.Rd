\name{preference}
\alias{preference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decision on preferences between loss distributions
}
\description{
This function implements the total ordering on losses, based on treating the moment sequences as hyperreal numbers, and returns the lesser of the loss distribution representatives in the hyperreal space.
}
\usage{
preference(x, y, verbose = FALSE, weights, points = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a loss, being either a number,a distribution or list of distributions (objects of class \code{mosg.lossdistribution})}
  \item{y}{
  a loss, being either a number,a distribution or list of distributions (objects of class \code{mosg.lossdistribution})}
  \item{weights}{
  a vector of \code{n = length(x) = length(y)} nonzero numbers (not necessarily summing up to 1), used only if \code{x} and \code{y} are lists of \code{mosg.lossdistribution} objects corresponding to \eqn{n>1} goals. In that case,
  the i-th goal gets assigend the weight (priority) \code{weights[[i]]}. Defaults to all goals having equal priority if the parameter is missing (\code{weights = rep(1/length(x), length(x))}).
}
  \item{verbose}{
if set to \code{TRUE}, the function returns the preferred of its arguments directly (thus, giving back \code{x} or \code{y}). If set to \code{FALSE} (default), then it returns the argument index (1 = \code{x}, 2 = \code{y}) or 0 in case that \eqn{x=y}.
}
  \item{points}{
the number of points at which the distributions are evaluated numerically to determine the preference.
}
}
\details{
Deciding the preference ordering defined in terms of moment sequence as proposed in (Rass, 2015). To avoid having to compute all moments up to an unknown order, this function decides by looking at the tails of the distribution, returning the one with faster decaying tail as the preferred distribution. This method delivers exact decisions for discrete distributions, but is only an approximate approach for continous densities.
}
\value{
the result is either a copy of the input parameter \code{x} or \code{y}, depending on which distribution is preferred. 
}
\references{
S. Rass. On Game-Theoretic Risk Management (Part One). Towards a Theory of
Games with Payoffs that are Probability-Distributions. ArXiv e-prints, June 2015.
http://arxiv.org/abs/1506.07368.
}
\author{
Stefan Rass
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lossDistribution}}, \code{lossDistribution.mosg}, \code{\link{print.mosg.lossdistribution}}
}
\examples{
# use data from CVSS risk assessments
cvss1base <- c(10,6.4,9,7.9,7.1,9)
cvss2base <- c(10,7.9,8.2,7.4,10,8.5,9,9,8.7)
ld1 <- lossDistribution(cvss1base)
ld2 <- lossDistribution(cvss2base)
lowerRisk <- preference(ld1, ld2) # get the result for later use
preference(ld1, ld2, verbose=TRUE) # view the detailed answer
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
