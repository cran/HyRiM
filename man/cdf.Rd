\name{cdf}
\alias{cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
(cumulative) loss distribution function}
\description{
returns the numeric values of the cumulative loss distribution \code{ld} evaluated at \code{x}, i.e., \eqn{\textrm{Pr}(X \leq x)}, where \eqn{X\sim \textrm{ld}}.
}
\usage{
cdf(ld, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ld}{
the loss distribution as obtained from \code{lossDistribution} or \code{mgss}.
}
  \item{x}{
the point at which the distribution function shall be evaluated (must be a numeric; vectors are not supported yet)
}
}
\details{
the function internally distinguishes discrete and continous disributions only in terms of rounding its argument to the largest integer less than x. Its value is obtained by numeric integration of the internal representation of the loss distribution (in the continuous case).

For discrete distributions, the function works on the internal probability mass function (which may be different from the empirical distribution in case that the loss distribution has been smoothed during its construction; see \code{\link{lossDistribution}}).
}
\value{
an approximation for the probability Pr(ld<=x).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Stefan Rass
}
\note{
in its current version, \code{cdf} does not vectorize, i.e., cannot be applied to vector arguments \code{x}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
suitable inputs for this function are provided by \code{\link{lossDistribution}} and \code{\link{mgss}}.%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cvss1base <- c(10,6.4,9,7.9,7.1,9)
ld <- lossDistribution(cvss1base)
cdf(ld, 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
