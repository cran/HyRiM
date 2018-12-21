\name{moment}
\alias{moment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute moments of loss distributions}
\description{
the moment of given order k is computed by numeric integration or summation (in case of discrete distributions)}
\usage{
moment(ld, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ld}{
the loss distribution as obtained from \code{lossDistribution} or \code{mgss}.
}
  \item{k}{
the order of the moment (must be an integer \eqn{\geq 1})
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
the k-th order moment of the given loss distribution
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Stefan Rass
}
\note{
In case of continuous distributions, the value returned is an approximation and based on the internal kernel density approximation. 

For categorical distributions, the function works on the internal probability mass function (which may be different from the empirical distribution in case that the loss distribution has been smoothed during its construction; see \code{\link{lossDistribution}}).

In its current version, \code{cdf} does not vectorize, i.e., cannot be applied to vector arguments \code{x}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
the methods \code{\link{mean}} and \code{\link{variance}} are based on this function. %% ~~objects to See Also as \code{\link{help}}, ~~~
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
