\name{variance}
\alias{variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the approximate variance of a loss distribution.
}
\description{
The computation is based on Steiner's theorem \eqn{\textrm{var}(X) = \textrm{E}(X^2) - (\textrm{E}(X))^2}, where the
respective first and second moments are computed using the \code{moment} function (from this package). Internally, these 
functions operate on the approximate kernel density estimation for both, continuous and categorical distributions
(see the \code{lossDistribution} function for details).
}
\usage{
variance(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{mosg.lossDistribution}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
the approximate variance value
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
the function works on the internal probability mass function (which may be different from the empirical distribution in case that the loss distribution has been smoothed during its construction; see \code{\link{lossDistribution}}). The function delivers only an approximate variance, whose error is due to numeric roundoff errors (known to occur in Steiner's formula), and the fact that the computation is done on an approximate density (rather than the 
empirical distribution).
%since the function acts on the loss distribution being a mere approximation of the %underlying data, 
%the returned value may be quite different from the empirical variance computed from the %data 
%directly (as shown in the example below)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{moment}}, \code{\link{lossDistribution}}
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (x) 
%{
%    moment(x, 2) - moment(x, 1)^2
%  }
%}

x <- c(10,6.4,9,7.9,7.1,9)
ld <- lossDistribution(x)
variance(ld)
var(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
