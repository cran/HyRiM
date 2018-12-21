\name{lossDistribution}
\alias{lossDistribution}
\alias{print.mosg.lossdistribution}
\alias{summary.mosg.lossdistribution}
\alias{plot.mosg.lossdistribution}
\alias{quantile.mosg.lossdistribution}
\alias{mean.mosg.lossdistribution}
\alias{density.mosg.lossdistribution}
\alias{print.summary.mosg.lossdistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
construction and handling of loss distributions
}
\description{
Loss distributions can be constructed from both, continuous and categorical data. In any case, the input data must be a list (vector) of at least two numeric values all being \eqn{\geq 1}. For discrete data, the function additionally takes the full range of categories, all being represented as integers (with the lowest category having the number 1). 
}
\usage{
# construct a loss distribution from data
lossDistribution(
  dat, 
  discrete = FALSE, 
  dataType = c("raw", "pdf", "cdf"),
  supp = NULL, 
  smoothing = c("none", "ongaps", "always"), 
  bw = NULL)
# get information about the loss distribution
\method{print}{mosg.lossdistribution}(x, ...) 
\method{summary}{mosg.lossdistribution}(object, ...)
\method{print}{summary.mosg.lossdistribution}(x, ...)
\method{plot}{mosg.lossdistribution}(x, points = 100, xlab = "", ylab = "",
           main = "", p = 0.999, newPlot = TRUE, cutoff = NULL, ...)  
# get quantiative information about the distribution 
\method{quantile}{mosg.lossdistribution}(x, p, eps = 0.001, ...)
\method{mean}{mosg.lossdistribution}(x, ...)
# evaluate the loss density function
\method{density}{mosg.lossdistribution}(x, t, ...)  
# for the cumulative distribution function, see the function 'cdf'
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a vector of at least two input observations (all \eqn{\geq 1} required)
}
  \item{discrete}{
defaults to \code{FALSE}. If set to \code{TRUE}, the loss distribution is constructed as discrete. In that case, a value for \code{supp} is required.
}
\item{dataType}{
applies only if \code{discrete=TRUE}, and specifies how the values in \code{dat} are to be interpreted. Defaults to \code{raw}, by which the data is taken as observations. Given as \code{pdf}, the values in \code{dat} are directly interpreted as a probability density (checked for nonnegativity and re-normalized if necessary). If the data type is specified as \code{cdf}, then the values in \code{dat} are taken as cumulative distribution function, i.e., checked to be non-decreasing, non-negative and re-normalized to 1 if necessary. 
}
  \item{supp}{
if the parameter \code{discrete} is set to \code{TRUE}, then this parameter must be set as a vector of two elements, specifying the minimal and maximal category, e.g. \code{supp=c(1,5)}. 
}
\item{bw}{
the bandwidth parameter (numeric value) for kernel smoothing. Defaults internally to the result of \link{bw.nrd0} if omitted.
}
\item{x}{
a loss distribution object returned by \code{lossDistribution} or \code{mgss}, or a value within the support of a loss distribution.
}
\item{t}{
a value within the support of \code{ld} or a summary object for a loss distribution.
}
\item{object}{
a loss distribution object
}
  \item{eps}{
the accuracy at which the quantile is approximated (see the details below).
}
\item{smoothing}{
  string; partially matched with "none" (default), "ongaps", and "always". If set to "always", then the function computes a discrete kernel density estimate (using a discretized version of a Gaussian density with a bandwidth as computed by \code{bw.nrd0} (Silverman's rule)), to assign categories with zero probability a positive likelihood. If set to "ongaps", then the smoothing is applied only if necessary (i.e., if the probability mass is zero on at least one category).
}

the function \code{plot.mosg.lossdistribution} takes the parameters:
\item{points}{the number of points at which loss densities are is evaluated (numerically) for plotting.}
\item{xlab}{a label for the x-axis in the plot.}
\item{ylab}{a label for the y-axis in the plot.}
\item{main}{a title for the plot}
\item{p}{a quantile that determines the plot range for the loss distribution}
\item{newPlot}{if set to \code{TRUE}, then a new plot is opened. Otherwise, the plot is added to the current plot window (typcially used by \code{plot.mosg} to visualize game matrices).}
\item{cutoff}{the cutoff point at which all densities shall be truncated before plotting (note that the mass functions are rescaled towards unit mass).}
\item{...}{further arguments passed to or from other methods}
}
\details{
The function internally computes a Gaussian kernel density estimator (KDE; using Silverman's rule of thumb for the bandwidth selection) on the continuous data. The distribution is truncated at the maximal observation supplied + 5*the bandwidth of the Gaussian KDE, or equivalently, at the right end of the support in case of discrete distributions.

For discrete distributions, missing observations are handled by smoothing the density (by convolution with a discretized Gaussian kernel). As an alternative, a re-definition of categories may be considered.

Degenerate distributions are not supported! The construction of classical games with real-valued payoffs works directly through \link{mosg} by supplying a list of values rather than loss distributions. See the example given with \link{mosg}.

The generic functions \code{quantile}, \code{mean} and \code{density} both distinguish discrete from continuous distributions in the way of how values are being computed. 

Quantiles are computed using the direct definition as an approximation \code{y} so that x = Pr(ld <= y). For continuous distributions, a bisective search is performed to approximate the inverse cumulative distribution function. For discret distributions, \code{quantile} works with cumulative sums. The accuracy parameter \code{eps} passed to \code{quantile} causes the bisective search to stop if the search interval has a length less than \code{eps}. In that case, the middle of the interval is returned. For discrete distributions, the computation is done by cumulative sums on the discrete probability mass function.

\code{mean} either invokes \code{moment(ld, 1)} to compute the first moment.

\code{density} is either a wrapper for the internal representation by the \code{function} object \code{lossdistr}, or directly accesses the probability mass function as internally stored in the field \code{dpdf} (see the 'values' section below).

For visualization, \code{plot} produces a bar plot for categorical distributions (over categories as specified by the \code{supp} field; see the 'values' section below), and for  continous distributions, a continuous line plot is returned on the range 1...max(\code{range} + 5*\code{bw}), where the values are described below. To ease comparison and a visual inspection of the game matrix, the default plot ranges can be overridden by supplying \code{xlim} and \code{ylim} for the plot function. 
}
\value{
The return values of \code{lossDistribution} is an object of class \code{mosg.lossdistribution}. The same goes for \code{lossDistribution.mosg}.
\item{observations}{carries over the data vector supplied to construct the distribution.}
\item{range}{the minimal and maximal loss observed, as a 2-element vector. For loss distributions induced by games, the range is the smallest interval covering the ranges of all distributions in the game.}
\item{bw}{the bandwidth used for the kernel density approximate.}
\item{lossdistr}{a \code{function} embodying the kernel density (probability mass function) as a spline function (for continuous densities only)}
\item{normalizationFactor}{the factor by which \code{lossdistr} must be multiplied (to normalize under the truncation at max(observations) + 5*bw.}
\item{is.mixedDistribution}{a flag indicating whether or not the distribution was constructed by a call to \code{lossDistribution} or the generic function \code{lossDistribution.mosg}.}
\item{is.discrete}{a flag set to \code{TRUE} if the distribution is over categories}
\item{dpdf}{if \code{is.discrete} is \code{TRUE}, then this is a vector of probability masses over the support (field \code{supp}).}
\item{supp}{if \code{is.discrete} is \code{TRUE}, then this is a 2-element vector specifying the minimal and maximal loss category (represented by integers).}

A \code{summary} returns an object of class \code{mosg.equilibrium.summary}, for which the generic \code{print} function can be applied, and which carries the following fields:
\item{range}{the minimal and maximal observation of the underlying data (if available), or the minimal and maximal losses anticipated for this distribution (e.g., in case of discrete distributions the common support).}
\item{mean}{the first moment as computed by \code{mean}.}
\item{variance}{the variance as computed by \code{variance}.}
\item{quantiles}{a 2x5-matrix of quantiles at levels 10\%, 25\%, 50\%, 75\% and 90\%.}
}

\note{If the plotting throws an error concerning too large figure margins, then adjusting the plot parameters using \link{par} may help, since the plot function does not override any of the current plot settings (e.g., issue \code{par(c(0,0,1,1) + 0.1))} before plotting to reduce the spacing close towards zero))

In some cases, plots may require careful customization to look well, so playing arourd with the other settings as offered by \link{par} can be useful.

If the distribution has been smoothed, then \code{mean}, \code{variance}, \code{quantile}, \code{density} and \code{cdf} will refer to the smoothed version of the distribution. In that case, the returned quantities are mere approximations of the analogous values obtained directly from the underlying data.

}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Stefan Rass
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mosg}}, \code{\link{mgss}}, \code{\link{cdf}}, \code{\link{variance}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cvss1base <- c(10,6.4,9,7.9,7.1,9)
ld <- lossDistribution(cvss1base)
summary(ld)
plot(ld)
# for further examples, see the documentation to 'mosg' and 'mosg.equilibrium'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
