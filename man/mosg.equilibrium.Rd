\name{mosg.equilibrium}
\alias{summary.mosg.equilibrium}
\alias{print.mosg.equilibrium}
\alias{plot.mosg.equilibrium}
\alias{print.mosg.equilibrium.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
embodies all information related to a lexicographic Nash equilibrium computed by the function \code{mgss}.
}
\description{
The generic functions \code{print} and \code{summary} provide brief, and detailed information about the lexicographic Nash equilibrium. The generic function \code{plot} can be used to visualize the equilibrium.
}
\usage{
\method{summary}{mosg.equilibrium}(object, ...)
\method{print}{mosg.equilibrium.summary}(x, ...)
\method{print}{mosg.equilibrium}(x, extended=FALSE, ...)
\method{plot}{mosg.equilibrium}(x, points=100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an \code{mgss} object as returned by the function \code{mgss}.
}
\item{object}{
an \code{mgss} object as returned by the function \code{mgss}.
}
for \code{print.mosg.equilibrium}, the following parameter can be supplied:

\item{extended}{if set to \code{TRUE}, then the individual assurances are printed as well.}

for \code{plot.mosg.equilibrium}, the following parameter can be supplied:

\item{points}{the number of points to evaluate the density function over its support for plotting}

\item{...}{further arguments passed to or from other methods.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
the result returned by the function \code{summary} carries the following fields:
\item{optimalDefense}{a discrete probability distribution over the action space for player 1 (the defender).}
\item{optimalAttacks}{a discrete probability distribution over the action space for player 2 (the attacker).}
\item{assurances}{an optimal loss distribution valid under the assumption that the defender plays \code{optimalDefense} as its mixed strategy. This is a list of \code{mosg.lossdistribution} objects, accessible through their assigend names (coming from the underlying game) or by indices.}
The action spaces for both players are defined in first place by the game for which the equilibrium was computed (via \code{mgss} on a game constructed by \code{mosg}).

\code{print} gives a shortened output restricted only to displaying the optimal defense for the defender and attack strategies per goal (as defined by the underlying game).

\code{summary} returns an object of class \code{summary.mosg.lossdistribution}, which has the fields:
"range"       "mean"        "variance"    "quantiles"   "is.discrete"
\item{range}{the minimal and maximal values of the loss (as anticipated by the observations)}
\item{mean}{the first moment as computed by \code{mean}}
\item{variance}{the variance as computed by \code{variance}}
\item{quantiles}{a 2x5-matrix of quantiles at the 10\%,25\%,50\%,75\% and 90\% level}
\item{is.discrete}{a Boolean flag being \code{TRUE} if the loss distribution is over categories}

\code{plot} displays a grid of plots, starting with the optimal defense behavior plotted as a discrete distribution on top of a (m x 2)-matrix of plots. Each line in this grid shows the discrete optimal attack strategy on the right side (as a bar plot), paired with the loss distribution (extracted from \code{x}) caused when the defender plays \code{optimalDefense} and the attacker plays the respective optimal attack strategy.

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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.mosg.equilibrium}}, \code{\link{mgss}}, \code{\link{mosg}}, \code{\link{lossDistribution}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## raw data (PURELY ARTIFICIAL, for demo purposes only)
# N=100 observations in each category
obs111<-c(rep(1,40),rep(3,20),rep(5,10),rep(7,20),rep(9,10));
obs112<-c(rep(1,50),rep(2,10),rep(4,10),rep(6,20),rep(8,10));
obs121<-c(rep(1,20),rep(4,30),rep(6,20),rep(8,10),rep(10,20));
obs122<-c(rep(1,40),rep(2.5,20),rep(5,20),rep(7.5,10),rep(9,10));
obs211<-c(rep(1,30),rep(2,30),rep(5,10),rep(8,10),rep(10,20));
obs212<-c(rep(1,10),rep(2,10),rep(4,20),rep(7,20),rep(10,40));
obs221<-c(rep(1,30),rep(3,30),rep(4,10),rep(7,20),rep(9,10));
obs222<-c(rep(1,10),rep(3,10),rep(5,50),rep(8,20),rep(10,10));
obs311<-c(rep(1,40),rep(2,30),rep(4,10),rep(7,10),rep(9,10));
obs312<-c(rep(1,20),rep(3,20),rep(4,20),rep(7,20),rep(10,20));
obs321<-c(rep(1,10),rep(3,40),rep(4,30),rep(7,10),rep(9,10));
obs322<-c(rep(1,10),rep(4,30),rep(5,30),rep(7,10),rep(10,20));

## compute payoff densities
f111<-lossDistribution(obs111)
f112<-lossDistribution(obs112)
f121<-lossDistribution(obs121)
f122<-lossDistribution(obs122)
f211<-lossDistribution(obs211)
f212<-lossDistribution(obs212)
f221<-lossDistribution(obs221)
f222<-lossDistribution(obs222)
f311<-lossDistribution(obs311)
f312<-lossDistribution(obs312)
f321<-lossDistribution(obs321)
f322<-lossDistribution(obs322)

payoffs<-list(f111,f112,f121, f122,f211,f212,f221,f222, f311,f312,f321,f322)
G <- mosg( n=2,
            m=2,
            payoffs,
            goals=3,
            goalDescriptions=c("g1", "g2", "g3"),
            defensesDescr = c("d1", "d2"),
            attacksDescr = c("a1", "a2"))
eq <- mgss(G,weights=c(0.25,0.5,0.25))
print(eq)
summary(eq)
plot(eq)

# access the loss distributions computed in the game
summary(eq$assurances$g1)
mean(eq$assurance$g1)    # get the average loss in goal "g1"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
