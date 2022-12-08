\name{disappointmentRate}
\alias{disappointmentRate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
computation of the disappointment rate
}
\description{
For a minimizing player, the \emph{disappointment rate} is the likelihood for the loss to exceed its expectation (thus disappoint the defender). For any random loss \eqn{X}, it is given by \eqn{Pr(X > E(X))}.
}
\usage{
disappointmentRate(d, x, y, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
a \link{lossDistribution} object or a \link{matrix}; typically the assurance from a previously computed equilibrium (see \code{\link{mgss}}). In that case, all other parameters are ignored.

Alternatively, one can provide a matrix of real values instead, to compute the disappointment rate in the so-specified zero-sum matrix game. In that case, the other parameters are also taken into considertion.
}
\item{x,y}{the mixed strategies under which the disappointment rate shall be computed. Usually, this would be an equilibrium of the (real-valued) matrix game. If only \code{x} or only \code{y} is supplied, the function computes a best response to the given (mixed) strategy. If both are omitted, the function internally computes an equilibrium by a call to \code{\link{mgss}}.
}
\item{verbose}{if set to \code{FALSE}, suppresses all messaging.
}
\item{...}{further parameters internally passed onwards to \code{\link{mgss}} to compute an equilibrium.
}
}
\details{
The disappointment rate can be taken as an auxiliary goal to optimize, though it is not supported for optimization in the current version of the package. Note that it does not make sense to consider this rate as an isolated (single) goal, since the optimal strategy would then be playing towards maximal losses (with explicit aid of the opponent) in order to minimize the mass to the left of the expected loss. However, it is a quantity of interest when the equilibrium has been computed, as it indicates how ``satisfying'' the equilibrium will be upon playing.
}
\value{
the likelihood to overshoot the expectation of the random loss \eqn{X} with distribution \code{d}, i.e., \eqn{Pr(X > E(X))}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
see for example, F. Gul: "A Theory of Disappointment Aversion", Econometrica, vol. 59, no. 3, p. 667, 1991.
}
\author{
Stefan Rass
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mgss}}
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

# get the disappointment rate for the first security goal g1
disappointmentRate(eq$assurances$g1)

############
# construct a game with one goal and related disappointment
payoffs <- list(f111,f112,f121,f122)
# note that from here onwards, the code is "generic", meaning that
# exactly the same procedure would apply to *any* kind of game that
# we want to play with disappointments, as long as the input data comes
# in the variable "payoffs" (as used in the code below)
expectations <- unlist(lapply(payoffs, mean))
disappointmentRates <- unlist(lapply(payoffs, disappointmentRate))
# put the two goals together in a game
gameWithDisappointment <- c(expectations, disappointmentRates)
G <- mosg( n=2,
          m=2,
          losses=gameWithDisappointment,
          goals=2,
          goalDescriptions=c("revenue", "disappointment"),
          defensesDescr = c("d1", "d2"),
          attacksDescr = c("a1", "a2"))
eq <- mgss(G,weights=c(0.1,0.9))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
