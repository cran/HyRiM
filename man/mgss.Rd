\name{mgss}
\alias{mgss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute a multi-goal security strategy}
\description{
The function applies fictitious play in a one-against all competition as described by (Sela, 1999) to approximate an equilibrium in the auxiliary game belonging to the input \code{mosg}, using the method proposed by (Rass and Rainer, 2014).
}
\usage{
mgss(G, weights, cutOff, ord = 5, eps, T = 1000, points = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
a multi-objective game constructed using \code{mosg}
}
  \item{weights}{
each goal in \code{G} can be assigned a weight to reflect its priority. If missing, the weights default to be all equal. The weights do not need to sum up to 1 (and are normalized towards a unit sum otherweise), but need to be all non-negative.
}
  \item{cutOff}{
(only used for continuous loss distributions) the maximal loss for which no events are expected or otherwise the risk of exceeding \code{cutOff} are accepted. If missing, this value defaults to the maximal observation on which the loss distributions were constructed (equivalently, the right end of their common support).
}
  \item{ord}{
the order up to which the loss distribution shall be approximated. This value may be set to high orders when it is necessary to distinguish distributions that are similar at the tails.
}
  \item{eps}{
the accuracy at which the fictitious play iteration is being stopped. The error is herein determined as the maximum-norm of the difference between approximations of the equilibrium distributions obtained in the current and the previous iteration of FP.
}
  \item{T}{
the maximal number of iterations after which FP is stopped (returning whatever approximation has been obtained).
}
  \item{points}{
the number of points at which the resulting equilibrium loss distributions are approximated.
}
}
\details{
The function assumes the payoffs to be Gaussian kernel density estimators (constructed using the function \code{lossDistribution}), and computes a Taylor-polynomial approximation at the \eqn{x} equal to \code{cutOff} for each distribution up to order \code{ord}. Preferences are decided using the method put forth in (Rass, 2015), using sign-alternating derivatives. 
}
\value{
An object of class \code{mosg.equilibrium}, containing the following fields:
\item{optimalDefense}{a discrete probability distribution over the action space of player 1 (defender)}
\item{optimalAttacks}{a discrete probability distribution over the action space of player 2 (attacker)}
\item{assurances}{a list of loss distributions valid under the assumption that player 1 adheres to the \code{optimalDefense} distribution in its randomized action choices. The list can be accessed by the names for each goal as specified through the input \code{mosg} object \code{G}. Each distribution within \code{assurances} is a mixed loss distribution constructed using \code{lossDistribution}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
S. Rass and B. Rainer. Numerical computation of multi-goal security strategies.
In Radha Poovendran and Walid Saad, editors, Decision and Game Theory for
Security, LNCS 8840, pages 118-133. Springer, 2014.

Aner Sela. Fictitious play in 'one-against-all' multi-player games. Economic Theory,
14:635-651, 1999.

S. Rass. On Game-Theoretic Risk Management (Part One). Towards a Theory of
Games with Payoffs that are Probability-Distributions. ArXiv e-prints, June 2015.
http://arxiv.org/abs/1506.07368.

D. Lozovanu, D. Solomon, and A. Zelikovsky. Multiobjective games and determining
pareto-nash equilibria. Buletinul Academiei de Stiinte a Republicii Moldova
Matematica, 3(49):115-122, 2005. ISSN 1024-7696.

% ~put references to the literature/web site here ~
}

\author{
Sandra Koenig, Stefan Rass
}
\note{
If both, the accuracy parameter \code{eps} and the maximal number of iterations \code{T} are specified, then FP stops upon whichever condition becomes satisfied first.

Note that the output loss distributions (accessible by the list \code{assurances}) cannot be used to construct a subsequent game (see \code{mosg}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
A brief info on the results can be obtained by \code{\link{print.mosg.equilibrium}}, and a more detailed summary (showing all loss distributions in detail) is obtained by \code{\link{summary.mosg.equilibrium}}.
}
\examples{
library(compare)
library(orthopolynom)
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
eq <- mgss(G,T=1000,weights=c(0.25,0.5,0.25))
print(eq)
summary(eq)
# construct another loss distribution from a given behavior in the game G
suboptimal <- lossDistribution.mosg(G, c(0.1,0.1,0.8), c(0.2,0.3,0.5))
plot(suboptimal)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
