\name{mgss}
\alias{mgss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute a multi-goal security strategy}
\description{
Finds security strategy that assures a maximal loss w.r.t. all goals of the given game, delivering a Pareto-efficient loss bound. Internally, it constructs an auxiliary one-against-all game and uses a sequence of linear programs to compute a Pareto-Nash equilibrium therein (Rass, 2013), using the methods described by (Lozovanu et al 2005; Rass 2015).
}
\usage{
mgss(G, weights, cutOff, ord = 5, fbr = FALSE, points = 512)
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
the order up to which a continuous loss distribution shall be approximated. This value may be set to high orders when it is necessary to distinguish distributions that are similar at the tails.
}
\item{fbr}{
if set to \code{TRUE}, instruct the function to additionally compute the  best replies regarding each goal individually, assuming that defender plays \code{optimalDefense} as a leader, and the attacker per goal follows (follower's best reply). These replies are always pure strategies.
}
  \item{points}{
the number of points at which the resulting equilibrium loss distributions are evaluated numerically.
}
}
\details{
For continuous loss distributions, the function uses a Gaussian kernel density approximation (constructed using the function \code{lossDistribution}), and computes a Taylor-polynomial approximation at the \eqn{x} equal to \code{cutOff} for each distribution up to order \code{ord}. Preferences are decided using the method put forth in (Rass, König, Schauer, 2016), using sign-alternating derivatives, representing a distribution by a vector with \code{ord} elements. Categorical distributions are represented likewise directly by the vector of their probability masses. In both cases, Theorem 3 and Lemma 4 in (Rass, König, Schauer 2016) allow a decision about the stochastic order between two distributions by a lexicographic comparison of the order between the vector-representation. Thus, the computed optima are lex-order optimal. Constructing a game using \code{mosg} with vectors in the payoff description can, consequently, allows to use \code{mgss} to compute lex-order optimal equilibria for multi-criteria games.
}
\value{
An object of class \code{mosg.equilibrium}, containing the following fields:
\item{optimalDefense}{a discrete probability distribution over the action space of player 1 (defender)}
\item{optimalAttacks}{a discrete probability distribution over the action space of player 2 (attacker). Note that this is \emph{not} a best-response to the player 1's \code{optimalDefense}, but rather the best that the attacker could
do if the game were \emph{just about the particular goal} that the attacker refers to. This worst-case scenario assumes that the defender would focus all its efforts to that single goal.}
\item{assurances}{a list of loss distributions valid under the assumption that player 1 adheres to the \code{optimalDefense} distribution in its randomized action choices, while the opponent plays its own zero-sum equilibrium strategy in the game that is only (and exclusively) about this particular goal. This value has to be interpreted with care, as it assumes that player 1 would put all efforts into a defense for the particular goal, but in reality, will have multiple criteria to simultaneously optimize.
This means that the attacker, in turn, could adapt to the \code{optimalDefense} of player 1, to cause more damage. The given assurance is thus only an upper bound of the worst-possible damage, under the assumption that player 1 would focus only on this particular goal.

The list can be accessed by the names for each goal as specified through the input \code{mosg} object \code{G}. Each distribution within \code{assurances} is a mixed loss distribution constructed using \code{lossDistribution}}
\item{br_to_optimalDefense}{This is a vector of best replies per goal for a leading defender playing the fixed strategy \code{optimalDefense}, and letting the adversary (player 2) follow. It is the (stochastically largest) damage among \eqn{optimalDefense^T\cdot A_p}, when \eqn{A_p} is the game structure for the \eqn{p}-th goal; the vector \code{br_to_optimalDefense} contains the indices of the individually best replies, pointing into the list of attack strategies.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
S. Rass, S. König, S. Schauer. Decisions with Uncertain Consequences-A Total Ordering on Loss-Distributions. PLoS ONE 11, e0168583. 2016, https://doi.org/10.1371/journal.pone.0168583

S. Rass. On Game-Theoretic Risk Management (Part One). Towards a Theory of
Games with Payoffs that are Probability-Distributions. June 2015.
http://arxiv.org/abs/1506.07368.

S. Rass. On Game-Theoretic Risk Management (Part Two). Algorithms to Compute Nash-Equilibria in Games with Distributions as Payoffs, 2015, arXiv:1511.08591v1 [q-fin.EC].

D. Lozovanu, D. Solomon, and A. Zelikovsky. Multiobjective games and determining
pareto-nash equilibria. Buletinul Academiei de Stiinte a Republicii Moldova
Matematica, 3(49):115-122, 2005. ISSN 1024-7696.

% ~put references to the literature/web site here ~
}

\author{
Sandra Koenig, Stefan Rass
}
\note{
The output loss distributions (accessible by the list \code{assurances}) cannot be used to construct a subsequent game (see \code{mosg}), since continuous distributions are represented as a sequence of points, rather than raw data or probability masses.

As of version 2.0.0 of the package, this function is no longer downwards compatible to earlier versions of itself, since the method of computation (formerly fictitious play) was replaced by linear programming to give exact solutions rather than approximations. Consequently, the parameters \code{T} (iteration count) and \code{eps} (accuracy) have become useless and have been removed after version 1.0.4.
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
eq <- mgss(G,weights=c(0.25,0.5,0.25))
print(eq)
summary(eq)

# construct another loss distribution from a given behavior in the game G
suboptimal <- lossDistribution.mosg(G, c(0.1,0.1,0.8), c(0.2,0.3,0.5))
plot(suboptimal)

# compute an equilibrium in a standard matrix game
#     [,1] [,2]
#[1,]    3    4
#[2,]    6    1
G <- mosg(n = 2, m = 2, goals = 1,
          losses = list(3,6,4,1), byrow=FALSE,
          attacksDescr = c("a1", "a2"))
mgss(G, fbr=TRUE)  # compute an equilibrium, including best replies if the adversary is a follower

# get best replies if there would be a following
# adversary per goal (taking the defender as a leader)
G$attacksDescriptions[eq$br_to_optimalDefense]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
