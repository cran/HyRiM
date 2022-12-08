\name{mosg}
\alias{mosg}
\alias{print.mosg}
\alias{plot.mosg}
\alias{lossDistribution.mosg}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construction and handling of multi-objective security games
}
\description{
this function takes a list of loss distributions construced using \code{lossDistribution}, along with a specification of the game's shape (number of strategies for both players and number of goals for the first player), and returns an object suitable for analysis by \code{mgss} to compute a multi-goal security strategy. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mosg( n,
      m,
      goals,
      losses,
      byrow = TRUE,
      goalDescriptions = NULL,
      defensesDescr = NULL,
      attacksDescr = NULL)

\method{print}{mosg}(x, ...)

\method{plot}{mosg}(x,
     goal = 1,
     points = 100,
     cutoff = NULL,
     largeGame = FALSE,
     subPlotWidth = 2,
     subPlotHeight = 2,
     cleanUp = TRUE, ...)

# construct a loss distribution by playing a given strategy in the game G
\method{lossDistribution}{mosg}(G, player1Strat, player2Strat, points = 512, goal = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of defense strategies (cardinality of the action space for player 1)
}
  \item{m}{
number of attack stratgies (cardinality of the action space for player 2)
}
  \item{goals}{
number of goals for player 1 (must be \eqn{\geq 1})
}
  \item{losses}{
a \code{list} with \code{n*m*goals} entries, which specifies a total of \code{goals} game matrices, each with shape \code{n}-by-\code{m}. The way in which the game matrices are filled from this list is controlled by the parameter \code{byrow}. Note that in every case, it is assumed that one matrix is specified after the other in the list.

Furthermore, the function assumes all loss distributions having a common support. This is only explicitly verified for discrete distributions (with errors reported), but implicitly assumed to hold for continuous distributions without further checks.

Typically, a game will be constructed from a list of loss distributions obtained by invocations of \code{lossDistribution}.

Games can be defined with real-valued (scalar) payoffs if a list of numbers is provided instead. Internally, the function converts these numbers into Bernoulli distributions; a scalar payoff \eqn{a} is converted into a Bernoulli random variable \eqn{X} having \eqn{\Pr(X=a)=p\propto a}. This conversion is equivalent to an invocation of \code{lossDistribution} with the parameters \code{dat=c(1-p, p)}, \code{discrete=TRUE}, \code{dataType="pdf"}, \code{smoothing="none"}, \code{bw = 1} and \code{supp=c(1,2)}.

If the list of losses comes as a list of vectors, \code{mosg} will construct a game assuming a lexicographic order on the loss vectors
(with the order being determined from left to right along the coordinates). To this end, \code{mosg} checks for all loss vectors to have the same length (otherwise, an error is reported). Negative and zero values in the loss vector \emph{are allowed}.
}
  \item{byrow}{
by default (\code{TRUE}), the game matrices are filled row-by-row from list \code{losses}. If set to \code{FALSE}, then the game matrices are filled column-by-column.
}
  \item{goalDescriptions}{
if specified, this can be any vector (e.g., textual descriptions) for the goals. Defaults to 1, 2, 3, ... if missing. The length must be equal to \code{goals}.
}
  \item{defensesDescr}{
if specified, this can be any vector (e.g., textual descriptions) for the defense strategies. Defaults to 1, 2, 3, ... if missing. The length must be equal to \code{n}.}
  \item{attacksDescr}{
if specified, this can be any vector (e.g., textual descriptions) for the attack strategies. Defaults to 1, 2, 3, ... if missing. The length must be equal to \code{m}.}

for the functions \code{print}, \code{summary} and \code{plot}

\item{x}{a game, object of class "mosg", as constructed by the function \code{mosg}}

The function \code{plot} additionally takes the following parameters:
\item{goal}{an integer referring to the goal of interest (for plotting or to construct a loss distribution for). Defaults to the first goal if omitted.}
\item{points}{The number of points at which the density is evaluated (for continuous losses); this parameter is ignored for categorical losses.}
\item{cutoff}{the cutoff point at which all densities shall be truncated before plotting (note that the mass functions are rescaled towards unit mass).}
The plot function overrides the following settings internally (so supplying these as parameters will raise an error): \code{xlab, ylab, main, type, names.arg} and \code{font.main} (applying differently for bar and line plots)

\item{largeGame}{if the plot exits with the error "figure margins too large", one can set this
parameter to TRUE, causing plot to write to a temporary SVG file (scalable vector graphics), to avoid the figure space issue and hence
the error. The price is a (potentially much) slower plotting, since the system creates the file, and loads it afterwards from
the harddisk (cleaning up the file after displaying it). The size of the plot is controllable by setting
the parameters \code{subPlotWidth} and \code{subPlotHeight}, see below.}

\item{subPlotWidth}{the width in inches for each payoff distribution in the game matrix. This parameter is
ignored when \code{largeGame} is set to FALSE (the default).}

\item{subPlotHeight}{the height in inches for each payoff distribution in the game matrix. This parameter is
ignored when \code{largeGame} is set to FALSE (the default).}

\item{cleanUp}{If the graph is to be used in other programs, one can supply \code{cleanUp = FALSE} to retain the temporary SVG file for subsequent use and prints a message where to find the file. By default, the temporary file gets deleted.}

The function \code{lossDistribution.mosg} can be used to play any (given) strategies for player 1 and player 2, and compute the resulting loss from the game.

\item{G}{a game constructed by \code{mosg} to deliver the loss distribution through its game matrices.}
\item{player1Strat}{a discrete distribution over the action space for the defending player 1 in the game \code{G}}
\item{player2Strat}{a discrete distribution over the action space for the attacking player 2 in the game \code{G}}
\item{...}{further arguments passed to or from other methods}
}
\details{
Upon input, the function does some consistency checks, such as testing the length of the parameter \code{losses} to be equal to n*m*goals. The loss distributions are checked for mutual consistency in terms of all being continuous or all being discrete (a mix is not allowed), and all being not mixed distributions (that is, the output distribution of a previous call to \code{mgss} cannot be used as input to this function).

The functions \code{print.mosg} gives a brief overview of the game, listing only the shape and strategies for both players. For detailed information, use \code{summary} on a specific loss distribution in the list for the game (field \code{losses}).

For plotting games, \code{plot.mosg} constructs an (n x m)-matrix of loss distributions with rows and columns in the grid being labeled by the values in \code{defensesDescr} and \code{attacksDescr}. The plot heading is the name for the specified goal. The function makes no changes to the plot parameters, so fine tuning can be done by changing the settings using the \link{par} function.

The function \code{lossDistribution.mosg} can be used to compute the distribution \eqn{x^T*A*y}, for the payoff distribution matrix \eqn{A}, and mixed strategies \eqn{x} (\code{player1strat}) and \eqn{y} (\code{player2strat}) in the game. The computation is by a pointwise addition of loss distributions, with the number of points being specifiable by the parameter \code{points}, which defaults to 512.
}
\value{
The function returns an object of class \code{mosg}, usable with the function \code{mgss} to determine a security strategy (i.e., a lexicographic Nash equilibrium assuming a zero-sum one-against-all competition). The fields returned in the \code{mosg} object are filled with the input values supplied. In detail, the fields are:
\item{nDefenses}{the value of the parameter \code{n}}
\item{nAttacks}{the value of the parameter \code{m}}
\item{dim}{the value of the parameter \code{goals}}
\item{attacksDescriptions, defensesDescriptions, goalDescriptions}{if supplied, then these are filled with the values of \code{goalDescr}, \code{defensesDescr} and \code{attacksDescr}; otherwise, they contain the default values described above.}
\item{maximumLoss}{the maximal loss taken over all specified loss distributions}
\item{loc}{a locus-function for accessing the list \code{losses} using a triple notation (goal,i,j), where goal addresses the game matrix and i,j are the row and column indices (starting from 1 as the smallest index). This function is used internally (only).}
}
%\references{
%% ~put references to the literature/web site here ~
%}

\note{It is important to remark that player 1 is always minimizing. To treat a maximizing player, one must reconstruct the game using regrets instead of losses, i.e., if the data for a specific loss distribution is \eqn{D}, then the game for a maximizing player 1 must be constructed from \code{(max(D) - D)} instead of \code{D}.}

\author{
Stefan Rass
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning }{
Games constructed with real-valued payoffs or payoff vectors over the reals are allowed with negative or zero values in the list of losses. In that case, embeds the loss values or vector into a \link{lossDistribution} object after shifting and scaling the values into the strictly positive range. This operation creates a strategically equivalent game, i.e., leaves the set of equilibria unchanged, yet the resulting \link{mosg} object \emph{is not} useful with the \link{lossDistribution.mosg}, \link{moment}, \link{cdf}, or any other member functions for \link{lossDistribution} objects obtained from equilibria. Those have to be computed manually. Be aware that there will be no warnings issued whatsoever in that case of misuse, since the \link{lossDistribution} objects constructed to carry the real or vector-valued payoffs of the original game carry no information about the semantics of the values or vectors that they have been created from. Hence, the \emph{computation of equilibria works correctly} using \link{mosg}, while any further analysis (including plots)  \emph{needs to be done manually}.
}
\seealso{
Security strategies for a \code{\link{mosg}} object can be obtained by calling \code{\link{mgss}}. The game itself can be constructed from the output of \code{\link{lossDistribution}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(compare)

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
print(G)
summary(G)
plot(G)

# construct and solve scalar valued (classical) game;
# losses are all numbers (degenerate distributions)
# the resulting matrix game has the payoff structure:
#     [,1] [,2]
#[1,]    3    4
#[2,]    6    1
G <- mosg(n = 2, m = 2, goals = 1, losses = list(3,6,4,1), byrow=FALSE)
mgss(G)  # compute a lexicographic Nash equilibrium
}

