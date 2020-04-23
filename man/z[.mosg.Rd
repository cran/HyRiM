\name{[.mosg}
\alias{[.mosg}
\alias{[<-.mosg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract or replace parts of a game's payoff matrix
}
\description{
Construct a new game by taking out a specified set of rows, columns and goals from a given game \eqn{G}. The new game inherits all descriptions (rows, cols and goals) from the \eqn{G}, and has its list of loss distributions organized in the same way (by rows or columns) as \eqn{G}.

The extraction or substitution works like as for data frames (see \code{\link{[.data.frame}}).
Strategies for both players, as well as goals, can equivalently be addressed by their string-names.
}
\usage{
\method{[}{mosg}(x,i,j,k=NULL)

\method{[}{mosg}(x,i,j,k=NULL) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a game of class \link{mosg}
}
  \item{i,j,k}{
a numeric value or numeric vector of row incides \code{i}, colum indices \code{j}, or goals \code{k}.
}
  \item{value}{
a \link{list} of \link{lossDistribution} objects, or a game object of class \code{\link{mosg}}.
}
}
\details{
For \code{[} extraction of elements from a payoff matrix, omitting any index dimension selects all elements in the respective dimension. Supplying negative values excludes the respective elements. For example, \code{G[c(1:3),1]} returns a game with only the rows 1..3 of \code{G}, but all column strategies that \code{G} had, and only the first of \code{G}'s goals retained.

For \code{[<-}, the list of substitute values needs to be of the same length as the number of elements addressed by the triple \code{(i,j,k)}, otherwise an error is returned. If the new elements come from another game object, say \code{G2}, only the loss distributions get replaced, but not the names of the strategies. The replacement checks if \code{G2} has its list of loss distributions organized in the same way as \code{G}, i.e., row-by-row or column-by-column. If there is a mismatch, the substitution is nonetheless done, but a warning about this issue is printed.
}
\value{
\code{[} returns a freshly constructed game object.
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

%% ~Make other sections like Warning with
\section{Warning}{For \code{[<-}, be aware that the replacement \emph{does not} semantically check if the newly incoming loss distributions make sense as elements of the new game (e.g., they can have different supports, or be discrete/continuous while the game was continuous/discrete in its payoffs). Respective errors may only subsequently come up when the modified or extracted game is used.}

\seealso{
\code{\link{[.data.frame}} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## raw data (PURELY ARTIFICIAL, for demo purposes only)
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
            m=3,
            payoffs,
            goals=2,
            goalDescriptions=c("g1", "g2"),
            defensesDescr = c("d1", "d2"),
            attacksDescr = c("a1", "a2", "a3"))

# modify the game by subsetting
G[,c(1,2),] # select only the first two strategies
G[,-3,] # exclude the third strategy (equivalent to before)

# replace a 2x2 subgame related to the second goal
# (replacement data is chosen arbitrarily here)
G2 <- mosg(n=2, m=2, goals=1, losses = list(f111,f112,f121, f122))
G[,c(1,2),1] <- G2  # replace the subgame

# construct another replacement game that is organized different (by column)
G2 <- mosg(n=2, m=2, goals=1, losses = list(f111,f112,f121, f122), byrow=FALSE)
G[,c(1,2),1] <- G2  # this will issue a warning

# plot a submatrix from the game
plot(G[-2,c(1,2),], goal=2)

}
