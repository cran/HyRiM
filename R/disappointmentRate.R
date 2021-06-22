# The all-in-one include header for the HyRiM R package
#
# Authors:         Sandra König, sandra.koenig@ait.ac.at
#                  Stefan Rass, stefan.rass@aau.at
#
# Copyright (C) 2014-2020 AIT Austrian Institute of Technology
# AIT Austrian Institute of Technology GmbH
# Giefinggasse 4 | 1210 Vienna | Austria
# http://www.ait.ac.at
#
# This file is part of the AIT HyRiM R Package.
# The AIT HyRiM R Package can be used for non-commercial and
# academic as well as evaluation purposes. For further information on
# commercial use, please contact the authors!
#
# The AIT HyRiM R Package is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# The AIT HyRiM R Package is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the AIT HyRiM R Package.
# If not, see <http://www.gnu.org/licenses/>.
#
disappointmentRate <-
function(d, x, y, verbose = TRUE, ...) {
  if (class(d) == "mosg.lossdistribution") {
    expectation <- moment(ld = d, k = 1)
    return(1 - cdf(ld = d, x = expectation));
  }
  else if (class(d) == "matrix") {
    n <- nrow(d)
    m <- ncol(d)
    # check if an equilibrium was provided,
    # and if not, compute one
    # cases 1+2: only one of the equilibria is missing, so we compute
    # the best response to what the other player did. This is necessarily
    # a pure strategy equilibrium
    if (missing(x) && !missing(y)) {
      if (verbose) { message("looking for equilibrium for given 2nd player strategy") }
      x <- rep(0, times = n)
      x[which.min(d %*% y)] <- 1
    }
    if (missing(y) && !missing(x)) {
      if (verbose) { message("looking for equilibrium for given 1st player strategy") }
      y <- rep(0, times = m)
      y[which.min(t(x) %*% d)] <- 1
    }
    # case 3: no equilibrium is given, so find ourselves one
    if (missing(x) && missing(y)) {
      if (verbose) { message("looking for equilibrium for given 1st player strategy") }
      auxG <- mosg(n, m, goals = 1, losses = as.vector(d), byrow = FALSE)
      eq <- mgss(auxG, ...)
      x <- as.vector(eq$optimalDefense)
      y <- as.vector(eq$optimalAttacks[,1])
    }
    if (verbose) {
      message("using the following equilibrium: ")
      message("x* = ")
      print(x)
      message("y* = ")
      print(y)
    }

    # now, replace all payoffs by their indicators for a minimizing
    # first player
    v <- (x%*%d%*%y)[1]
    A <- d
    A[A < v] <- 1
    A[A >= v] <- 0
    return((t(x) %*% A %*% y)[1])
  }
  stop("disappointment rate only defined for classes 'mosg.lossdistribution' and 'matrix'")
}
