mosg <-
function(n,
           m,
           goals,
           losses,
           byrow = TRUE,
           goalDescriptions = NULL,
           defensesDescr = NULL,
           attacksDescr = NULL) {
## The all-in-one include header for the HyRiM R package
#
# Authors:         Sandra König, sandra.koenig@ait.ac.at
#                  Stefan Rass, stefan.rass@aau.at
#
# Copyright (C) 2014-2017 AIT Austrian Institute of Technology
# AIT Austrian Institute of Technology GmbH
# Donau-City-Strasse 1 | 1220 Vienna | Austria
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
     mosg <- NULL
    if (missing(n) || missing(m) || missing(goals)) {
      stop("shape of the game (n,m,goals) is not completely specified")
    }
    mosg$nDefenses <- n
    mosg$nAttacks <- m
    if (n * m * goals > length(losses)) {
      stop("too few losses given")
    }
    if (byrow) {
      mosg$loc <- function(goal, i, j) {
        (goal - 1) * n * m + (i - 1) * m + j
      }
    }
    else {
      mosg$loc <- function(goal, i, j) {
        (goal - 1) * n * m + (j - 1) * n + i
      }
    }

    # check for the special case of the game being constructed over real values only
    if (all(unlist(lapply(losses, function(x) { is.numeric(x) && length(x)==1L })))) {
      losses <- unlist(losses)
      minPayoff <- min(losses)
      # normalize all payoffs to be >= 1
      losses <- losses - minPayoff + 1
      losses <- losses / (10^(1+max(ceiling(log10(losses)))))
      # emulate a scalar payoff by a Bernoulli distribution
      losses <- lapply(losses, function(x) {
        lossDistribution(c(1-x, x),
                         discrete=TRUE,
                         dataType="pdf",
                         smoothing="none",
                         bw = 1,
                         supp=c(1,2))
      })
    }

    mosg$losses <- losses
    mosg$dim <- goals
    if (is.null(defensesDescr)) {
      mosg$defensesDescriptions <- rep(1:n)

    }
    else {
      if (length(defensesDescr) != n) {
        stop("incompatible defense description vector")
      }
      mosg$defensesDescriptions <- defensesDescr

    }
    if (is.null(attacksDescr)) {
      mosg$attacksDescriptions <- rep(1:m)

    }
    else {
      if (length(attacksDescr) != m) {
        stop("incompatible attack description vector")
      }
      mosg$attacksDescriptions <- attacksDescr

    }
    if (is.null(goalDescriptions)) {
      mosg$goalDescriptions <- rep(1:goals)

    }
    else {
      if (length(goalDescriptions) != goals) {
        stop("incompatible goal description vector")
      }
      mosg$goalDescriptions <- goalDescriptions

    }

    # check if the list of payoff distributions has elements of the proper type
    #for(ld in losses) {
    #  if (class(ld) != "mosg.lossdistribution") {
    #    stop("improper elements in list of losses found (must all be of class 'mosg.lossdistribution'")
    #  }
    #}
    if (!all(unlist(lapply(losses, "class"))=="mosg.lossdistribution")) {
      stop("improper elements in list of losses found (must all be of class 'mosg.lossdistribution')")
    }

    # check if all losses have the same range (applies only to discrete distributions)
    N <- length(losses)
    discr <- rep(FALSE, N)
    for (i in 1:N) {
      discr[i] <- losses[[i]]$is.discrete
    }
    if (any(discr == FALSE) && any(discr == TRUE)) {
      stop("a mix of discrete and continuous distributions is not allowed")
    }
    else {
      if (discr[1]) {
        # the support is only required for discrete distributions
        supp <- losses[[1]]$supp
        if (length(losses) > 1) {
          for (i in 2:length(losses)) {
            if (compare::compare(losses[[i]]$supp, supp)$result != TRUE) {
              stop("discrete loss distributions must all be supported on the same range")
            }
          }
        }
        # for categorical distributions, we need to update the maximal loss to be the right end of the
        # support
        mosg$maximumLoss <- max(supp)
      }
      else {
        # game is over continuous distributions
        # determine maximal loss as the default cut off point
        # (saves resources, since the FP procedure does not have to compute its default value upon each call)
        maxloss <- 0
        for (p in losses) {
          maxloss <- max(maxloss, max(p$observations))
        }
        mosg$maximumLoss <- maxloss
      }
    }

    class(mosg) <- "mosg"
    mosg
  }
