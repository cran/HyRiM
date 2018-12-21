lossDistribution.mosg <-
function(G,
           player1Strat,
           player2Strat,
           points = 512,
           goal = 1) {
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
     n <- length(G$losses)
    low <- Inf
    up <- -Inf
    maxbw <- -Inf
    nm <- G$nAttacks * G$nDefenses
    for (i in ((goal - 1) * nm + 1):(goal * nm)) {
      ld <- G$losses[[i]]
      low <- min(low, ld$range[1])
      up <- max(up, ld$range[2])
      maxbw <- max(maxbw, ld$bw)
    }
    isDiscrete <-
      G$losses[[1]]$is.discrete  # if the first distribution is discrete, then all others are too (has been checked by the 'mosg' function)
    
    result <- list()
    if (isDiscrete) {
      # note that all discrete distributions must have the same support
      ran <- G$losses[[1]]$range
      v <- rep(0, length = ran[2] - ran[1] + 1)
      for (i in which(player1Strat > 0)) {
        for (j in which(player2Strat > 0)) {
          ld <- G$losses[[G$loc(goal, i, j)]]
          v <-
            v + player1Strat[i] * player2Strat[j] * ld$dpdf
        }
      }
      result[[1]] <- stepfun(1:length(v), c(v, 0))
    }
    else {
      h <- seq(low, up + 5 * maxbw, length = points)
      v <- rep(0, length = points)
      
      for (i in which(player1Strat > 0)) {
        for (j in which(player2Strat > 0)) {
          ld <- G$losses[[G$loc(goal, i, j)]]
          v <-
            v + player1Strat[i] * player2Strat[j] * ld$lossdistr(h)
        }
      }
      result[[1]] <- function(z) {
        w <- splinefun(h, v, method = "natural")(z)
        w[w < 0] <- 0
        w
      }
    }
    
    result[[2]] <- maxbw
    result[[3]] <- c(low, up)
    names(result) <- c("lossdistr", "bw", "range")
    result$is.mixedDistribution <- TRUE
    class(result) <- "mosg.lossdistribution"
    result$is.discrete = isDiscrete
    if (isDiscrete) {
      result$dpdf <- v
    }
    else {
      leftLimit <- result$range[1] - 5*maxbw
      # for some reason, integration up to infinity runs into numeric roundoff troubles, so we take the integration limit up to 10 times the bandwidth, from the last (maximal) observation
      result$normalizationFactor <- integrate(result$lossdistr, lower = max(1,leftLimit), upper = result$range[2] + 10*maxbw)$value
      # ideally, here we should never end up with a zero mass returned by "integrate",
      # because this would elsewhere have been detected (hopefully). Nonetheless,
      # we issue a warning here (just in case)
      if (result$normalizationFactor < .Machine$double.eps) {
        warning("mixed loss distribution cannot be normalized to unit mass (numeric integration error); density will be returned non-normalized")
        result$normalizationFactor <- 1
      }
      result$normalizationFactor <- 1 / result$normalizationFactor
    }
    result$supp <- seq(from=result$range[1], to=result$range[2])
    return(result)
  }
