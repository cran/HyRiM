preference <-
function(x,
                       y,
                       verbose = FALSE,
                       weights,
                       points = 512) {
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
   fcall <- match.call(expand.dots = FALSE)
  # trivial case: compare two degenerate distributions
  if (class(x) == "numeric" && class(y) == "numeric") {
    if (verbose) { return(min(x,y)) }
    else {
      if (x < y) { return(1) }
      else if (x > y) { return(2) }
      return(0) # upon identity
    }
  }
  
  # simple case: compare a distribution to a number
  if ((class(x) == "mosg.lossdistribution" && class(y) == "numeric") ||
      (class(x) == "numeric" && class(y) == "mosg.lossdistribution")) {
    
    if (class(x) == "numeric" && class(y) == "mosg.lossdistribution") {
      # switch x and y accordingly so that x is the distribution and y is the number
      z <- x
      x <- y
      y <- z
      bSwitch <- TRUE
      xName <- fcall[3]
      yName <- fcall[2]
    }
    else {
      bSwitch <- FALSE
      xName <- fcall[2]
      yName <- fcall[3]
    }
    # we add the bandwidth only for continuous distributions
    if (x$is.discrete) {
      xRange <- x$supp[2]
    }
    else {
      xRange <- x$range + 5*x$bw
    }
    if (xRange < y) {
      if (verbose) {
        #cat(paste("\npreferred distribution:", xName, "\n\n"))
        return(x) # return full distribution object(s)
      }
      return(ifelse(bSwitch, 2, 1))  # return the argument index (reversing the argument switch if it happened)
    }
    else {
      if (verbose) {
        #cat(paste("\npreferred distribution:", yName, "\n\n"))
        return(y) # return full distribution object(s)
      }
      return(ifelse(bSwitch, 1, 2)) # return argument index only (reversing the argument switch if it happened)
    }
  }
  
  # to handle multiple goals, including the case of two distributions only, 
  # we create a list and flatten it. In case that x and y are loss distributions,
  # the resulting list will be singleton. In case of x and y being lists, the resulting
  # list will be a plain list again. So, in any case, we can iterate over the elements easily
  if (class(x) == "mosg.lossdistribution" && class(y) == "mosg.lossdistribution") {
    xList <- list(x)
    yList <- list(y)
  }
  else {
    xList <- x
    yList <- y
  }
  
  n <- length(xList)
  if (n != length(yList)) {
    stop("number of criteria differs between 'x' and 'y'")
  }
  if (missing(weights)) {
    weights <- rep(1/n, n)
  }
  else {
    if (length(weights) != n) {
      stop("number of weights must equal number of criteria")
    }
    if (any(weights <= 0)) {
      stop("weights must all be > 0")
    }
  }
  for(i in 1:length(xList)) {
    xi <- xList[[i]]
    yi <- yList[[i]]
    # nontrivial case: compare distributions to one another
    if (class(xi) != "mosg.lossdistribution" ||
        class(yi) != "mosg.lossdistribution") {
      stop("preferences can only be computed between loss distributions")
    }
    if (xi$is.discrete != yi$is.discrete) {
      stop("comparison between categorical and continuous distributions is not supported")
    }
    if (xi$is.discrete) {
      if (any(xi$supp - yi$supp != 0)) {
        stop("distributions xi and yi must be supported on the same set")
      }
      h <- xi$supp[1]:xi$supp[2]
    }
    else {
      rx <- xi$range + 5 * c(-xi$bw, xi$bw)
      ry <- yi$range + 5 * c(-yi$bw, yi$bw)
      h <-
        seq(
          from = min(rx[1], ry[1]),
          to = max(rx[2], ry[2]),
          length.out = points
        )
    }
    if (i == 1) { delta <- weights[[i]] * (density(xi, h) - density(yi, h)) }
    else {
      delta <- delta + weights[[i]] * (density(xi, h) - density(yi, h))
    }
  }
  
  if (all(delta == 0)) {
    if (verbose) {
      #cat("\ndistributions are identical\n")
      return(x)
    }
    return(0)  # distributions are identical
  }
  
  lastNonzero <- max(which(delta != 0))
  if (delta[lastNonzero] < 0) {
    if (verbose) {
      #cat(paste("\npreferred distribution:", fcall[2], "\n\n"))
      return(x) # return the full object(s)
    }
    return(1) # return the argument index only
  }
  if (verbose) {
    #cat(paste("\npreferred distribution:", fcall[3], "\n"))
    return(y) # return the full object(s)
  }
  return(2)
}
