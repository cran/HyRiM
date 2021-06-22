lossDistribution <-
function(dat,
           discrete = FALSE,
           dataType = c("raw", "pdf", "cdf"),
           supp = NULL,
           smoothing = c("none", "ongaps", "always"),
           bw = NULL) {
## The all-in-one include header for the HyRiM R package
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
     
    # note that we need the bandwidth perhaps for smoothing, so we compute it here
    if (is.null(bw)) {
      bw <- bw.nrd0(dat)
    }
    
    dens <- list()
    dens[[1]] <- NULL # reserve an empty element to be filled later
    dens[[2]] <- bw
    # TODO: perhaps rename this data field to name other possible content more accurately
    dens$observations <- dat
    
    names(dens)[c(1, 2)] <- c("lossdistr", "bw")
    dens$is.discrete <- discrete
    if (!discrete && !is.null(supp)) {
      warning("argument 'supp' is ignored unless the distribution is specified as discrete")
    }
    
    if (missing(dataType)) {
      dataType <- "raw"
    }
    if (discrete) {
      if (!is.na(pmatch("p", dataType))) {
        # the input data is a probability mass function; so directly take its values to the output
        if (is.null(supp)) {
          # default if no support is supplied externally
          dens$supp <- c(1, length(dat))
        }
        else {
          dens$supp <- supp
        }
        if (any(dat < 0)) {
          stop("negative probability masses not allowed")
        }
        if (sum(dat) != 1) {
          warning("renormalizing probability mass function")
          dat <- dat / sum(dat)
        }
        dens$range <- c(1, max(dens$supp))
        dens$dpdf <- rep(x = 0, dens$range[2])
        dens$dpdf[1:length(dat)] <- dat
      }
      if (!is.na(pmatch("c", dataType))) {
        # the input data is a cumulative distribution function; so derive its mass function
        if (is.null(supp)) {
          # default if no support is supplied externally
          dens$supp <- c(1, length(dat))
        }
        else {
          dens$supp <- supp
        }
        if (any(diff(dat) < 0)) {
          stop("cumulative distribution function must not be decreasing")
        }
        if (any(dat < 0)) {
          stop("cumulative distribution function cannot take negative values")
        }
        if (max(dat) > 1) {
          warning("renormalizing cumulative distribution function")
          dat <- dat / max(dat)
        }
        dens$range <- c(1, max(dens$supp))
        dens$dpdf <- rep(x = 0, dens$range[2])
        dens$dpdf[1:length(dat)] <- c(dat[1], diff(dat))
      }
      if (!is.na(pmatch("r", dataType))) {
        # note that, however, extra treatment is required for plotting, since we want to present a
        # nicely looking bar plot to the user (rather than the "ugly" internal approximation)
        if (is.null(supp)) {
          # default if no support is supplied externally
          dens$supp <- range(dat)
        }
        else {
          dens$supp <- supp
        }
        if (min(dens$supp) < 1) {
          stop("loss distribution must be supported within the interval [1,infinity)")
        }
        dens$dpdf <-
          hist(dat,
               plot = FALSE,
               breaks = seq(
                 from = min(dens$supp) - 0.5,
                 to = max(dens$supp) + 0.5,
                 by = 1
               ))$density
        dens$range <- range(dat)
      }
      # if the distribution contains gaps, then see how the user wants them to be treated
      if (length(smoothing) == 1) {
        # if omitted, use the default, which is "none"
        if (!is.na(pmatch("a", smoothing)) ||
            (!is.na(pmatch("o", smoothing)) && any(dens$dpdf == 0))) {
          
          # local functions for smoothing using discretized Gaussian kernels
          dkde <- function(n, h) {
            return(pnorm(
              q = n + 0.5,
              mean = 0,
              sd = h
            ) - pnorm(
              q = n - 0.5,
              mean = 0,
              sd = h
            ))
          }
          dkdsmooth <- function(x, h) {
            n <- length(x)
            y <- dkde((-5 * n):(5 * n), h)
            L <- 5 * n + 1 # = (length(y) - 1) / 2 + 1
            s <-
              convolve(x, y, type = "open")
            s <-
              s[L:(L + length(x) - 1)] # cut out the original region
            return(s / sum(s)) # re-normalization
          }
          
          dens$dpdf <- dkdsmooth(dens$dpdf, bw)
          dens$range <- c(1,max(which(dens$dpdf > 0)))
        }
      }
      # the loss distribution is merely the constructed probability mass function
      dens$lossdistr <- function(x) {
        return(dens$dpdf[x])
      }
      
    }
    else {
      # the loss distribution is continous
      if (min(dat) < 1) {
        stop("all observations must be >= 1")
      }
      # construct kernel density estimate explicitly
      loss <- function(z, dat) {
        n <- length(dat)
        lossfactor <- 1 / (n * bw * sqrt(2 * pi))
        loss <- 0
        N <- length(dat)
        bwsqh <- -0.5 / (bw ^ 2)
        for (i in 1:N) {
          loss <- loss + exp(bwsqh * (z - dat[i]) ^ 2)
        }
        loss <- lossfactor * loss
        loss
      }
      dens$lossdistr <- function(x) {
        return(loss(x, dat))
      }
      # normalization factor is required, since we truncate the density's lower tail at x=1
      leftLimit <- min(dat) - 10*bw
      dens$normalizationFactor <- integrate(dens$lossdistr,
                      lower = max(1,leftLimit), #lower = 1,
                      upper = max(dat) + 10 * bw)$value
      if (dens$normalizationFactor < .Machine$double.eps) {
        stop("numerical problem: loss distribution assigns approximately zero mass up to thespecified support (consider shrinking the interval on which the loss distributions are defined")
      }
      dens$normalizationFactor <- 1 / dens$normalizationFactor 
      dens$range <- range(dat)
    }
    
    dens$is.mixedDistribution <- FALSE
    class(dens) <- "mosg.lossdistribution"
    dens
  }
