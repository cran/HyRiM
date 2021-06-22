plot.mosg.lossdistribution <-
function(x,
           points = 100,
           xlab = "",
           ylab = "",
           main = "",
           p = 0.999,
           newPlot = TRUE,
           cutoff = NULL,
           ...) {
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
     r <- x$range
    if (!newPlot) {
      par(new = TRUE)
    }
    if (x$is.discrete) {
      # auxiliary normalization factor for (second) truncation at cutoff point
      y <- x$dpdf # no truncation so far
      if (!is.null(cutoff)) {
        cutoff <- min(cutoff, length(x$dpdf))
        aux <- sum(x$dpdf[1:cutoff])
        y <- x$dpdf / aux
        if (cutoff < length(x$dpdf)) {
          y[(cutoff + 1):length(x$dpdf)] <- 0
        }
      }

      # if xlim and ylim are passed from outside, we shall use these
      # values to replace the default
      args <- list(...)
      if ("xlim" %in% attributes(args)$names) {
        xlimit <- args$xlim
      }
      else {
        xlimit <- c(0, length(x$dpdf) + 1)
      }
      if ("ylim" %in% attributes(args)$names) {
        ylimit <- args$ylim
      }
      else {
        ylimit <- c(0, max(y))
      }
      args[['height']] <- y;
      args[['names.arg']] <- 1:length(x$dpdf)
      args[['xlab']] <- xlab
      args[['ylab']] <- ylab
      args[['main']] <- main
      args[['xlim']] <- xlimit
      args[['ylim']] <- ylimit
      args[['font.main']] <- 1
      
      do.call(barplot, args)
    }
    else {
      r[2] <- quantile.mosg.lossdistribution(x, p)
      # auxiliary normalization factor for (second) truncation at cutoff point
      if (!is.null(cutoff)) {
        leftLimit <- x$range[1] - 5*x$bw
        rightLimit <- x$range[2] + 5*x$bw
        aux <- integrate(f = x[[1]],
                         lower = max(1,leftLimit), #lower = 1,
                         upper = min(cutoff,rightLimit))$value #upper = cutoff)$value
        if (aux < .Machine$double.eps) {
          warning(
            paste("numerical problem: loss distribution assigns approximately zero mass up to the cutoff point ", 
                  cutoff,
                  "\ndensity function will not be normalized for the plot")
          )
          aux <- 1
        }
        M <- cutoff
      }
      else {
        aux <- 1
        M <- r[2] + 5 * x$bw
      }
      h <- seq(r[1], M, length = points)
      y <- x[[1]](h) * x$normalizationFactor / aux
      
      # if xlim and ylim are passed from outside, we shall use these
      # values to replace the default
      args <- list(...)
      if ("xlim" %in% attributes(args)$names) {
        xlimit <- args$xlim
      }
      else {
        xlimit <- c(r[1], r[2] + 5 * x$bw)
      }
      if ("ylim" %in% attributes(args)$names) {
        ylimit <- args$ylim
      }
      else {
        ylimit <- range(y)
      }
      args[['x']] <- h
      args[['y']] <- y
      args[['type']] <- "l"
      args[['xlab']] <- xlab
      args[['ylab']] <- ylab
      args[['main']] <- main
      args[['xlim']] <- xlimit
      args[['ylim']] <- ylimit
      args[['font.main']] <- 1
      do.call(plot,args)

      if (!is.null(cutoff)) {
        abline(v = cutoff, lty = 2)
      }
    }
  }
