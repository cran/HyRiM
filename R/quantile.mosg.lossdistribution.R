quantile.mosg.lossdistribution <-
function(x, p, eps = 0.001, ...) {
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
   if (x$is.discrete) {
    cs <- cumsum(x$dpdf)
    return(which(cs >= p)[1])
  }
  else {
    r <- x$range
    # we do a bisective search over the support of the distribution
    d <- r[2] - r[1]
    leftLimit <- max(1, r[1] - 5*x$bw)
    while (d > eps) {
      m <- (r[1] + r[2]) / 2
      A <-
        integrate(x$lossdistr, lower = leftLimit, upper = m)$value * x$normalizationFactor
      if (A > p) {
        r[2] <- m
      }
      else {
        r[1] <- m
      }
      d <- r[2] - r[1]
    }
    return((r[2] + r[1]) / 2)
  }
}
