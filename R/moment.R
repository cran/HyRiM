moment <-
function(ld, k) {
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
   if (as.integer(k) < 1) {
    stop("moment order must be an integer >= 1")
  }
  if (ld$is.discrete) {
    return(sum(seq(
      from = ld$supp[1], to = ld$supp[2]
    ) ^ k * ld$dpdf))
  }
  else {
    f <- function(x) {
      x ^ k * ld$lossdistr(x) * ld$normalizationFactor
    }
    leftLimit <- ld$range[1] - 5*ld$bw
    # "integrate" occasionally has issues when the limit is taken as infinite,
    # so we avoid this by setting the integration limit up to the maximal observation + 10 times the bandwidth
    # (bearning in mind that the underlying kernel densities are Gaussian)
    return(integrate(f, lower = max(1,leftLimit), upper = ld$range[2] + 10*ld$bw)$value)
    #return(integrate(f, lower = max(1,leftLimit), upper = Inf)$value)
  }
}
