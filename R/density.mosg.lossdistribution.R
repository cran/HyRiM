density.mosg.lossdistribution <-
function(x, t, ...) {
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
    # set everything larger the support to zero and omit (crop) categories with indices <= 0 
    if (min(t) <= 0) {
      warning("categories with numbers < 1 have zero probability mass by default and will be cropped off from the output")
      # cut off any indices <= 0
      t <- t[t>0]
    }
    if (length(t) == 0) { return(NA) }
    if (max(t) > length(x$dpdf)) {
      warning("the given range exceeds the suppport (zero mass values are returend for these categories)")
    }
    if (any(floor(t) < ceiling(t))) {
      stop("density is only defined for integer arguments (categories)")
    }
    else
      result <- x$dpdf[as.integer(t)]
      result[is.na(result)] <- 0   # replace NAs by zero when they exist
      return(result)
  }
  else {
      result <- x$lossdistr(t) * x$normalizationFactor
      # set everything outside the support to zero
      result[which(!((x$range[1] - 5 * x$bw <= t) & (t <= x$range[2] + 5 * x$bw)))] <- 0
      return(result)
  }
}
