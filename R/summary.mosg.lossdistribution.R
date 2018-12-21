summary.mosg.lossdistribution <-
function(object, ...) {
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
   S <- NULL
  S$range <- object$range
  S$mean <- mean(object)
  S$variance <- variance(object)
  q <-
    matrix(
      c(
        quantile(object, 0.1),
        quantile(object, 0.25),
        quantile(object, 0.5),
        quantile(object, 0.75),
        quantile(object, 0.9)
      ),
      nrow = 1,
      ncol = 5
    )
  colnames(q) <- c("10%", "25%", "50%", "75%", "90%")
  S$quantiles <- q
  S$is.discrete = object$is.discrete
  class(S) <- "summary.mosg.lossdistribution"
  S
}
