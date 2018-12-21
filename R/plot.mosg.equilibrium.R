plot.mosg.equilibrium <-
function(x, points = 100, ...) {
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
   d <- length(x$assurances)
  goalNames <- names(x$assurances)
  plot.new()
  layout(matrix(c(1, 1, 2:(2 * d + 1)), d + 1, 2, byrow = TRUE))
  # plot optimal defense strategy
  barplot(x$optimalDefense[, 1], ylab = "prob.", sub = "optimal defense")
  for (p in 1:d) {
    barplot(
      x$optimalAttacks[, p],
      ylab = "prob.",
      sub = paste("worst case for goal '", goalNames[p], "'", sep = "")
    )
    plot(x$assurances[[p]], ...)#, sub=paste("loss distribution; goal '", goalNames[p], "'", sep=""))
  }
}
