print.mosg.equilibrium <-
function(x, extended = FALSE, ...) {
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
  cat("\n\tequilibrium for multiobjective security game (MOSG)\n\n")
  cat("optimal defense strategy:\n")
  print(x$optimalDefense, ...)
  cat("\nworst case attack strategies per goal:\n")
  print(x$optimalAttacks, ...)
  if (extended) {
    cat("\nassurances:\n")
    print(x$assurances, ...)
  }
  if (!is.na(x$br_to_optimalDefense)) {
    tmp <- matrix(x$br_to_optimalDefense, nrow=1, ncol=length(x$br_to_optimalDefense))
    if (!is.null(x$goalDescriptions)) {
      colnames(tmp) <- x$goalDescriptions
    }
    else {
      colnames(tmp) <- 1:ncol(tmp)
    }
    cat(paste("\nbest replies per goal, following the (fixed) optimal defense strategy:\n"))
    print(tmp, ...)
  }
}
