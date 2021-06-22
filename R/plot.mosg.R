plot.mosg <-
function(x,
           goal = 1,
           points = 100,
           cutoff = NULL,
           largeGame = FALSE,
           subPlotWidth=2,
           subPlotHeight=2,
           cleanUp = TRUE,
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
    n <- x$nDefenses
    m <- x$nAttacks
    if (largeGame == TRUE) {
      tmpSVG <- tempfile(fileext = ".svg")
      svg(filename = tmpSVG, width = n * subPlotHeight, height = m * subPlotWidth)
    }
    par(mfrow = c(n, m))
    for (i in 1:n) {
      for (j in 1:m) {
        k <-
          x$loc(goal, i, j)  # account for specification of game matrix "by row" or "by column"

        if (j == 1)
          rowLabel <- x$defensesDescriptions[i]
        else
          rowLabel <- ""
        if (i == 1)
          colLabel <- x$attacksDescriptions[j]
        else
          colLabel <- ""

        plot(
          x$losses[[k]],
          xlab = paste("loss(", i, ",", j, ")", sep = ""),
          ylab = rowLabel,
          main = colLabel,
          points,
          cutoff = cutoff,
          ...
        )
        title(main = x$goalDescriptions[[goal]], outer = TRUE)
      }
    }
    if (largeGame == TRUE) {
      dev.off()
      pic <- readPicture(tmpSVG)
      grid.picture(pic)
      # cleanup temp file (no longer needed)
      if (cleanUp && file.remove(tmpSVG) == FALSE) {
        warning(paste("failed to cleanup temporary file ", tmpSVG))
      }
      if (!cleanUp) {
        message("temporary file left for further use: ", tmpSVG)
      }
    }
  }
