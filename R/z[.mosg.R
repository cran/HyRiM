`[.mosg` <-
function(x,i,j,k=NULL) {
  rows <- 1:x$nDefenses
  cols <- 1:x$nAttacks
  goals <- 1:x$dim
  if (!missing(i)) {
    # take care of subsetting "by name"
    if (class(i) == "character") {
      i <- match(i, x$defensesDescr)
      if (anyNA(i)) {
        stop("unknown defense strategy")
      }
    }
    rows <- rows[i]
    if (anyNA(rows)) {
      stop("unknown defense strategy")
    }
  }
  if (!missing(j)) {
    if (class(j) == "character") {
      j <- match(j, x$attacksDescr)
      if (anyNA(j)) {
        stop("unknown attack strategy")
      }
    }
    cols <- cols[j]
    if (anyNA(cols)) {
      stop("unknown attack strategy")
    }
  }
  if (!is.null(k)) {
    if (class(k) == "character") {
      k <- match(k, x$goalDescriptions)
      if (anyNA(i)) {
        stop("unknown goal(s)")
      }
    }
    goals <- goals[k]
    if (anyNA(goals)) {
      stop("unknown goal(s)")
    }
  }

  sublist <- list()
  for(k in goals) {
    for(i in rows) {
      for(j in cols) {  # set up the new list row by row from the given game matrix
        sublist <- append(sublist, list(x$losses[[x$loc(k, i, j)]]))
      }
    }
  }
  # determine if the game x has its payoff list enumerated by row or by
  # column; this piece of information is not carried over with x, but is easy
  # to find out from a test call to the "loc" function
  # the returned game will be organized in the same fashion
  byrow = (x$loc(1, 1, 1) + 1 == x$loc(1,1,2))
  Gsub <- mosg(n = length(rows), m = length(cols), goals = length(goals),
               losses = sublist,
               byrow = byrow,
               goalDescriptions = x$goalDescriptions[goals],
               defensesDescr = x$defensesDescriptions[rows],
               attacksDescr = x$attacksDescriptions[cols])
  return(Gsub)
}
