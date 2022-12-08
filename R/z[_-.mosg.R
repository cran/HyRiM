`[<-.mosg` <-
function(x,i,j,k=NULL,value) {
  rows <- 1:x$nDefenses
  cols <- 1:x$nAttacks
  goals <- 1:x$dim
  if (!missing(i)) {
    # take care of subsetting "by name"
    #if (class(i) == "character") {
    if (is(i, "character")) {
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
    #if (class(j) == "character") {
    if (is(j, "character")) {
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
    #if (class(k) == "character") {
    if (is(k, "character")) {
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

  # determine if the game x has its payoff list enumerated by row or by
  # column; this piece of information is not carried over with x, but is easy
  # to find out from a test call to the "loc" function
  byrow = (x$loc(1, 1, 1) + 1 == x$loc(1,1,2))
  #if (class(value) == "mosg") {
  if (is(value, "mosg")) {
    # if we substitute using another game object,
    # check if its payoff list is organized as that of x
    if ((value$loc(1,1,1) + 1 == value$loc(1,1,2)) != byrow) {
      warning("payoff lists are differently structured (one is 'by rows', the other 'by columns'). Values will be inserted, but check the game structure afterwards!")
    }
    value <- value$losses
  }
  #if (class(value) != "list") {
  if (!is(value, "list")) {
    stop("need a 'list' object for insertion")
  }

  if (length(value) != length(rows) * length(cols) * length(goals)) {
    stop(paste("list of replacement values has not the correct size; should be", rows*cols*goals))
  }

  idx <- 0
  if (byrow) { # enumeration by row
    for(k in goals) {
      for(i in rows) {
        for(j in cols) {  # set up the new list row by row from the given game matrix
          idx <- idx + 1
          x$losses[[x$loc(k, i, j)]] <- value[[idx]]
        }
      }
    }
  } else {  # list by columns
    for(k in goals) {
      for(j in cols) {  # set up the new list col by col from the given game matrix
        for(i in rows) {
          idx <- idx + 1
          x$losses[[x$loc(k, i, j)]] <- value[[idx]]
        }
      }
    }
  }
  return(x)
}
