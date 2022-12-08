# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# internal helper functions, not exported as part of the API
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to avoid the following issue when checking "mgss":
#mgss: no visible binding for global variable 'x'
#Undefined global functions or variables:
#  x
utils::globalVariables(c("x"))  # to avoid issues with

# predicate x "lexicographically greater than" y
.lexgt <- function(x, y) {
  tolerance <- 1e-10  # just to account for numeric roundoff errors
  if (length(x) == 1) {
    # the tolerance is only to account for anyway negligible roundoff errors...
    return(x > y + tolerance)
  }
  return((x[1] > y[1]) || (x[1] == y[1] && .lexgt(x[-1], y[-1])))
}

# predicate x "lexicographically less than" y
.lexlt <- function(x, y) {
  tolerance <- 1e-10  # just to account for numeric roundoff errors
  if (length(x) == 1) {
    # the tolerance is only to account for anyway negligible roundoff errors...
    return(x < y - tolerance)
  }
  return((x[1] < y[1] - tolerance) || (x[1] == y[1] && .lexlt(x[-1], y[-1])))
}

# convert a game's payoff structure into a set of vectors so that the distributions become
# lex-order comparable using this data. This means that:
# 1) discrete distributions go into the vector directly with their probability masses, only in reversed
#    order
# 2) continuous distributions are replaced by the derivatives of their density function, up to order
#    "ord", supplied as a parameter
# Additionally, the function truncates the distributions at x=M as a preprocessing step
.vectorize_lexcomparable <- function(G,M,ord) {
  d <- G$dim
  n <- G$nDefenses
  m <- G$nAttacks

  if (!G$losses[[1]]$is.discrete) {
    ###### compute Taylor approximations for loss distributions
    #hpoly <-
    #  orthopolynom::hermite.h.polynomials(ord, normalized = FALSE)# Hermite polynomials needed to compute derivatives of normal kernel
    # o-th derivative
    # NOTE: due to package checking issues with "orthopolynom", the relevant code from this package has been transferred hereto (with fixes)
    derivdens <- function(z, obs, o, bw) {
      N <- length(obs)
      derivdensfactor <-
        1 / (N * bw) / sqrt(2 * pi) * (-1 / (sqrt(2) * bw)) ^
        o
      derivdens <- 0
      hpoly <- .hermite_h_polynomials(ord, normalized = FALSE)
      # for efficiency, pull out factors that do not depend on i in the following summation
      c1 <- -1 / (2 * bw ^ 2)
      c2 <- 1 / (sqrt(2) * bw)
      H <- as.function(hpoly[[o + 1]])
      for (i in 1:N) {
        derivdens <-   derivdens + exp(c1 * (z - obs[i]) ^ 2) * H((z - obs[i]) *
                                                                    c2)
      }
      derivdens <- derivdensfactor * derivdens
      derivdens
    }
    #### transform densities to Taylor approximations (player p)
    Ai <- list()
    for (p in 1:d) {
      Ai[[p]] <- list()

      for (i in 1:n) {
        Ai[[p]][[i]] <- list()
        for (j in 1:m) {
          # one scenario (one specific choice of PS1 & PS2)
          # compute index of function in the loss vector, i.e. losses[[k]][[1]] is density for scenario (i,j), losses[[k]][[2]] is its estimated bandwidth
          k <-
            G$loc(p, i, j)  # account for specification of loss matrix "by row" or "by column"
          Ai[[p]][[i]][[j]] <-
            c(rep(0, ord + 1)) # will contain derivatives up to order o
          bw <- G$losses[[k]][[2]]
          # adapt the normmalization factor to account for the truncation of the mass function at x=M

          # note that we need to be careful with the numeric integration here to let it cover only
          # the support of the distribution, but not an unnecessarily wide interval overlapping it
          # (otherwise, the integral may incorrectly come back as 0, although the full mass would be
          # located inside the specified region). Hence, we adapt the lower and upper limits with the
          # support of the distribution, respectively.
          supp <- G$losses[[k]]$range
          leftLimit <- supp[1] - 5*G$losses[[k]]$bw
          rightLimit <- supp[2] + 5*G$losses[[k]]$bw
          truncationFactor <- integrate(G$losses[[k]]$lossdistr,
                                        lower = max(1,leftLimit),
                                        upper = min(M,rightLimit))$value
          # account for the integration possibly returning a zero (in which case we end up with NaNs later)
          if (truncationFactor < .Machine$double.eps) {
            warning(
              paste("numerical problem: loss distribution assigns approximately zero mass up to the cutoff point ", M, " (consider a different point to cut off or try eliminating dominated strategies)")
            )
          }
          for (o in 1:(ord + 1)) {
            if (truncationFactor < .Machine$double.eps) {
              Ai[[p]][[i]][[j]][o] <- Inf
            }
            else {

              Ai[[p]][[i]][[j]][o] <-
                derivdens(M, G$losses[[k]]$obs, o - 1, bw) / truncationFactor # (o-1)-th derivative at M
            }
          }
        }
      }
    }
  } # end of preprocessing for continuous distributions
  else {
    # for discrete distributions, the preprocessing is much simpler
    N <-
      min(M, G$losses[[1]]$range[2])  # all distributions are assured to have the same support
    # note that by this point, we have assured that all supports are equal
    # (so that the shape of Ai is consistent)
    Ai <- list()
    for (p in 1:d) {
      Ai[[p]] <- list()
      for (i in 1:n) {
        Ai[[p]][[i]] <- list()
        for (j in 1:m) {
          # one scenario (one specific choice of PS1 & PS2)
          # compute index of function in the loss vector, i.e. losses[[k]][[1]] is density for scenario (i,j), losses[[k]][[2]] is its estimated bandwidth
          k <-
            G$loc(p, i, j)  # account for specification of loss matrix "by row" or "by column"

          # adapt the normmalization factor to account for the truncation of the mass function at x=M
          truncationFactor <- 1 / sum(G$losses[[k]]$dpdf[1:N])

          # note that we reverse the mass vector here to check the
          # lex-ordering from the first coordinate (below)
          Ai[[p]][[i]][[j]] <-
            rev(G$losses[[k]]$dpdf[1:N]) * truncationFactor
        }
      }
    }
    ord <-
      N - 1  # take care for the code to assume the lists to have a length = ord + 1
  } # end of preprocessing for discrete distributions
  return(Ai)
}

################################################################################
# the following code was directly copied from "orthopolynom" package,
# version 1.0-5, licensed under GPL (>=2)
# this is for fixing issues with the orthopolynom package, as published on https://cran.r-project.org/web/checks/check_results_orthopolynom.html
.hermite_h_recurrences <- function (n, normalized = FALSE)
{
  if (n < 0)
    stop("negative highest polynomial order")
  if (n != round(n))
    stop("highest polynomial order is not an integer")
  np1 <- n + 1
  r <- data.frame(matrix(nrow = np1, ncol = 4))
  names(r) <- c("c", "d", "e", "f")
  j <- 0
  k <- 1
  if (normalized) {
    while (j <= n) {
      r[k, "c"] <- 1
      r[k, "d"] <- 0
      r[k, "e"] <- sqrt(2/(j + 1))
      if (j == 0) {
        r[k, "f"] <- 0
      }
      else {
        if (k == 1) {
          r[k, "f"] <- 0
        }
        else {
          r[k, "f"] <- sqrt(j/(j + 1))
        }
      }
      j <- j + 1
      k <- k + 1
    }
    return(r)
  }
  else {
    while (j <= n) {
      r[k, "c"] <- 1
      r[k, "d"] <- 0
      r[k, "e"] <- 2
      r[k, "f"] <- 2 * j
      j <- j + 1
      k <- k + 1
    }
    return(r)
  }
  return(NULL)
}

.orthogonal_polynomials <- function (recurrences)
{
  #require(polynom)   # this was required to be removed (see the above URL)
  np1 <- nrow(recurrences)
  n <- np1 - 1
  c <- recurrences$c
  d <- recurrences$d
  e <- recurrences$e
  f <- recurrences$f
  polynomials <- as.list(rep(NULL, np1))
  p.0 <- polynomial(c(1))
  polynomials[[1]] <- p.0
  j <- 0
  while (j < n) {
    cj <- c[j + 1]
    dj <- d[j + 1]
    ej <- e[j + 1]
    fj <- f[j + 1]
    monomial <- polynomial(c(dj, ej))
    if (j == 0) {
      p.jp1 <- (monomial * p.0)/cj
    }
    else {
      p.jm1 <- polynomials[[j]]
      p.j <- polynomials[[j + 1]]
      p.jp1 <- (monomial * p.j - fj * p.jm1)/cj
    }
    polynomials[[j + 2]] <- p.jp1
    j <- j + 1
  }
  return(polynomials)
}

.hermite_h_polynomials <- function (n, normalized = FALSE)
{
  recurrences <- .hermite_h_recurrences(n, normalized)
  # we do not need the polynomials to be normalized in this code
#  if (normalized) {
#    h.0 <- sqrt(pi)
#    p.0 <- polynomial(c(1/sqrt(h.0)))
#    polynomials <- .orthonormal_polynomials(recurrences, p.0)
#  }
#  else
  polynomials <- .orthogonal_polynomials(recurrences)
  return(polynomials)
}
