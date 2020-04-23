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
    hpoly <-
      orthopolynom::hermite.h.polynomials(ord, normalized = FALSE)# Hermite polynomials needed to compute derivatives of normal kernel
    # o-th derivative
    derivdens <- function(z, obs, o, bw) {
      N <- length(obs)
      derivdensfactor <-
        1 / (N * bw) / sqrt(2 * pi) * (-1 / (sqrt(2) * bw)) ^
        o
      derivdens <- 0
      hpoly <- orthopolynom::hermite.h.polynomials(ord, normalized = FALSE)
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
