
mgss <-
  function(G,
           weights,
           cutOff,
           ord = 5,
           fbr = FALSE,
           points = 512,
           tol = 0.0) {
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
    d <- G$dim
    if (missing(cutOff)) {
      M <- G$maximumLoss
    }
    else {
      M <- cutOff
    }
    n <- G$nDefenses
    m <- G$nAttacks # loss matices have dimension nxm
    if (missing(weights)) {
      w <- rep(1 / G$dim, G$dim)
    }
    else {
      if (any(weights <= 0)) {
        stop("weights must all be > 0")
      }
      if (length(weights) < G$dim) {
        stop("less weights than goals specified")
      }
      if (sum(weights) == 0) {
        stop("weights cannot all be zero")
      }
      w <- weights / sum(weights)
    }

    #check for all distributions being not mixed
    if (any(unlist(lapply(G$losses, as.function(
      alist(x = , x$is.mixedDistribution)
    ))))) {
      stop("games constructed from mixed distributions are (currently) not supported")
    }

    # if the game is over categorical distributions, verify that no empty bins are in there
    # note that can intentionally avoid the more elegant "for(ld in losses)" notation here, and instead
    # resort to the slower version invoking the loc-function for a more informative warning to the user indicating
    # *where* in the game the faulty loss distribution was found
    for (ld in G$losses) {
      if (any(ld$dpdf == 0)) {
        stop(cat("categorical distributions with empty categories are not allowed. Consider reorganizing the game by smoothing the loss distributions."))
      }
    }

    if (G$losses[[1]]$is.discrete) {
      N <- min(M, G$losses[[1]]$range[2])  # all distributions are assured to have the same support
      ord <- N - 1
    }
    Ai <- .vectorize_lexcomparable(G,M,ord)

    #### construct auxiliary game
    # new loss matrix A for defender (averaged losses, single goal = weighted sum of all original goals)
    A <- list()
    for (i in 1:n) {
      A[[i]] <- list()
      for (j in 1:m) {
        A[[i]][[j]] <- w[1] * Ai[[1]][[i]][[j]]
        if (d > 1) {
          for (p in 2:d) {
            A[[i]][[j]] <- A[[i]][[j]] + w[p] * Ai[[p]][[i]][[j]]
          }
        }
      }
    }

    # below, we are going to loop from 2:n, and 2:m respectively, but we need to prevent
    # R from running a downward loop for degenerate games having either m = 1 or n = 1
    if (m > 1) { mRange <- 2:m } else { mRange <- NULL }
    if (n > 1) { nRange <- 2:n } else { nRange <- NULL }


    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # compute security strategy by linear programming
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    o <- ord+1
    # re-arrange the flattened list of payoffs into a list of matrices per coordinate for a lexicographic optimization
    Ui <- purrr::map(1:o,   # loop over all coordinates (in variable k)
                     function(k) {   # take out only the k-th coordinate to form a distinct game matrix
                       matrix(sapply(
                         unlist(A, recursive = FALSE), function(x) { x[k] }),   # extract the k-th coordinate
                         nrow = n, ncol = m, byrow = TRUE)   # construct the payoff matrix; note that "A" is *always* filled by row in the code above
                     })

    # shift all values into the positive range to avoid errors with the linear optimization
    delta <- abs(min(sapply(Ui, min))) + 1  # get the smallest value and
    Ui <- lapply(Ui, "+", delta)

    vopt <- rep(0, times=o)
    n <- nrow(Ui[[1]])
    m <- ncol(Ui[[1]])

    f.con <- matrix(rep(1, times=(n+1)*(m+1)), nrow=m+1, ncol=n+1)
    f.con[m+1, 1] <- 0
    f.obj <- c(1, rep(0, times=n))
    f.dir <- c(rep(">=", times = m), "==")
    f.rhs <- c(rep(0, times=m), 1)

    for(p in 1:o) {
      f.con[1:m, 2:(n+1)] <- -t(Ui[[p]])
      if (p > 1) {
        # update the LP for the next game
        A2 <- cbind(rep(0, times=m), t(Ui[[p-1]]))
        f.con <- rbind(f.con, A2)
        # objective remains the same ... no change to f.obj
        f.dir <- c(f.dir, rep("<=", times = m))
        #f.rhs <- c(f.rhs, rep(vopt[p-1], times = m))
        f.rhs <- c(f.rhs, rep(vopt[p-1] + tol, times = m))
      }

      result <- Rglpk::Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, max = FALSE,
                                      control = list("canonicalize_status" = FALSE))
      if (result$status != 5) {
        stop(paste("internal error: solver for linear program failed with GLPK status",
                   result$status,
                   ". Perhaps re-try with a tolerance (e.g., tol=1.0e-5)!?"))
      }
      vopt[p] <- result$optimum
    }
    optimalDefense <- as.matrix(result$solution[-1])   # this was the computed security strategy

    # having found the security strategy for the defender,
    # we compute the worst-case attacks as individually optimal relies (using LP again)
    # NOTE that this is *not* a best reply to the actual defense of player 1, but rather the
    # best that the attacker could achieve in engaging in its own zero-sum game with the defender
    # as if it were the defender's *only* goal to protect itself against this adversary. Since
    # in reality, the defender is caring for several goals at the same time, this worst-case attack is
    # actually understood as the best that the attacker could achieve, given that the defender would
    # focus all its efforts on protecting this particular goal. This is essentially different to a
    # best reply to "optimalDefense" as computed above

    optimalAttacks <- matrix(rep(0, times = m*d), nrow=m, ncol=d)
    for (p in 1:d) {  # run over all goals (= opponents)
      # like for the defender, re-arrange the flattened list of payoffs into a list of matrices per coordinate
      Ui <- purrr::map(1:o,   # loop over all coordinates (in variable k)
                       function(k) {   # take out only the k-th coordinate to form a distinct game matrix
                         t(  # transpose to switch roles of the players (for the LP only)
                           matrix(sapply(
                             unlist(Ai[[p]], recursive = FALSE), function(x) { x[k] }),   # extract the k-th coordinate
                             nrow = n, ncol = m, byrow = TRUE)   # construct the payoff matrix; note that "A" is *always* filled by row in the code above
                         )
                       })
      # shift all values into the positive range, since the LP solver implicitly
      # assumes all variables to be >= 0 (and would otherwise report a failure on
      # finding a feasible solution to start with)
      delta <- abs(min(sapply(Ui, min))) + 1  # get the smallest value and shift accordingly
      Ui <- lapply(Ui, "+", delta)

      # reconstruct the LP like above, only swapping the variables n and m, and letting
      # the adversary be a maximizer (thus reversing the inequalities in the constraints
      # and the optimization goal)
      f.con <- matrix(rep(1, times=(m+1)*(n+1)), nrow=n+1, ncol=m+1)
      f.con[n+1, 1] <- 0
      f.obj <- c(1, rep(0, times=m))
      f.dir <- c(rep("<=", times = n), "==")
      f.rhs <- c(rep(0, times=n), 1)

      for(j in 1:o) {   # run over all coordinates for that opponent's game
        f.con[1:n, 2:(m+1)] <- -t(Ui[[j]])
        if (j > 1) {
          # update the LP for the next game
          A2 <- cbind(rep(0, times=n), t(Ui[[j-1]]))
          f.con <- rbind(f.con, A2)
          # objective remains the same ... no change to f.obj
          f.dir <- c(f.dir, rep(">=", times = n))
          f.rhs <- c(f.rhs, rep(vopt[j-1] - tol, times = n))
        }

        result <- Rglpk::Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, max = TRUE,
                                        control = list("canonicalize_status" = FALSE))  # attacker maximizes
        if (result$status != 5) {  # status 5 = optimal solution
          # occasionally, the failure was due to numeric roundoff errors below the 6th digit after
          # the comma, so we re-try with a bit of a tolerance subtracted (=> worked in all test cases)
          stop(paste("internal error: LP failed (for attacker); GLPK status was",
                     result$status,
                     ". Perhaps re-try with a tolerance (e.g., tol=1.0e-5)!?"))
        }
        vopt[j] <- result$optimum
      }
      optimalAttacks[,p] <- result$solution[-1]
    }

    ## Compile result object *********************

    ## estimated probabilities for defender (mixed strategies)
    ## compile equilibrium object for further use (with summary, print, etc.)
    equilibrium <- NULL
    colnames(optimalDefense) <- "prob."
    rownames(optimalDefense) <- G$defensesDescriptions
    rownames(optimalAttacks) <- G$attacksDescriptions
    colnames(optimalAttacks) <- G$goalDescriptions
    equilibrium$optimalDefense <- optimalDefense
    equilibrium$optimalAttacks <- optimalAttacks

    hybridRiskMetric <- list()
    assurance <- list()
    for (p in 1:d) {
      # loop over all goals, and construct the respective loss distribution per goal
      hybridRiskMetric[[p]] <-
        lossDistribution.mosg(G,
                              optimalDefense,
                              optimalAttacks[, p],
                              goal = p,
                              points = points)
    }
    names(hybridRiskMetric) <- G$goalDescriptions
    equilibrium$assurances <- hybridRiskMetric

    #######################################
    # skip this computation upon a flag with default value
    # additionally, compute the individually best replies to the fixed "optimalDefense"
    # this item is made available, but not printed explicitly
    br <- NA
    if (fbr) {
      br <- rep(0, times=d)
      zeroMassVector <- rep(0, length(Ai[[1]][[1]][[1]]))  # needed below (in both methods)

      for (p in 1:d) {

        # run over all rows
        Vmax <- zeroMassVector
        for (i in 1:n) {
          # note that there is no weight w[p] needed in the sum here
          # (as this player is interacting only with the defender)
          Vmax <- Vmax + optimalDefense[i] * Ai[[p]][[i]][[1]]  # take column 1 as first maximum-candidate
        }
        col <- 1
        for(j in mRange) { # run over all columns j=2,3,... for the p-th opponent (column j=1 was done in the loop above)
          # given the so-far recorded behavior of the defender,
          # determine the best attack (maximizing the loss)
          colPayoff <- zeroMassVector
          for (i in 1:n) { # run over all rows i
            colPayoff <- colPayoff + optimalDefense[i] * Ai[[p]][[i]][[j]]
          }
          if (.lexgt(colPayoff, Vmax)) {
            Vmax <- colPayoff # update the optimum
            # store the location of the current optimum
            # only if it *did change*
            col <- j
          }
        }

        br[p] <- col # store best reply for p-th opponent as index
      }
    }
    equilibrium$br_to_optimalDefense <- br
    #######################################

    class(equilibrium) <- "mosg.equilibrium"
    equilibrium
  }
