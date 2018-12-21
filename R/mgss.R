mgss <-
function(G,
           weights,
           cutOff,
           ord = 5,
           eps,
           T = 1000,
           points = 512) {
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
     if (!missing(eps)) {
      if (eps < 0)
        stop("invalid precision (eps > 0 required)")
    }
    if (!missing(T)) {
      if (T < 0)
        stop("invalid precision (iteration count > 0 required)")
    }
    
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
    #for(k in G$dim) {
    #  for(i in G$nDefenses) {
    #    for(j in G$nAttacks) {
    #      ld <- G$losses[[G$loc(goal = k, i = i, j = j)]]
          if (any(ld$dpdf == 0)) {
            #stop(cat("categorical distributions with empty categories are not allowed (found at (", i, ",", j, ") for goal ", k, ". Consider reorganizing the game by smoothing the loss distributions."))
            stop(cat("categorical distributions with empty categories are not allowed. Consider reorganizing the game by smoothing the loss distributions."))
          }
    #    }
    #  }
    }
    
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
            
            
            Ai[[p]][[i]][[j]] <-
              rev(G$losses[[k]]$dpdf[1:N]) * truncationFactor
          }
        }
      }
      ord <-
        N - 1  # take care for the code to assume the lists to have a length = ord + 1
    } # end of preprocessing for discrete distributions
    
    
    ### compare two densities g,h represented by their Taylor expansion and find min / max
    mindens <- function(g, h) {
      if (compare::compare(g, h)$result == TRUE) {
        return(g) # if identical take any of them as min
      }
      else{
        ord = length(g)
        for (i in 1:(ord + 1)) {
          if (g[i] < h[i])
            return(g)
          if (g[i] > h[i])
            return(h)
        }
      }
    }
    maxdens <-
      function(g, h) {
        ord = length(g)
        if (compare::compare(g, h)$result == TRUE) {
          return(g) # if identical take any of them as min
        }
        else{
          for (i in 1:(ord + 1)) {
            if (g[i] < h[i])
              return(h)
            if (g[i] > h[i])
              return(g)
          }
        }
      }
    
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
    
    ####### find Nash equilibrium of aux. game -> use single-goal algo for FP (iterate part of opponent d times)
    #### initialization
    x <- rep(0, n)
    y <- list()
    for (p in 1:d) {
      y[[p]] <- rep(0, m)
      
    }
    
    ##### first step for each player
    ### defender
    ## max density of each row
    smax <-
      rep(list(), n) # each density stored in a list (length ord+1)
    for (i in 1:n) {
      smax[[i]] <- A[[i]][[1]]
      for (j in mRange) {
        smax[[i]] <- maxdens(smax[[i]], A[[i]][[j]])
      }
    }
    # min of row-minimas
    vlow <- smax[[1]]
    
    for (i in nRange) {
      vlow <- mindens(vlow, smax[[i]])
    }
    # position of min
    vlowpos <- c(rep(FALSE, n))
    for (i in 1:n) {
      vlowpos[i] <- identical(smax[[i]], vlow)
    }
    row <- min(which(vlowpos))
    U <- rep(list(), n)
    
    for (i in 1:n) {
      U[[i]] <- A[[1]][[1]] # same format
    }
    for (i in 1:(ord + 1)) {
      for (j in 1:n) {
        U[[j]][[i]] <- 0
      }
    }
    
    ##### first step for each opponent ##############
    ### p-th player
    ## min density of each column
    for (p in 1:d) {
      smin <- Ai[[p]][[1]]
      for (j in 1:m) {

        for (i in nRange) {
          smin[[j]] <- mindens(smin[[j]], Ai[[p]][[i]][[j]])
        }
      }
      # max of col-minimas
      vup <- smin[[1]]
      
      for (j in mRange) {
        vup <- maxdens(vup, smin[[j]])
      }
      # position of max
      vuppos <- c(rep(FALSE, m))
      for (j in 1:m) {
        vuppos[j] <- identical(smin[[j]], vup)
      }
      col <-
        min(which(vuppos)) # which density equals vup (i.e. where is result of compare TRUE for the first time?)
      for (i in 1:n) {
        U[[i]] <- U[[i]] + w[p] * Ai[[p]][[i]][[col]]
      }
      y[[p]][[col]] <- y[[p]][[col]] + 1
    }
    #################################################
    
    ### initialization for opponents
    V <- list()
    for (p in 1:d) {
      V[[p]] <- A[[1]]
      for (j in 1:m) {
        V[[p]][[j]] <- rep(0, (ord + 1))
      }
    }
    
    ################### play game until time T
    if (missing(eps)) {
      eps <- -1
    }
    delta <- Inf
    itnumb <- 0
    while (delta > eps && itnumb < T) {
      itnumb <- itnumb + 1
      
      # min of U
      Umin <- U[[1]]
      
      for (i in nRange) {
        Umin <- mindens(Umin, U[[i]])
      }
      # position of min
      Uminpos <- c(rep(FALSE, n))
      for (i in 1:n) {
        Uminpos[i] <- compare::compare(U[[i]], Umin)$result
      }
      row <- min(which(Uminpos, TRUE))
      Uminscal <- Umin / itnumb # rescale densities
      vup <- maxdens(Uminscal, vup)
      
      xbefore <- x
      x[row] <- x[row] + 1
      delta <- max(abs(xbefore - x)) / itnumb
      
      ### optimize goal p 
      for (p in 1:d) {
        for (j in 1:m) {
          V[[p]][[j]] <- V[[p]][[j]] + w[p] * Ai[[p]][[row]][[j]]
        }
        Vmax <- V[[p]][[1]]
        
        for (j in mRange) {
          Vmax <- maxdens(Vmax, V[[p]][[j]])
        }
        # position of max
        Vmaxpos <- c(rep(FALSE, m))
        for (j in 1:m) {
          Vmaxpos[j] <- compare::compare(V[[p]][[j]], Vmax)$result
        }
        col <- min(which(Vmaxpos, TRUE))
        Vmaxscal <- Vmax / itnumb # rescale densities
        vlow <- mindens(Vmaxscal, vlow)
        for (i in 1:n) {
          U[[i]] <- U[[i]] + w[p] * Ai[[p]][[i]][[col]]
        }
        y[[p]][[col]] <- y[[p]][[col]] + 1 # end optimization goal p
      }
      
      
    } # end iteration over itnumb
    
    ## Compile result object *********************
    
    ## estimated probabilities for defender (mixed strategies)
    ## compile equilibrium object for further use (with summary, print, etc.)
    equilibrium <- NULL
    optimalDefense <- as.matrix(x / sum(x))
    colnames(optimalDefense) <- "prob."
    rownames(optimalDefense) <- G$defensesDescriptions
    optimalAttacks <-
      matrix(unlist(y) / (itnumb + 1), ncol = d, byrow = F)
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
    
    class(equilibrium) <- "mosg.equilibrium"
    equilibrium
  }
