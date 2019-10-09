# Aim:
# Simulate SDE using Euler approximation.

# Input:
# - breaks
# - spline (a, b, c) values for each break
# - gamma

# Output:
# - simulated SDE path
# - associated time
# - associtaed spline for mu

ouSplineEulerSimulation <- function(x0, g, spline, breaks, s=1, t) {
        # approximate Orstein-Uhlenbeck process
        # dX = -\gamma X dt + \sigma dW
        # Y = X + m(t)
        # dY = -gamma (Y - m(t)) dt + m'(t) dt + sigma dW
        # by Euler-Maruyama approximation
        # t is time discretisation in [t0, T]
        # t0 = t_0 < t_1 < ... < t_N = T
        n <- length(breaks) - 1
        numBreaks <- nrow(spline)
        if (n != numBreaks) stop("breaks is not compatible with spline")
        maxBreak <- breaks[numBreaks + 1]
        N <- length(t) - 1 
        y <- numeric(length = (N + 1))
        y[1] <- x0
        maxT <- max(t)
        deltaT <- t[2:(N+1)] - t[1:N]
        deltaW <- rnorm(N, mean = 0, sd = sqrt(deltaT))
        # split t up into relevant breaks
        breakMultiplier <- max(t) / maxBreak
        fullBreaks1 <- rep(breaks, length = (length(breaks) + 1) * breakMultiplier)
        fullBreaks2 <- rep(0:ceiling(breakMultiplier), each = length(breaks), 
                          length = (length(breaks) + 1) * breakMultiplier) * 
                                (maxBreak + 1)
        fullBreaks <- fullBreaks1 + fullBreaks2 
        fullBreaks <- fullBreaks[1:min(which(fullBreaks > maxT))]
        whichTBreaks <- lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                  bF <- fullBreaks[i+1]
                                  bB <- fullBreaks[i]
                                  #tMod <- t %% maxBreak
                                  whichT <- which(t < bF & t >= bB)
                                  whichT
                        })
        degree <- c(2, 1, 0)
        derivativeDegree <- c(1, 0)
        muTBreaks <- lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                  w <- whichTBreaks[[i]]
                                  tmpSeq <- t[w]
                                  tmpSeqMod <- tmpSeq %% maxBreak
                                  xQuadratic <- outer(tmpSeqMod, degree, "^")
                                  xDerivative <- outer(tmpSeqMod,
                                                       derivativeDegree, "^")
                                  xDerivative[, 1] <- 2 * xDerivative[, 1]
                                  iMod <- i %% maxBreak
                                  if (iMod == 0) iMod <- maxBreak
                                  splineTmp <- spline[iMod,]
                                  muTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                                  muDerivativeTmp <- as.vector(tcrossprod(splineTmp[-3], 
                                                                          xDerivative))
                                  #yTmp
                                  list(t = tmpSeq, mu = muTmp,
                                        muDerivative = muDerivativeTmp)
                        })
        muT <- unlist(lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                     muTBreaks[[i]]$mu
                        }))
        muDerT <- unlist(lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                     muTBreaks[[i]]$muDerivative
                        }))
        ###
        for(j in 1:N) {
                termOne <- -g * (y[j] - muT[j] - muDerT[j] / g) * deltaT[j]
                termTwo <- s * deltaW[j]
                y[j + 1] <- y[j] + termOne + termTwo
        }
        list(t = t, y = y, mu = muT)
}
