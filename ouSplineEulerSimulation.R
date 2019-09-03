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

ouSplineEulerSimulation <- function(x0, gamma, spline, breaks, sigma=1, t) {
        # approximate Orstein-Uhlenbeck process
        # dX = \gamma (\mu - X) dt + \sigma dW
        # by Euler-Maruyama approximation
        # t is time discretisation in [t0, T]
        # t0 = t_0 < t_1 < ... < t_N = T
        n <- length(breaks) - 1
        numBreaks <- nrow(spline)
        maxBreak <- breaks[numBreaks + 1]
        N <- length(t) - 1 
        y <- numeric(length = (N + 1))
        y[1] <- x0
        deltaT <- t[2:(N+1)] - t[1:N]
        deltaW <- rnorm(N, mean = 0, sd = sqrt(deltaT))
        # split t up into relevant breaks
        breakMultiplier <- max(t) / maxBreak
        fullBreaks1 <- rep(breaks, length = (length(breaks) + 1) * breakMultiplier)
        fullBreak2 <- rep(0:ceiling(breakMultiplier), each = length(breaks), 
                          length = (length(breaks) + 1) * breakMultiplier) * 
                                (maxBreak + 1)
        fullBreak <- fullBreaks1 + fullBreak2 
        whichTBreaks <- lapply(seq_len(length(fullBreak) - 1), function(i) {
                                  bF <- fullBreak[i+1]
                                  bB <- fullBreak[i]
                                  #tMod <- t %% maxBreak
                                  whichT <- which(t < bF & t >= bB)
                                  whichT
                        })
        degree <- c(2, 1, 0)
        derivativeDegree <- c(1, 0)
        muTBreaks <- lapply(seq_len(length(fullBreak) - 1), function(i) {
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
        muT <- unlist(lapply(seq_len(length(fullBreak) - 1), function(i) {
                                     muTBreaks[[i]]$mu
                        }))
        muDerT <- unlist(lapply(seq_len(length(fullBreak) - 1), function(i) {
                                     muTBreaks[[i]]$muDerivative
                        }))
        ###
        for(j in 1:N) {
                termOne <- -gamma * (y[j] - muT[j] - muDerT[j] / gamma) * deltaT[j]
                termTwo <- sigma * deltaW[j]
                y[j + 1] <- y[j] + termOne + termTwo
        }
        list(t = t, y = y, mu = muT)
}

if(FALSE) {
        # model splines (for plotting)
        lenSeq <- 100
        xy <- lapply(seq_len(n), function(i) {
                             pTmp <- breaks[i+1]
                             pBackTmp <- breaks[i]
                        tmpSeq <- seq(from = pBackTmp, to = pTmp, length = lenSeq)
                        tmpSeq <- tmpSeq[-lenSeq] 
                        tmpSeqMod <- tmpSeq %% maxBreak
                        xQuadratic <- outer(tmpSeqMod, degree, "^")
                        iMod <- i %% maxBreak
                        if (iMod == 0) iMod <- maxBreak
                        splineTmp <- spline[iMod,]
                        yTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                        #yTmp
                        list(x = tmpSeq, y = yTmp)
                                             })
        xSpline <- unlist(lapply(seq_len(n), function(i) xy[[i]]$x))
        ySpline <- unlist(lapply(seq_len(n), function(i) xy[[i]]$y))
}
