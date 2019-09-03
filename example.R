#! /usr/bin/env Rscript

# ouSplines R function required: see splineFunction.R avaliable at
# https://github.com/pws3141/splineEstimation
# source("../sdeOUSplineMu/splineFunction.R")

# example data is mean of montly data over N years
exampleMeans <- c(-42.007931, -50.501747, -47.866878, -32.840013,  -8.159889,  21.157662,
                        47.859793, 59.509856, 52.442058, 28.548208, -1.298951, -26.842167)

# partition includes both end points
partition <- (0:12) 
spline <- ouSplines(means = exampleMeans, breaks = partition)
#checks <- splineCheck(optPar = spline$par, means = exampleMeans, breaks = partition)

# plot
# splinePlot(spline = spline, r = 6, savePlot = TRUE)

# simulate
tau <- seq(from = 0, to = 120, length = 6000)
ouSim <- ouSplineEulerSimulation(x0 = exampleMeans[1] + 20, gamma = 1, 
                                 spline = spline$par, breaks = partition, 
                                 sigma = 10, t = tau)                                              

#plotRange <- range(c(ouSim$y, ouSim$mu))
#plot(ouSim$t, ouSim$y, t = 'l', lwd = 2, ylim = plotRange)
#lines(ouSim$t, ouSim$mu, col = rgb(0, 0, 0, .5))

