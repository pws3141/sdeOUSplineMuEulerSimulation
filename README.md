# Simulate SDE using Euler–Maruyama method

SDE given by:
\[
        dX_t = -\gamma (X_t - m(t) - \frac{m'(t)}{gamma}) dt +
                    \sigma dB_t
                    \]

Here, m(t) is a set of quadratic functions found using splines over the partitioned
space. See [this repo](https://github.com/pws3141/splineEstimation) for the function
to estimate quadratic periodic splines given the means of data in each partition.
