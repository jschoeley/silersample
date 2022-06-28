# Sample from Siler Distribution

Jonas Schöley
June 28, 2022

The Siler Distribution is characterized by a hazard of the form
_h_(_x_) = *a*<sub>1</sub>exp (−*b*<sub>1</sub>_x_) + *a*<sub>2</sub> + *a*<sub>3</sub>exp (_b_<sub>3</sub>_x_),
i.e. a competing risks model with 3 independent “causes” of death. Thus,
one can sample survival times from the model via
min (_X_<sub>1</sub>,_X_<sub>2</sub>,_X_<sub>3</sub>), where
_X_<sub>_i_</sub> are the latent survival times sampled for each of the
three causes of death. The lowest of those times is the time of death.
The hazard implies that survival times are distributed as negative
Gompertz for cause 1, as exponential for cause 2, and as Gompertz for
cause 3. Survival times from these distributions can readily be using
`base` R and the `flexsurv` package.

```r
rsiler <- function (
    n = 100, a1 = 0.001, b1 = -0.2, a2 = 0.0003, a3 = 0.0001, b3 = 0.14
) {
  require(flexsurv)
  X1 <- rgompertz(n = n, rate = a1, shape = b1)
  X2 <- rexp(n = n, rate = a2)
  X3 <- rgompertz(n = n, rate = a3, shape = b3)
  X <- pmin(X1, X2, X3)
  return(X)
}
```

We implement the Siler density function to graphically check if the
sample converges to the true density.

```r
dsiler <- function (
    x, a1 = 0.001, b1 = -0.2, a2 = 0.0003, a3 = 0.0001, b3 = 0.14
) {
  require(flexsurv)
  # survival function of Siler is S1(x)*S2(x)*S3(x)
  S1 <- pgompertz(x, rate = a1, shape = b1, lower.tail = FALSE)
  S2 <- pexp(x, rate = a2, lower.tail = FALSE)
  S3 <- pgompertz(x, rate = a3, shape = b3, lower.tail = FALSE)
  S <- S1*S2*S3
  # Siler hazard
  h <- a1*exp(b1*x) + a2 + a3*exp(b3*x)
  # Siler density
  d = h*S
  return(d)
}

hist(rsiler(n = 10000), freq = FALSE)
lines(dsiler(0:70))
```

![](./ass/unnamed-chunk-2-1.png)<!-- -->
