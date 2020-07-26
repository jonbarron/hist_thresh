csum <- function(z) cumsum(z[1:(length(z) - 1)])

dsum <- function(z) rev(cumsum(rev(z)))[2:length(z)]

argmax <- function(x, f) mean(x[1:(length(x) - 1)][f == max(f)])

clip <- function(z) pmax(1e-30, z)

preliminaries <- function(n, x = NULL) {
  stopifnot(all(n >= 0))
  if (is.null(x)) {
    x <- seq(0, length(n) - 1)
  }
  stopifnot(all(x[2:length(x)] >= x[1:(length(x) - 1)]))
  w0 <- clip(csum(n))
  w1 <- clip(dsum(n))
  p0 <- w0 / (w0 + w1)
  p1 <- w1 / (w0 + w1)
  mu0 <- csum(n * x) / w0
  mu1 <- dsum(n * x) / w1
  d0 <- csum(n * x ^ 2) - w0 * mu0 ^ 2
  d1 <- dsum(n * x ^ 2) - w1 * mu1 ^ 2
  list(x, w0, w1, p0, p1, mu0, mu1, d0, d1)
}

GHT <- function(n, x = NULL, nu = 0, tau = 0, kappa = 0, omega = 0.5) {
  stopifnot(nu >= 0)
  stopifnot(tau >= 0)
  stopifnot(kappa >= 0)
  stopifnot(omega >= 0 && omega <= 1)
  p <- preliminaries(n, x)
  x <- p[[1]]
  w0 <- p[[2]]
  w1 <- p[[3]]
  p0 <- p[[4]]
  p1 <- p[[5]]
  d0 <- p[[8]]
  d1 <- p[[9]]
  v0 <- clip((p0 * nu * tau ^ 2 + d0) / (p0 * nu + w0))
  v1 <- clip((p1 * nu * tau ^ 2 + d1) / (p1 * nu + w1))
  f0 <- -d0 / v0 - w0 * log(v0) + 2 * (w0 + kappa *      omega)  * log(w0)
  f1 <- -d1 / v1 - w1 * log(v1) + 2 * (w1 + kappa * (1 - omega)) * log(w1)
  list(
    argmax(x, f0 + f1),
    f0 + f1
  )
}
