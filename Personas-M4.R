# settings
rm(list=ls())

# datos
load("Personas.RData")

# datos y estadísticos
y  <- log(dat$ingtot)
n  <- length(y)
yb <- mean(y)

# hiperparámetros
mu0 <- mean(y)
t20 <- (0.25*mu0)^2
al0 <- 1
be0 <- 1/var(y)
# valores iniciales
set.seed(4)
kappa <- 3
theta <- rnorm(n = 1, mean = mu0, sd = sqrt(t20))
sig2  <- rgamma(n = 1, shape = 0.5*al0, rate = 0.5*be0)
vsi2  <- 1/rgamma(n = n, shape = 0.5*kappa, rate = 0.5*kappa*sig2)
# almacenamiento
THETA <- NULL
SIG2  <- NULL
LL    <- NULL
# cadena
B <- 11000
n_burn <- 1000
set.seed(4)
for (b in 1:B) {
  # actualizar theta
  vtheta <- 1/(1/t20 + sum(1/vsi2))
  theta  <- rnorm(n = 1, mean = vtheta*(mu0/t20 + sum(y/vsi2)), sd = sqrt(vtheta)) 
  # actualizar sigma^2
  sig2 <- rgamma(n = 1, shape = 0.5*(al0 + kappa*n), rate = 0.5*(be0 + kappa*sum(1/vsi2)))
  # actualizar upsilon^2
  vsi2 <- 1/rgamma(n = n, shape = 0.5*(kappa + 1), rate = 0.5*(kappa*sig2 + (y - theta)^2))
  # almacenar & log-verosimilitud
  if (b > n_burn) {
    THETA <- rbind(THETA, theta)
    SIG2  <- rbind(SIG2, sig2)
    LL    <- rbind(LL, sum(metRology::dt.scaled(x = y, df = kappa, mean = theta, sd = sqrt(sig2), log = T)))
  }
}
M4 <- list(THETA = THETA, SIG2 = SIG2, LL = LL)
save(M4, file = "Peronas-M4.RData")