# settings
rm(list=ls())

# datos
load("Personas.RData")

# datos y estadísticos
y  <- log(dat$ingtot)
n  <- length(y)
yb <- mean(y)
s2 <- var(y)

# hiperparámetros
mu0 <- mean(y)
t20 <- (0.25*mu0)^2
nu0 <- 1
s20 <- var(y)
# valores iniciales
set.seed(1)
theta <- rnorm(n = 1, mean = mu0, sd = sqrt(t20))
sig2  <- 1/rgamma(n = 1, shape = 0.5*nu0, rate = 0.5*nu0*s20 )
# almacenamiento
THETA <- NULL
SIG2  <- NULL
LL    <- NULL
# cadena
B <- 11000
n_burn <- 1000
set.seed(1)
for(b in 1:B) {
  # actualizar theta
  v2 <- 1/(1/t20 + n/sig2)      
  theta <- rnorm(n = 1, mean = v2*(mu0/t20 + n*yb/sig2), sd = sqrt(v2))
  # actualizar sig2
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0 + n), rate = 0.5*(nu0*s20 + (n-1)*s2 + n*(yb - theta)^2))
  # almacenar & log-verosimilitud
  if (b > n_burn) {
    THETA <- rbind(THETA, theta)
    SIG2  <- rbind(SIG2, sig2)
    LL    <- rbind(LL, sum(dnorm(x = y, mean = theta, sd = sqrt(sig2), log = T)))
  }
}
M1 <- list(THETA = THETA, SIG2 = SIG2, LL = LL)
save(M1, file = "Peronas-M1.RData")