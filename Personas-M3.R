# settings
rm(list=ls())
suppressMessages(suppressWarnings(library(dplyr))) 

# datos
load("Personas.RData")

# datos y estadísticos
m <- length(table(dat$dominio))
n <- sum(table(dat$dominio))
y <- log(dat$ingtot)
Y <- vector(mode = "list", length = m)
g <- rep(NA, n)
for (j in 1:m) {
  idx <- dat$dominio == unique(dat$dominio)[j]
  g[idx] <- j
  Y[[j]] <- y[idx]
}
estadisticos <- dat %>% 
  group_by(dominio) %>% 
  summarise(dominio = unique(dominio), nj = n(), yb = mean(log(ingtot)), s2 = var(log(ingtot)))
nj <- estadisticos$nj
yb <- estadisticos$yb
s2 <- estadisticos$s2

# hiperparámetros
mu0  <- mean(y)
g20  <- (0.25*mu0)^2
eta0 <- 1  
t20  <- var(y)
al0  <- 1
be0  <- 1/var(y)
# valores iniciales
theta <- yb
sig2  <- s2   # sigma_j^2
mu    <- mean(theta)
tau2  <- var (theta)
nu    <- 1
ups2  <- 100  # sigma^2
# almacenamiento
THETA <- NULL
SIG2  <- NULL
MU    <- NULL
TAU2  <- NULL
UPS2  <- NULL
LL    <- NULL
# cadena
B <- 11000
n_burn <- 1000
set.seed(3)
for (b in 1:B) {
  # actualizar theta
  vtheta <- 1/(1/tau2 + nj/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
  # actualizar sigma_j^2
  sig2 <- 1/rgamma(n = m, shape = 0.5*(nu + nj), rate = 0.5*(nu*ups2 + (nj-1)*s2 + nj*(yb - theta)^2))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu))
  # actualizar tau2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # actualizar sigma^2
  ups2 <- rgamma(n = 1, shape = 0.5*(al0 + m*nu), rate = 0.5*(be0 + nu*sum(1/sig2)))
  # almacenar & log-verosimilitud
  if (b > n_burn) {
    THETA <- rbind(THETA, theta)
    SIG2  <- rbind(SIG2, sig2)
    MU    <- rbind(MU, mu)
    TAU2  <- rbind(TAU2, tau2)
    UPS2  <- rbind(UPS2, ups2)
    LL    <- rbind(LL, sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(rep(sig2, nj)), log = T)))
  }
}
M3 <- list(THETA = THETA, SIG2 = SIG2, MU = MU, TAU2 = TAU2, UPS2 = UPS2, LL = LL)
save(M3, file = "Peronas-M3.RData")