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
nu0  <- 1  
s20  <- var(y)
# valores iniciales
theta <- yb
sig2  <- mean(s2)
mu    <- mean(theta)
tau2  <- var (theta)
# almacenamiento
THETA <- NULL
SIG2  <- NULL
MU    <- NULL
TAU2  <- NULL
LL    <- NULL
# cadena
B <- 11000
n_burn <- 1000
set.seed(2)
for (b in 1:B) {
  # actualizar theta
  vtheta <- 1/(1/tau2 + nj/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
  # actualizar sigma^2
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0 + n), rate = 0.5*(nu0*s20 + sum((nj-1)*s2 + nj*(yb - theta)^2)))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
  # actualizar tau^2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # almacenar valores & log-verosimilitud
  if (b > n_burn) {
    THETA <- rbind(THETA, theta)
    SIG2  <- rbind(SIG2, sig2)
    MU    <- rbind(MU, mu)
    TAU2  <- rbind(TAU2, tau2)
    LL    <- rbind(LL, sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(sig2), log = T)))
  }
}
M2 <- list(THETA = THETA, SIG2 = SIG2, MU = MU, TAU2 = TAU2, LL = LL)
save(M2, file = "Peronas-M2.RData")