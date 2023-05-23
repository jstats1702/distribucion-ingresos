# settings
rm(list=ls())

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
set.seed(5)
kappa <- 3
mu    <- rnorm(n = 1, mean = mu0, sd = sqrt(g20))
tau2  <- 1/rgamma(n = 1, shape = 0.5*eta0, rate = 0.5*eta0*t20)
theta <- rnorm(n = m, mean = mu, sd = sqrt(tau2))
sig2  <- rgamma(n = 1, shape = 0.5*al0, rate = 0.5*be0)
vsi2  <- 1/rgamma(n = n, shape = 0.5*kappa, rate = 0.5*kappa*sig2)
# almacenamiento
THETA <- NULL
SIG2  <- NULL
MU    <- NULL
TAU2  <- NULL
LL    <- NULL
# cadena
B <- 11000
n_burn <- 1000
set.seed(5)
for (b in 1:B) {
  # actualizar theta
  for (j in 1:m) {
    vtheta   <- 1/(1/tau2 + sum(1/vsi2[g==j]))
    theta[j] <- rnorm(n = 1, mean = vtheta*(mu/tau2 + sum(y[g==j]/vsi2[g==j])), sd = sqrt(vtheta)) 
  }
  # actualizar sig2
  sig2 <- rgamma(n = 1, shape = 0.5*(al0 + kappa*n), rate = 0.5*(be0 + kappa*sum(1/vsi2)))
  # actualizar vsi2
  vsi2 <- 1/rgamma(n = n, shape = 0.5*(kappa + 1), rate = 0.5*(kappa*sig2 + (y - rep(theta, nj))^2))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
  # actualizar tau^2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # almacenar & log-verosimilitud
  if (b > n_burn) {
    THETA <- rbind(THETA, theta)
    SIG2  <- rbind(SIG2, sig2)
    MU    <- rbind(MU, mu)
    TAU2  <- rbind(TAU2, tau2)
    LL    <- rbind(LL, sum(metRology::dt.scaled(x = y, df = kappa, mean = rep(theta, nj), sd = sqrt(sig2), log = T)))
  }
}
M5 <- list(THETA = THETA, SIG2 = SIG2, MU = MU, TAU2 = TAU2, LL = LL)
save(M5, file = "Peronas-M5.RData")