# settings
rm(list=ls())
setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/caso2")
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

# modelos ajustados
load("Peronas-M1.RData")
load("Peronas-M2.RData")
load("Peronas-M3.RData")
load("Peronas-M4.RData")
load("Peronas-M5.RData")
load("Peronas-M6.RData")

# log-verosimilitud
windows()
yrange <- range(M1$LL, M2$LL, M3$LL, M4$LL, M5$LL, M6$LL) + c(0,10)
plot (M1$LL, type = "p", pch = ".", cex = 1.1, col = 1, ylim = yrange, xlab = "Iteración", ylab = "Log-verosimilitud", main = "")
lines(M2$LL, type = "p", pch = ".", cex = 1.1, col = 2)
lines(M3$LL, type = "p", pch = ".", cex = 1.1, col = 3)
lines(M4$LL, type = "p", pch = ".", cex = 1.1, col = 4)
lines(M5$LL, type = "p", pch = ".", cex = 1.1, col = 5)
lines(M6$LL, type = "p", pch = ".", cex = 1.1, col = 6)
legend("top", legend = paste0("M",1:6), fill = 1:6, border = 1:6, bty = "n", horiz = T)

# DIC
lpyth_m1 <- sum(dnorm(x = y, mean = mean(M1$THETA), sd = sqrt(mean(M1$SIG2)), log = T))
lpyth_m2 <- sum(dnorm(x = y, mean = rep(colMeans(M2$THETA), nj), sd = sqrt(mean(M2$SIG2)), log = T))
lpyth_m3 <- sum(dnorm(x = y, mean = rep(colMeans(M3$THETA), nj), sd = sqrt(rep(colMeans(M3$SIG2), nj)), log = T))
lpyth_m4 <- sum(metRology::dt.scaled(x = y, df = 3, mean = mean(M4$THETA), sd = sqrt(mean(M4$SIG2)), log = T))
lpyth_m5 <- sum(metRology::dt.scaled(x = y, df = 3, mean = rep(colMeans(M5$THETA), nj), sd = sqrt(mean(M5$SIG2)), log = T))
lpyth_m6 <- sum(metRology::dt.scaled(x = y, df = 3, mean = rep(colMeans(M6$THETA), nj), sd = sqrt(rep(colMeans(M6$SIG2), nj)), log = T))

pDIC_m1 <- 2*(lpyth_m1 - mean(M1$LL))
pDIC_m2 <- 2*(lpyth_m2 - mean(M2$LL))
pDIC_m3 <- 2*(lpyth_m3 - mean(M3$LL))
pDIC_m4 <- 2*(lpyth_m4 - mean(M4$LL))
pDIC_m5 <- 2*(lpyth_m5 - mean(M5$LL))
pDIC_m6 <- 2*(lpyth_m6 - mean(M6$LL))

dic_m1 <- -2*lpyth_m1 + 2*pDIC_m1 
dic_m2 <- -2*lpyth_m2 + 2*pDIC_m2 
dic_m3 <- -2*lpyth_m3 + 2*pDIC_m3 
dic_m4 <- -2*lpyth_m4 + 2*pDIC_m4 
dic_m5 <- -2*lpyth_m5 + 2*pDIC_m5 
dic_m6 <- -2*lpyth_m6 + 2*pDIC_m6

tab <- cbind(c(pDIC_m1,pDIC_m2,pDIC_m3,pDIC_m4,pDIC_m5,pDIC_m6),
             c( dic_m1, dic_m2, dic_m3, dic_m4, dic_m5, dic_m6))
rownames(tab) <- paste0("M", 1:6)
colnames(tab) <- c("pDIC","DIC")
round(tab, 1)

# estadísticos de prueba: Bogotá
tsf <- function(x) as.numeric(c(mean(x), median(x), sd(x), sd(x)/mean(x), diff(range(x)), diff(quantile(x,c(0.25,0.75)))))
i   <- 3
yi  <- y[g==i]
ni  <- length(yi)
# estadísticos de prueba
B <- 10000
TS <- array(data = NA, dim = c(B,6,6))
set.seed(123456)
for (b in 1:B) {
  # M1
  yrep <- rnorm(n = ni, mean = M1$THETA[b], sd = sqrt(M1$SIG2[b]))
  TS[b,1,] <- tsf(yrep)
  # M2
  yrep <- rnorm(n = ni, mean = M2$THETA[b,i], sd = sqrt(M2$SIG2[b]))
  TS[b,2,] <- tsf(yrep)
  # M3
  yrep <- rnorm(n = ni, mean = M3$THETA[b,i], sd = sqrt(M3$SIG2[b,i]))
  TS[b,3,] <- tsf(yrep)
  # M4
  yrep <- metRology::rt.scaled(n = ni, df = 3, mean = M4$THETA[b], sd = sqrt(M4$SIG2[b]))
  TS[b,4,] <- tsf(yrep)
  # M5
  yrep <- metRology::rt.scaled(n = ni, df = 3, mean = M5$THETA[b,i], sd = sqrt(M5$SIG2[b]))
  TS[b,5,] <- tsf(yrep)
  # M6
  yrep <- metRology::rt.scaled(n = ni, df = 3, mean = M6$THETA[b,i], sd = sqrt(M6$SIG2[b,i]))
  TS[b,6,] <- tsf(yrep)
}
# tabla
TS0 <- tsf(yi)
PPP <- matrix(NA, nrow = 6, ncol = 6)
for (j in 1:6) {
  for (k in 1:6) {
    PPP[j, k] <- mean(TS[,j,k] > TS0[k])
  }
}
rownames(PPP) <- paste0("M_", 1:6)
colnames(PPP) <- c("Media","Mediana","DE","CV","R","RI")
round(PPP,2)
# viz
TS_display <- c("Media","Mediana","DE","CV","R","RI")
windows(width = 7.5, height = 5)
par(mfrow = c(2,3), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
for (k in 1:6) {
  hist(x = TS[,6,k], freq = F, nclass = 20, col = "gray90", border = "gray90", ylab = "Densidad", xlab = TS_display[k], main = TS_display[k])
  abline(v = quantile(TS[,6,k], c(0.025,0.975)), lty = 2)
  abline(v = TS0[k], col = 2)
}

# ranking Bayesiano
ids  <- estadisticos$dominio
that <- colMeans(M6$THETA)
ic   <- apply(X = M6$THETA, MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
# ranking
ranking <- order(that) 
ids  <- ids [ ranking]
that <- that[ ranking]
ic   <- ic  [,ranking]
# colores
c0 <- log(1014980) # SMLMV
colo <- rep(2, m)
colo[which(ic[1,] > c0)] <- 3
colo[which(ic[2,] < c0)] <- 1
colo <- c("darkred","black","darkgreen")[colo]
# viz
windows()
par(mfrow = c(1,1), mar = c(4,8,1.5,1.5), mgp = c(2.5,0.75,0))
plot(NA, NA, xlab = "Ingresos (log)", ylab = "", main = "Ranking", xlim = range(ic), ylim = c(1,m), cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:m, labels = ids, las = 2)
abline(v = c0,  col = "gray", lwd = 3)
abline(h = 1:m, col = "lightgray", lwd = 1)
for (j in 1:m) {
  segments(x0 = ic[1,j], y0 = j, x1 = ic[2,j], y1 = j, col = colo[j])
  lines(x = that[j], y = j, type = "p", pch = 16, cex = 0.8, col = colo[j])
}

# top 5
ids   <- estadisticos$dominio
that  <- colMeans(M6$THETA)
s2hat <- colMeans(3/(3-2)*M6$SIG2)
ranking <- order(that) 
ids   <- ids [ranking]
that  <- that[ranking]
s2hat <- s2hat[ranking]
tab <- cbind(25:1,exp(that),sqrt(exp(that)^2*s2hat), sqrt(s2hat)*100)
colnames(tab) <- c("Ranking","Media","DE","CV(%)")
rownames(tab) <- ids
round(tab[25:21,],1)

# agrupamiento
# https://uc-r.github.io/kmeans_clustering
library(factoextra)
library(cluster)
X <- scale(x = cbind(colMeans(M6$THETA), colMeans(sqrt(M6$SIG2))))
rownames(X) <- dominios
colnames(X) <- c("theta","sigma")
# agrupamiento jerárquico
clust <- agnes(x = X, method = "ward")
labs <- cutree(as.hclust(clust), k = 4)
colo <- c("blue","green","red","black")
# viz
windows(width = 7.5, height = 5)
par(mfrow = c(1,2), mar = c(4,3,1.4,1.4), mgp = c(1.75,0.75,0))
pltree(clust, cex = 0.6, hang = -1, main = "") 
rect.hclust(clust, k = 4, border = colo[c(1,4,2,3)])
plot(X, col = colo[labs], pch = c(15:18)[labs], cex = 1.1, xlab = expression(theta), ylab = expression(sigma))
