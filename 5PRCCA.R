library(epiR)
library(deSolve)

setwd("C:/Users/simon/Documents/Doutorado/1_Tese/5caravan/5caravan")
source("5MODEL_PRCCA.R")

# Assuming uniform distribution
set.seed(1234)

death <- runif(n=100, min=0, max=0.01/365)
w1 <- runif(n=100, min=0.1/365, max=0.95/365)
w2 <- runif(n=100, min=0.1/365, max=0.95/365)
w3 <- runif(n=100, min=0.1/365, max=0.95/365)
mum <- runif(n=100, min=0.4/30, max=1/30)
muh <- runif(n=100, min=0.001/30, max=0.0015/30)
a <- runif(n=100, min=5/30, max=10/30)
b <- runif(n=100, min=0.2, max=0.5)
c <- runif(n=100, min=0.2, max=0.5)
epsilonm <- runif(n=100, min=1/10, max=1/4)
epsilonh <- runif(n=100, min=1/10, max=1/4)
gamah <- runif(n=100, min=1/10, max=1/4)
alfah <- runif(n=100, min=0.01/30, max=0.05/30)

# Defining a vector for Infh1 #----
Infh1 <- numeric(length(a))

for (i in 1:length(a)){
  prcca1 <- model1(a[i],b[i],c[i],muh[i],mum[i],gamah[i],epsilonh[i],
                 epsilonm[i],alfah[i],w1[i])
  Infh1[i] <- prcca1$Infh1
}

# Sensitivity Analysis PRCC for Infh1 ----
dat1 <- data.frame(cbind(X1=a, X2=b, X3=c, X4=muh, X5=mum, X6=gamah,
                         X7=epsilonh, X8=epsilonm, X9=alfah, X10=w1, Y=Infh1))
colnames(dat1) <- c("biting rate","Fraction of infective bites",
                    "mosquito susceptibility","Human natural mortality rate",
                    "mosquito mortality rate",
                    "Human recovery rate", "latency rate in humans", 
                    "latency rate in mosquitoes", 
                    "Disease-induced mortality rate", 
                    "vaccination rate in the endemic community", "cases")
est1 <- epi.prcc(dat1, sided.test = 2, conf.level = 0.95)
est1


# PLOTTING PRCCA 1----

par(mfrow=c(1,1))

bp1 <- barplot(est1$est, ylim=c(-1,1), main = "Endemic community",
              names.arg=c("a", "b", "c",expression(mu[h]), expression(mu[m]), 
                          expression(gamma[h]),expression(epsilon[h]),
                          expression(epsilon[m]), expression(alpha[h]), "w1"),
              ylab='Coefficient', axes=FALSE)
axis(2)

mtext(text='Parameters', side=1, line=2.5)
box()
#for(i in 1:8) lines(c(i,i),c(summary[i,4], summary[i,5])) #the lines for CI
abline(h=0)
text(bp1, -0.9, paste("r = ", round(est1$est,2), sep=""),cex=0.5, pos=3) 
grid(nx = 11, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.8)      # Grid line width


# Sensitivity Analysis PRCC for Infh2 ----
source("5MODEL_PRCCA.R")

Infh2 <- numeric(length(a))

for (i in 1:length(a)){
  prcca2 <- model2(a[i],b[i],c[i],muh[i],mum[i],gamah[i],epsilonh[i],
                   epsilonm[i],alfah[i],w2[i], death[i])
  Infh2[i] <- prcca2$Infh2
}

# Sensitivity Analysis PRCC for Infh1 ----
dat2 <- data.frame(cbind(X1=a, X2=b, X3=c, X4=muh, X5=mum, X6=gamah,
                         X7=epsilonh, X8=epsilonm, X9=alfah, X10=w2, 
                         X11=death, Y=Infh2))
colnames(dat2) <- c("biting rate","Fraction of infective bites",
                    "mosquito susceptibility","Human natural mortality rate",
                    "mosquito mortality rate",
                    "Human recovery rate", "latency rate in humans", 
                    "latency rate in mosquitoes", 
                    "Disease-induced mortality rate", 
                    "vaccination rate in the caravan", "additional mortality", 
                    "cases")
est2 <- epi.prcc(dat2, sided.test = 2, conf.level = 0.95)
est2

# PLOTTING 2 ----
par(mfrow=c(1,1))

bp2 <- barplot(est2$est, ylim=c(-1,1), main = "Caravan",
               names.arg=c("a", "b", "c",expression(mu[h]), expression(mu[m]), 
                           expression(gamma[h]),expression(epsilon[h]),
                           expression(epsilon[m]), expression(alpha[h]), "w2",
                           expression(theta[c])),
               ylab='Coefficient', axes=FALSE)
axis(2)

mtext(text='Parameters', side=1, line=2.5)
box()
#for(i in 1:8) lines(c(i,i),c(summary[i,4], summary[i,5])) #the lines for CI
abline(h=0)
text(bp2, -0.9, paste("r = ", round(est2$est,2), sep=""),cex=0.5, pos=3) 
grid(nx = 12, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.8)      # Grid line width


# Sensitivity Analysis PRCC for Infh3 ----

source("5MODEL_PRCCA.R")

# Defining a vector for Infh3 #
Infh3 <- numeric(length(a))

for (i in 1:length(a)){
  prcca3 <- model3(a[i],b[i],c[i],muh[i],mum[i],gamah[i],epsilonh[i],
                 epsilonm[i],alfah[i],w3[i])
  Infh3[i] <- prcca3$Infh3
}


dat3 <- data.frame(cbind(X1=a, X2=b, X3=c, X4=muh, X5=mum, X6=gamah,
                         X7=epsilonh, X8=epsilonm, X9=alfah, X10=w3, Y = Infh3))
colnames(dat3) <- c("biting rate","Fraction of infective bites",
                    "mosquito susceptibility",
                    "Human natural mortality rate","mosquito mortality rate",
                    "Human recovery rate", "latency rate in humans", "latency rate in mosquitoes", 
                    "Disease-induced mortality rate",
                    "vaccination rate in the disease-free community", "cases")
est3 <- epi.prcc(dat3, sided.test = 2, conf.level = 0.95)
est3

# PLOTTING PRCCA 3----

par(mfrow=c(1,1))

bp3 <- barplot(est3$est, ylim=c(-1,1), main="Disease-free community", 
              names.arg=c("a", "b", "c", expression(mu[h]),
                          expression(mu[m]), expression(gamma[h]),expression(epsilon[h]),
                          expression(epsilon[m]), expression(alpha[h]), "w3"),
              ylab='Coefficient', axes=FALSE)
axis(2)

mtext(text='Parameters', side=1, line=2.5)
box()
#for(i in 1:8) lines(c(i,i),c(summary[i,4], summary[i,5])) #the lines for CI
abline(h=0)
text(bp3, -0.9, paste("r = ", round(est3$est,2), sep=""),cex=0.5, pos=3) 
grid(nx = 11, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.8)      # Grid line width


# PLOTTING ALL 3 PLOTS TOGETHER
par(mfrow=c(1,3))
# go to bp1 and bp3 again


#stop here for now, keep going later


#### UNCERTANTY ANALYSIS OF b (MONTE CARLO SIMULATION) ----
source("PRCCA_model.R")

b <- runif(n=100, min=0.2, max=0.5)
Infh2 <- numeric(length(b))

for (i in 1:length(b)){
  Infh2[i] <- inf.bites(b[i])$Infh2
}

par(mfrow=c(1,2))
plot(b, Infh2,
     xlab= "parameter b",
     ylab= "Cases in host community")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width


# Probability interval (95% PI) 
quantile(Infh2, probs=c(0.025,0.975))


#### UNCERTANTY ANALYSIS OF MIGRATION RATE deltah_12 (MONTE CARLO SIMULATION) ----
source("PRCCA_model.R")

deltah_12 <- runif(n=100, min=0.1/365, max=0.9/365)
Infh2 <- numeric(length(deltah_12))

for (i in 1:length(deltah_12)){
  Infh2[i] <- delta(deltah_12[i])$Infh2
}

plot(deltah_12, Infh2,
     xlab= "migration rate",
     ylab= "Cases in host community")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width

# Probability interval (95% PI) 
quantile(Infh2, probs=c(0.025,0.975))