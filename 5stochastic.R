# STOCHASTIC 3 VACCINATION - SCENARIO 2

#change between Sh0 = 1000 to 10000, and also change Ih2 and Ih3

# loading packages ####
library(deSolve)
library(GillespieSSA)
# https://cran.r-project.org/web/packages/GillespieSSA/GillespieSSA.pdf

### 1.1 PARAMETERS ###

a	<- 10/30   # Average biting rate	(Haemagogus and Aedes)
b <- 0.25   # Fraction of bites actually infective (Haem = 0.25 Aedes = 0.2)
c <- 0.4    # Vector susceptibility to the virus (Haemagogus = 0.4, Aedes = 0.25)
death <- 6200/(898000/9)/9/12/30 # Human death rate in caravans	(=3.967e-05/day)
muh	<- 0.00119/30 # Human natural mortality rate	(=3.967e-05/day)
mum <- 0.46/30   # Natural mortality rate of mosquitoes	(Haemagogus= 0.46; Aedes= 0.913/30)
gamah <- 1/7 # Human recovery rate	(4.0 month-1 ~7 days)
epsilonh <- 1/6 # YF latency rate for humans (~6 days)
epsilonm <- 1/7 #	Latency rate in mosquitoes	(4.0 month-1 ~7 days)
alfah <- 0.0244/30 #	Disease-induced mortality rate (0.0244 month-1)

# migration parameters
#deltah <- 1/40 #distance walked daily

# Vaccine parameters
pvac1 <- 0.95 # proportion of effective vaccination coverage
w1 <- pvac1/365 # daily vaccination rate

pvac2 <- 0.95
w2 <- pvac2/365

pvac3 <- 0.95 # proportion of effective vaccination coverage
w3 <- pvac3/365 # daily vaccination rate

par.model <- c(a=a, b=b, c=c, epsilonh=epsilonh, epsilonm=epsilonm, 
               muh=muh, gamah=gamah, alfah=alfah, 
               death=death,mum=mum, w1=w1, w2=w2, w3=w3)


### 1.2. VARIABLES AND INITIAL CONDITIONS ###

Sh0_1 <- 1000; Sm0_1 <- Sh0_1*1.5; m1 <- Sm0_1/Sh0_1
R0_1 <- a*a*b*c*m1*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_1
# caravan
Sh0_2 <- 1000; Sm0_2 <- Sh0_2*1.5; m2 <- Sm0_2/Sh0_2
R0_2 <- a*a*b*c*m2*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_2
# community 3
Sh0_3 <- 1000; Sm0_3 <- Sh0_3*1.5; m3 <- Sm0_3/Sh0_3
R0_3 <- a*a*b*c*m3*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_3


# Initial conditions
state.model <- c(Sh1=Sh0_1-(Sh0_1*pvac1), Eh1=0, Ih1=0, Rh1=0, 
                   Vh1=Sh0_1*pvac1, Sm1=Sm0_1, Em1=0, Im1=1,
                   Sh2=Sh0_2, Eh2=0, Ih2=0.01987244, Rh2=0, Vh2=0,
                   Sm2=Sm0_2, Em2=0, Im2=0,
                   Sh3=Sh0_3, Eh3=0, Ih3=0.005947227, Rh3=0, 
                   Vh3=0, Sm3=Sm0_3, Em3=0, Im3=0) #also change Ih2 and Ih3

### 1.3 STOCHASTIC MODELLING ###

# Propensity vector

prop.model <- c("a*b*Im1*(Sh1/(Sh1+Eh1+Ih1+Rh1+Vh1))","muh*Eh1","muh*Ih1",
                "muh*Rh1","muh*Vh1","w1*Sh1","alfah*Ih1","epsilonh*Eh1",
                "gamah*Ih1","a*c*Sm1*(Ih1/(Sh1+Eh1+Ih1+Rh1+Vh1))","mum*Em1",
                "mum*Im1","epsilonm*Em1","a*b*Im2*(Sh2/(Sh2+Eh2+Ih2+Rh2+Vh2))",
                "muh*Eh1","muh*Ih2","muh*Rh2","muh*Vh2","w2*Sh2","alfah*Ih2",
                "epsilonh*Eh2","gamah*Ih2","a*c*Sm2*(Ih2/(Sh2+Eh2+Ih2+Rh2+Vh2))",
                "mum*Em1","mum*Im2","epsilonm*Em2","death*Sh2","death*Eh2",
                "death*Ih2","death*Rh2","death*Vh2","a*b*Im3*(Sh3/(Sh3+Eh3+Ih3+Rh3+Vh3))",
                "muh*Eh1","muh*Ih3","muh*Rh3","muh*Vh3","w3*Sh3","alfah*Ih3",
                "epsilonh*Eh3","gamah*Ih3","a*c*Sm3*(Ih3/(Sh3+Eh3+Ih3+Rh3+Vh3))",
                "mum*Em3","mum*Im3","epsilonm*Em3")

# State change matrix
state.change.mat <- matrix(c(-1,1,1,1,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             1,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,-1,0,0,0,-1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,1,1,1,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,-1,1,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,1,1,1,-1,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,-1,1,-1,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,1,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,-1,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1),
                           nrow=24,byrow=T)



tf <- 1*365  # Final time
simName <- "model_stoch"


### Binomial tau-leap Method (Approximated Method) ----
# Stochastic Simulations = 50

set.seed(1234)
out <- ssa(state.model,prop.model,state.change.mat,
           par.model,tf,method=ssa.btl(),simName,verbose=F,
           consoleInterval=1)
outdf <- as.data.frame(out$data)


# PLOT Infected in the caravan
par(mfrow=c(1,2))

plot(mod2$t,mod2$Ih2,lwd=2, type = "l", ylim=c(0, 100),
     xlab="Time (days)", ylab="Number of cases",
     sub=paste("Vaccination =", pvac2))
lines(outdf$t, outdf$Ih2, type="l", 
     col="gray")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width
title(main="B)")
lines(mod2$t,mod2$Ih2,lwd=2)

for (i in 1:99){
  out <- ssa(state.model,prop.model,state.change.mat,
             par.model,tf,method=ssa.btl(),simName,verbose=F,
             consoleInterval=1)
  outdf <- as.data.frame(out$data)
  
  lines(outdf$t, outdf$Ih2, type="l", col="gray")
}

### ADD LINES FROM DETERMINISTIC MODEL.
lines(mod2$t,mod2$Ih2,lwd=2)
