# CARAVAN SEIR


# Loading packages
library(deSolve)
library(ggplot2)

# 1. PARAMETERS in daily rates
# From Esteva et al, 2019: https://onlinelibrary.wiley.com/doi/full/10.1002/cmm4.1059

a	<- 10/30   # Average biting rate	(Haemagogus and Aedes)
b <- 0.25   # Fraction of bites actually infective (Haem = 0.25 Aedes = 0.2)
c <- 0.4    # Vector susceptibility to the virus (Haemagogus = 0.4, Aedes = 0.25)
death <- 6200/(898000/9)/9/12/30 # Human death rate in caravans
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

# 2. INITIAL CONDITIONS #----

# community 1
Sh0_1 <- 1000; Sm0_1 <- Sh0_1*1.5; m1 <- Sm0_1/Sh0_1
R0_1 <- a*a*b*c*m1*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_1
# caravan
Sh0_2 <- 1000; Sm0_2 <- Sh0_2*1.5; m2 <- Sm0_2/Sh0_2
R0_2 <- a*a*b*c*m2*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_2
# community 3
Sh0_3 <- 1000; Sm0_3 <- Sh0_3*1.5; m3 <- Sm0_3/Sh0_3
R0_3 <- a*a*b*c*m3*epsilonm/(mum*(mum+epsilonm)*(muh+alfah+gamah)); R0_3


#### MODEL 1: RESIDENT ENDEMIC COMMUNITY ####----
# Efective Vaccinated proportion:
(1-(1/R0_1))*100

state.endemic <- c(Sh1=Sh0_1-(Sh0_1*pvac1), Eh1=0, Ih1=0, Rh1=0, 
                   Infh1=0, Vh1=Sh0_1*pvac1, Sm1=Sm0_1, Em1=0, Im1=1)

# Time of simulation
tsim <- 1*365 # For how long lasts a caravan? 1000h walking = 62 days/16h/day
Dt <- 1 #But not in the Haemagogus area


# Frequency-dependent transmission
model.endemic <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # Community
    dSh1 <- -a*b*Im1*(Sh1/(Sh1+Eh1+Ih1+Rh1+Vh1)) + muh*(Eh1+Ih1+Rh1+Vh1) -
      w1*Sh1 + alfah*Ih1
    dEh1 <- a*b*Im1*(Sh1/(Sh1+Eh1+Ih1+Rh1+Vh1)) - epsilonh*Eh1 - muh*Eh1 
    dIh1 <- epsilonh*Eh1 - gamah*Ih1 - muh*Ih1 - alfah*Ih1
    dRh1 <- gamah*Ih1 - muh*Rh1
    dInfh1 <- a*b*Im1*(Sh1/(Sh1+Eh1+Ih1+Rh1+Vh1))
    dVh1 <- w1*Sh1 - muh*Vh1 
    
    # Mosquitoes vectors in community 1
    dSm1 <- -a*c*Sm1*(Ih1/(Sh1+Eh1+Ih1+Rh1+Vh1)) + mum*(Em1 + Im1)
    dEm1 <- a*c*Sm1*(Ih1/(Sh1+Eh1+Ih1+Rh1+Vh1)) - (epsilonm + mum)*Em1
    dIm1 <- epsilonm*Em1 - mum*Im1
    
    
    # return the output of the model
    return(list(c(dSh1, dEh1, dIh1, dRh1, dInfh1, dVh1,
                  dSm1, dEm1, dIm1)))
    
  })
}

tempo <- seq(from=0,to=tsim,by=Dt)

mod1 <- ode(y = state.endemic, times = tempo, func = model.endemic, 
            parms = par.model, method = "lsoda")

mod1 <- as.data.frame(mod1)

names(mod1) <- c("t","Sh1","Eh1","Ih1","Rh1","Infh1","Vh1","Sm1","Em1","Im1")

#STOP
#### MODEL 2: CARAVAN ####----

Ih0_2 <- mod1$Ih1[mod1$t==20]; Ih0_2

state.caravan <- c(Sh2=Sh0_2, Eh2=0, Ih2=Ih0_2, Rh2=0, Infh2=0,Vh2=0,
               Sm2=Sm0_2, Em2=0, Im2=0)

# Time of simulation
tsim <- 1*356 # For how long lasts a caravan? 1000h walking = 62 days/16h/day
Dt <- 1 #But not in the Haemagogus area


# Frequency-dependent transmission
model.caravan <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    dSh2 <- -a*b*Im2*(Sh2/(Sh2+Eh2+Ih2+Rh2+Vh2)) + muh*(Eh2+Ih2+Rh2+Vh2) - 
      death*Sh2 - w2*Sh2
    dEh2 <- a*b*Im2*(Sh2/(Sh2+Eh2+Ih2+Rh2+Vh2)) - epsilonh*Eh2 - muh*Eh2 - death*Eh2
    dIh2 <- epsilonh*Eh2 - gamah*Ih2 - muh*Ih2 - alfah*Ih2 - death*Ih2
    dRh2 <- gamah*Ih2 - muh*Rh2 - death*Rh2
    dInfh2 <- a*b*Im2*(Sh2/(Sh2+Eh2+Ih2+Rh2+Vh2))
    dVh2 <- w2*Sh2 - muh*Vh2 
    
    # Mosquitoes vectors in community 2
    dSm2 <- -a*c*Sm2*(Ih2/(Sh2+Eh2+Ih2+Rh2+Vh2)) + mum*(Em2 + Im2)
    dEm2 <- a*c*Sm2*(Ih2/(Sh2+Eh2+Ih2+Rh2+Vh2)) - (epsilonm + mum)*Em2
    dIm2 <- epsilonm*Em2 - mum*Im2
    
    # return the output of the model
    return(list(c(dSh2, dEh2, dIh2, dRh2, dInfh2, dVh2,
                  dSm2, dEm2, dIm2)))
    
  })
}

tempo <- seq(from=0,to=tsim,by=Dt)

mod2 <- ode(y = state.caravan, times = tempo, func = model.caravan, 
           parms = par.model, method = "lsoda")

mod2 <- as.data.frame(mod2)

names(mod2) <- c("t","Sh2","Eh2","Ih2","Rh2","Infh2","Vh2",
                 "Sm2","Em2","Im2")

#STOP

#### MODEL 3: DISEASE FREE RESIDENT ----

Ih0_3 <- mod2$Ih2[mod2$t==20];Ih0_3
#Ih0_3 <- Ih0_2

# Considering currrent vaccination coverage

state.resident <- c(Sh3=Sh0_3, Eh3=0, Ih3=Ih0_3, Rh3=0, 
                    Infh3=0, Vh3=0, Sm3=Sm0_3, Em3=0, Im3=0)

# Time of simulation
tsim <- 1*365
Dt <- 1

# 3. DETERMINISTIC ROSS-MACDONALD MODEL #

# Frequency-dependent transmission
model.resident <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
# Human hosts in community 3
dSh3 <- -a*b*Im3*(Sh3/(Sh3+Eh3+Ih3+Rh3+Vh3)) + muh*(Eh3+Ih3+Rh3+Vh3) -
  w3*Sh3
dEh3 <- a*b*Im3*(Sh3/(Sh3+Eh3+Ih3+Rh3+Vh3)) - epsilonh*Eh3 - muh*Eh3
dIh3 <- epsilonh*Eh3 - gamah*Ih3 - muh*Ih3 - alfah*Ih3
dRh3 <- gamah*Ih3 - muh*Rh3
dInfh3 <- a*b*Im3*(Sh3/(Sh3+Eh3+Ih3+Rh3+Vh3))
dVh3 <- w3*Sh3-muh*Vh3

# Mosquitoes vectors in community 3
dSm3 <- -a*c*Sm3*(Ih3/(Sh3+Eh3+Ih3+Rh3+Vh3)) + mum*(Em3 + Im3)
dEm3 <- a*c*Sm3*(Ih3/(Sh3+Eh3+Ih3+Rh3+Vh3)) - (epsilonm + mum)*Em3
dIm3 <- epsilonm*Em3 - mum*Im3

# return the output of the model
return(list(c(dSh3, dEh3, dIh3, dRh3, dInfh3,dVh3,
              dSm3, dEm3, dIm3)))

  })
}

tempo <- seq(from=0,to=tsim,by=Dt)

mod3 <- ode(y = state.resident, times = tempo, func = model.resident, 
           parms = par.model, method = "lsoda")

mod3 <- as.data.frame(mod3)

names(mod3) <- c("t","Sh3","Eh3","Ih3","Rh3","Infh3","Vh3",
                 "Sm3","Em3","Im3")

#STOP
# 4. PLOTTING ----

par(mfrow=c(1,3))

plot(mod1$Ih1, main = "A) Endemic",
     type = "l", lwd=2, col = "grey",
     sub=paste("Vac. cover.=", pvac1,
               "Total cases =",round(max(mod1$Infh1))),
         xlab= "Time (days)",
     ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width

plot(mod2$t,mod2$Ih2, main = "B) Caravan",
     type = "l", lwd=2, col = "grey",
     sub=paste("Vac. cover.=", pvac2,
       "Total cases =",round(max(mod2$Infh2))),
     xlab= "Time (days)",
     ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width

matplot(mod3$Ih3, type = "l", lwd=2, col = "grey",
        main = "C) Disease-free",
        sub=paste("Vac. cover.=", pvac3,
                  "Total cases =",round(max(mod3$Infh3))),
        xlab= "Time (days)",
        ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width


#STOP

# Plotting vaccinating Caravans scenario ----

par(mfrow=c(1,3))

plot(mod2$Ih2, main = paste("C) Vaccination =", pvac2),
     type = "l", lwd=2, col = "grey",
     sub=paste("Total cases =",round(max(mod2$Infh2))),
     xlab= "Time (days)",
     ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width


# Ploting realistic data ----

par(mfrow=c(1,3))

plot(mod1$Ih1, main = "A) Endemic",
     type = "l", lwd=2, col = "grey",
     sub=paste(
               "Total cases =",round(max(mod1$Infh1))),
     xlab= "Time (days)",
     ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width

plot(mod2$t,mod2$Ih2, main = "B) Caravan",
     type = "l", lwd=2, col = "grey",
     sub=paste("(Lasting 31 days)",
               "Total cases =",round(max(mod2$Infh2))),
     xlab= "Time (days)",
     ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width
abline(v=20, size=2)

matplot(mod3$Ih3, type = "l", lwd=2, col = "grey",
        main = "C) Disease-free",
        sub=paste(
                  "Total cases =",round(max(mod3$Infh3))),
        xlab= "Time (days)",
        ylab= "Number of cases")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.2)      # Grid line width
