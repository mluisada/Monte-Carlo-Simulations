###############################
#   Monte-Carlo Simulations   #
#    Asian Options Pricing    # 
#           May 2014          #
#                             #
#       ENSAE ParisTech       #
#        Group project:       #
#      COLUMELLI Capucine     #
#     LUISADA Marie-Laura     #
#        TANGUY Charles       #
###############################


##############
#   Set-up   #
##############

### Library to generate pseudo / quasi random numbers
#install.packages("randtoolbox")
library(randtoolbox)



#####################
#  Basic Functions  #
#####################


#### CREATE A FUNCTION THAT GENERATES A MATRIX SO THAT
#### EACH ROW IS A BROWNIAN MOTION 
#### WITH k timestamps between 0 and T (= expiration date)

# function arguments : k = nb of dates, T = expiration date, n_sim = nb of simulations
Brownian_Motion <- function(k,T,n_sim){
  # initialize 
  BM <- array(0,dim=c(n_sim,k))
  # for each simulation
  for (simul in 1:n_sim){
    # define brownian motion
    BM[simul,2:k] <- cumsum(rnorm(k-1,0,1)*sqrt(T/k))
  }
  return(BM)
}
Brownian_Motion(5,1,10)



#### CREATE A FUNCTION WHICH TAKES SOME DATA AS AN INPUT
#### AND GENERATE A MATRIX OF UNDERLYING ASSETS
#### BASED ON BLACK AND SCHOLES MODEL

# N.B. Comparing the underlying to the strike enables to characterize the option

# function arguments : S_0 = initial underlying, mu = mean, sigma = standard deviation, 
#                      T = expiration date, data
Black_scholes <- function(S_0, mu, sigma, T, data){
  n_sim <- dim(data)[1]
  k <- dim(data)[2]
  # Underlying assets
  S <- array(0,dim=dim(data))
  for (simul in 1:n_sim){
    for (date in 1:k){
      # Black and Scholes
      S[simul,date] <- S_0*exp((mu-(sigma^2/2))*(date*T/k) + sigma*data[simul,date])
    }
  }
  return(S)
}



#### CREATE A FUNCTION WHICH TAKES SOME DATA AS AN INPUT
#### AND GENERATE A MATRIX OF DISCRETIZED UNDERLYING ASSETS
#### BASED ON EULER METHOD

# function arguments : S_0 = initial underlying, mu = mean, sigma = standard deviation, 
#                      T = expiration date, data
Black_scholes_disc <- function(S_0, mu, sigma, T, data){
  n_sim <- dim(data)[1]
  k <- dim(data)[2]
  S <- array(0,dim=dim(data))
  for (simul in 1:n_sim){
    S[simul,1] <- S_0
    for (date in 2:k){
      # Euler method
      S[simul,date] <- S[simul,date - 1]*(1+ mu*T/k + sigma*(data[simul,date]-data[simul,date - 1]))
    }
  }
  return(S)
}


#### CREATE A FUNCTION WHICH CHARACTERIZES THE OPTION
#### BY COMPARYING THE STRIKE TO THE UNDERLYING

# function arguments : r = interest rate, strike's value, underlying's value, T = expiration date
CalculC <- function(r, strike, underlying, T){
  # nb of simulation is the number of underlying assets
  n_sim <- dim(underlying)[1]
  C <- array(0,dim=c(n_sim,1))
  for (simul in 1:n_sim){
    C[simul] <- max(0,exp(-r*T)*(mean(underlying[simul,]) - strike))
  }
  return(C)
}

#### CREATE A FUNCTION WHICH PRINTS A HISTOGRAM OF SIMULATIONS,
#### PLOTS THE DENSITY FUNCTION 
#### AND DISPLAYS THE MEAN AND STANDARD DEVIATION

Estim <- function(vecteur){
  hist(vecteur, breaks = 30, freq=FALSE, col="lightblue", xlab = "estimator")
  lines(density(vecteur),col="red")
  print(mean(vecteur))
  print(length(vecteur)/(length(vecteur)-1)*mean((vecteur-mean(vecteur))^2))
}


#### APPLICATION 
set.seed(1234)
testwt <- Brownian_Motion(50,1,1000)
testS <- Black_scholes(10,0,0.3,1,testwt)
testc <- CalculC(0.005,5,testS,1)
Estim(testc)

testwt <- Brownian_Motion(50,1,1000)
testS_disc <- Black_scholes_disc(10,0,0.3,1,testwt)
testc_disc <- CalculC(0.005,5,testS_disc,1)
Estim(testc_disc)



##########################
#  Antithetic variables  #
##########################

#### BUILD A FUNCTION WHICH RETURNS
#### ANTITHETIC VARIABLES

Var_anti <- function(k,T,n_sim){
  VA <- array(0,dim=c(n_sim*2,k))
  VA[1:n_sim,] <- Brownian_Motion(k,T,n_sim)
  VA[(n_sim+1):(n_sim*2),] <- (-1)*VA[1:n_sim,]
  return(VA)
}

### APPLICATION
set.seed(1234)
testVA <- Var_anti(50,1,1000)
testSVA <- Black_scholes(10,0,0.3,1,testVA)
testcVA <- CalculC(0.005,5,testSVA,1)
Estim(testcVA)

set.seed(1234)
testVA <- Var_anti(50,1,1000)
testSVA_disc <- Black_scholes_disc(10,0,0.3,1,testVA)
testcVA_disc <- CalculC(0.005,5,testSVA_disc,1)
Estim(testcVA_disc)


#######################
#  Control variables  #
#######################

# Continuous
control_var <- function(mu,sigma,intrate,strike,k,T,n_sim,S_0){
  testwt <- Brownian_Motion(k,T,n_sim)
  testS <- Black_scholes(S_0,mu,sigma,T,testwt)
  testc <- CalculC(intrate,strike,testS,T)
  testexpwt <- array(0,dim=c(n_sim,1))
  for (i in 1:n_sim){
    testexpwt[i] <- mean(exp(testwt[i,]))
  }
  coeff <- cov(testc,testexpwt)/var(testexpwt)
  print(coeff)
  for (simul in 1:n_sim){                                  
    testc[simul] <- max(testc[simul] - coeff*(testexpwt[simul]-mean(testexpwt)),0)
  }
  return(testc)
}

# Discrete
control_var_disc <- function(mu,sigma,intrate,strike,k,T,n_sim,S_0){
  testwt <- Brownian_Motion(k,T,n_sim)
  testS <- Black_scholes_disc(S_0,mu,sigma,T,testwt)
  testc <- CalculC(intrate,strike,testS,T)
  testexpwt <- array(0,dim=c(n_sim,1))
  for (i in 1:n_sim){
    testexpwt[i] <- mean(exp(testwt[i,]))
  }
  coeff <- cov(testc,testexpwt)/var(testexpwt)
  print(coeff)
  for (simul in 1:n_sim){                                  
    testc[simul] <- max(testc[simul] - coeff*(testexpwt[simul]-mean(testexpwt)),0)
  }
  return(testc)
}


set.seed(1234)
testCV <- control_var(0,0.3,0.005,5,50,1,1000,10)
Estim(testCV)

testCV_disc <- control_var_disc(0,0.3,0.005,5,50,1,1000,10)
Estim(testCV_disc)


#######################
#  Quasi Montecarlo   #
#######################

QMC <- function(k,T,n_sim){
  BMQMC <- array(0,dim=c(n_sim,k))
  for (i in 1:n_sim){
    halton <- sqrt(T/k)*halton(k, dim = 1, init = TRUE, normal = TRUE, usetime = FALSE)    
    for (date in 1:k){
      BMQMC[i,] <- cumsum(halton)
    }
  }
  return(BMQMC)
}

set.seed(1234)
testQMC <- QMC(50,1,1000)
testSQMC <- Black_scholes(10,0,0.3,1,testQMC )
testcQMC <- CalculC(0.005,5,testSQMC ,1)
#Estim(testcQMC)


set.seed(1234)
testQMC <- QMC(50,1,1000)
testSQMC_disc <- Black_scholes_disc(10,0,0.3,1,testQMC )
testcQMC_disc <- CalculC(0.005,5,testSQMC_disc ,1)
#Estim(testcQMC_disc )

expectationMC<-c()
expectationAV<-c()
expectationCV<-c()

j=1
for (i in seq(1,10000,by=500)){
  testMC<- Brownian_Motion(50,1,i)
  testMC <- Black_scholes(10,0,0.3,1,testMC)
  testMC <- CalculC(0.005,5,testMC,1)
  testAV<- Var_anti(50,1,i)
  testAV<- Black_scholes(10,0,0.3,1,testAV)
  testAV<- CalculC(0.005,5,testAV,1)
  testCV<- control_var(0,0.3,0.005,5,50,1,i,10)
  expectationMC[j]<-mean(testMC)
  expectationAV[j]<-mean(testAV)
  expectationCV[j]<-mean(testCV)
  j<-j+1
}

plot(expectationMC,col="red",type="l")
lines(expectationAV,col="blue")
lines(expectationCV,col="green")

