##########################################################################
#
# BLACK-BROWED ALBATROSS SURVIVAL ANALYSIS 2006-2017
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, June 2018
# requested by Andy Stanworth, Falklands Conservation

## simple model run on 11 June yielded implausible survival for successful breeders - need to include an at-sea state
## update 12 June 2018: simplified survival to not distinguish between failed and successful breeders
## update 12 June 2018: removed dead state because only 1 bird was ever observed dead

## update 31 July 2018: included temporally varying survival
## after chat with Andy Stanworth: take breeding success info with a grain of salt, as only failures reliably recorded

library(tidyverse)
library(jagsUI)
library(data.table)


##### LOAD FORMATTED RINGING DATA ###########
setwd("A:/RSPB/UKOT/Falklands/BBAL")
setwd("C:/STEFFEN/RSPB/UKOT/Falklands/BBAL")
raw<- fread("BBAL_ringing_data.csv", header=T)
raw<-raw %>% filter(!is.na(ring))
head(raw)


#### FORMAT FOR MULTISTATE MODEL ############

CH<-raw %>% gather(key="year",value="state",starts_with("20"))
tail(CH)

multistateCH<-CH %>% filter(grepl('fate',year)) %>%
  mutate(BREED=ifelse(is.na(as.numeric(state)), "fail","succ")) %>%
  mutate(YEAR=as.numeric(substr(year,1,4))) %>%
  dplyr::select(ring,YEAR,BREED)%>%
  arrange(ring,YEAR)
dim(multistateCH)
simpleCH<-CH %>% filter(!(grepl('fate',year))) %>%
  mutate(STATE=ifelse(state=="1","alive",ifelse(state=="dead","dead","notSeen"))) %>%
  mutate(YEAR=as.numeric(substr(year,1,4))) %>%
  dplyr::select(ring,YEAR,STATE)%>%
  arrange(ring,YEAR)
dim(simpleCH)
CHfull<-simpleCH
CHfull$BREED<-multistateCH$BREED
dim(CHfull)

## CHECK DUPLICATES ##
CHfull %>% mutate(count=1) %>% group_by(ring,YEAR) %>%
  summarise(dupl=sum(count)) %>% filter (dupl>1)


## CREATE INPUT FOR MODELLING ##
rawinput<- CHfull %>% mutate(state=ifelse(STATE=="alive",BREED,STATE)) %>%
  mutate(state=ifelse(is.na(state),"notSeen",state)) %>%
  select(ring, YEAR, state) 
head(rawinput)

rawinput %>% filter(is.na(state))

table(rawinput$state)    ## only 1 bird ever observed in 'dead' state!!










#########################################################################
# Specify very basic CJS model to replicate Andy's analysis
#########################################################################

sink("BBAL_CJS.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability for breedingadults
    # p: recapture probability when breeding
    # -------------------------------------------------

    # Priors and constraints
    #for (i in 1:nind){
      for (t in 1:(n.occasions-1)){
        phi[t] ~ dunif(0.7, 1) 
        p[t] ~ dunif(0, 1) 
      } #t
    #} #i
    
    # for (t in 1:(n.occasions-1)){
    #   ann.phi[t] ~ dunif(0.7, 1)         # Prior for mean survival
    #   ann.p[t] ~ dunif(0, 1)           # Prior for mean recapture
    # } #t
    
    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
          # State process
          z[i,t] ~ dbern(mu1[i,t])
          mu1[i,t] <- phi[t-1] * z[i,t-1]
          # Observation process
          y[i,t] ~ dbern(mu2[i,t])
          mu2[i,t] <- p[t-1] * z[i,t]
        } #t
      } #i
    }
    ",fill = TRUE)
sink()



#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################
head(CHfull)
unique(CHfull$STATE)

#### REPLACE STATE NAMES WITH NUMBERED STATES
input<- CHfull %>% mutate(state=ifelse(STATE=="alive",1,0)) %>%
  mutate(state=ifelse(is.na(state),0,state)) %>%
  dplyr::select(ring, YEAR, state) %>% 
  group_by(ring) %>%
  spread(key=YEAR, value=state)
head(input)
rCH<-as.matrix(input[,-1])
head(rCH)
tail(rCH)

# ELIMINATE CAPTURE HISTORIES OF BIRDS ONLY RINGED IN THE LAST OCCASION
del <- apply(rCH[,1:11], 1, sum)
dim(rCH)
rCH<-rCH[!(del==0),]
dim(rCH)

# ELIMINATE TRANSIENTS ONLY OBSERVED IN A SINGLE YEAR
del <- apply(rCH, 1, sum)
dim(rCH)
rCH<-rCH[!(del==1),]
dim(rCH)


# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(rCH, 1, get.first)

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values 
inits <- function(){list(phi = runif(11, 0.7, 1),
                         pp = runif(11, 0, 1))}
 

# Parameters monitored
parameters <- c("phi", "p")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 5000
nc <- 4

# Call JAGS from R
#BBALsurv <- jags(jags.data, inits, parameters, "S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Falklands\\BBAL\\BBAL_surv_simple.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
BBALsurv <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Falklands\\BBAL\\BBAL_CJS.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)




#########################################################################
# PRODUCE OUTPUT TABLE
#########################################################################

out<-as.data.frame(BBALsurv$summary)
out$parameter<-row.names(BBALsurv$summary)
write.table(out,"BBAL_Falklands_Survival_estimates.csv", sep=",", row.names=F)



#########################################################################
# PRODUCE OUTPUT GRAPH
#########################################################################


pdf("BBAL_survival_Falklands.pdf", width=11, height=7)
out[1:11,] %>% select(c(1,5,2,3,7)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl')) %>%
  mutate(Year=colnames(rCH)[1:11]) %>%
  
  ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  ylab("Annual adult survival probability") +
  scale_y_continuous(breaks=seq(0.5,1,0.1), limits=c(0.5,1))+
  #scale_x_continuous(breaks=seq(2006,2017,1))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()








#########################################################################
# ATTEMPT TO FIT MULTI-STATE MODEL WITH BREEDING SUCCESS
#########################################################################
### this was discarded on 31 July 2018
### estimates are nonsensical (survival <0.85)
### likely issue is that breeding success not conclusively recorded
### non-breeding state missing?


sink("BBAL_surv_time.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability for adults
    # breedsucc: probability of breeding successfully (regardless of previous success)
    # pS: recapture probability when breeding successfully
    # pF: recapture probability when not breeding or failed
    # -------------------------------------------------
    # States (S):
    # 1 successful breeder
    # 2 unsuccessful breeder
    # 3 dead
    # Observations (O):
    # 1 seen as successful 
    # 2 seen as unsuccessful or nonbreeder
    # 3 not seen
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival and recapture: uniform
    
    for (t in 1:(n.occasions-1)){
    phi[t] ~ dunif(0.7, 1)           ## excluded unrealistic survival probabilities
    breedsucc[t] ~ dunif(0, 1)
    }
    pS ~ dunif(0.8, 1)             ## very high probability of detecting a breeder
    pF ~ dunif(0, 0.9)               ## lower probability of detecting failed and nonbreeders
    
    # Define state-transition and observation matrices 	
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- phi[t] * breedsucc[t]
    ps[1,i,t,2] <- phi[t] * (1-breedsucc[t])
    ps[1,i,t,3] <- (1-phi[t])
    
    ps[2,i,t,1] <- phi[t] * breedsucc[t]
    ps[2,i,t,2] <- phi[t] * (1-breedsucc[t])
    ps[2,i,t,3] <- (1-phi[t])
    
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- pS
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- (1-pS)
    
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pF
    po[2,i,t,3] <- (1-pF)
    
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture   
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()



#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

#### REPLACE STATE NAMES WITH NUMBERED STATES
input<- CHfull %>% mutate(state=ifelse(STATE=="alive",BREED,STATE)) %>%
  mutate(state=ifelse(is.na(state),"notSeen",state)) %>%
  dplyr::select(ring, YEAR, state) %>% 
  #mutate(state=ifelse(state=="succ",1,ifelse(state=='fail',2,ifelse(state=="notSeen",4,3)))) %>%      ### removed state 3 'dead'
  mutate(state=ifelse(state=="succ",1,ifelse(state=='fail',2,3))) %>%
  group_by(ring) %>%
  spread(key=YEAR, value=state)
head(input)
rCH<-as.matrix(input[,-1])
head(rCH)
tail(rCH)
head(z)

# ELIMINATE CAPTURE HISTORIES OF BIRDS NEVER OBSERVED
eliminate <- function(x) min(x)<3
del <- apply(rCH, 1, eliminate)
dim(rCH)
rCH<-rCH[del==TRUE,]
dim(rCH)

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=3))
f <- apply(rCH, 1, get.first)


# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for not seen
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)     ## only states 1 and 2 are an option as 3 = dead
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values 
inits <- function(){list(phi = runif(11, 0.7, 1),
                         breedsucc = runif(11, 0, 1),
                         pS = runif(1, 0.8, 1),
                         pF = runif(1, 0, 0.9),
                         z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("phi", "breedsucc", "pF", "pS")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 5000
nc <- 4

# Call JAGS from R
#BBALsurv <- jags(jags.data, inits, parameters, "S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Falklands\\BBAL\\BBAL_surv_simple.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
BBALsurv <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Falklands\\BBAL\\BBAL_surv_time.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)



#########################################################################
# PRODUCE OUTPUT TABLE
#########################################################################

out<-as.data.frame(BBALsurv$summary)
out$parameter<-row.names(BBALsurv$summary)
write.table(out,"BBAL_Falklands_Survival_estimates.csv", sep=",", row.names=F)










