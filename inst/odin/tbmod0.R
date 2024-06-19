################## just TB as a warm up

################## parameters

## dimensions

## number of PPD states
NP <- 5
## 1 = remand
## 2 = short stay
## 3 = long stay
## 4 = open
## 5 = previously detained

## number of TPT states
NT <- 3
## 1 = never TPT
## 2 = on TPT
## 3 = previous TPT


## PPD transitions
inflow <- user()          # inflow rate (less recidivists)
remand_short <- user()    # remand -> short:   1->2
remand_long <- user()     # remand -> long:    1->3
remand_release <- user()  # remand -> release: 1->5
long_short <- user()      # long -> short:     3->2
short_release <- user()   # short -> release:  2->5
short_open <- user()      # short -> open:     2->4
long_release <- user()    # long -> release:   3->5
open_release <- user()    # open -> release:   4->5
previous_remand <- user() # previous -> remand 5->1


## TB splits on arrival
parm_frac_U <- user() #uninfected
parm_frac_E <- user() #early latent
parm_frac_L <- user() #late latent
parm_frac_SD <- user() #subclinical disease
parm_frac_CD <- user() #clinical disease
parm_frac_ATT <- user() #AntiTB treatment
parm_frac_epTB <- user() #early post-TB
parm_frac_lpTB <- user() #late post-TB

## print("parm_frac_lpTB : {parm_frac_lpTB}")

## initial TB states
parm_ifrac_U <- user() #uninfected
parm_ifrac_E <- user() #early latent
parm_ifrac_L <- user() #late latent
parm_ifrac_SD <- user() #subclinical disease
parm_ifrac_CD <- user() #clinical disease
parm_ifrac_ATT <- user() #AntiTB treatment
parm_ifrac_epTB <- user() #early post-TB
parm_ifrac_lpTB <- user() #late post-TB
parm_ifrac_prevTPT[] <- user() #TPT state


## initial PPD split
parm_init_PPD[] <- user()

## multiplier for transmission once out
m <- user(1.0)


################## initial state
initial(U[1:NP,1:NT]) <- parm_ifrac_U * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##uninfected
initial(E[1:NP,1:NT]) <- parm_ifrac_E * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##early latent
initial(L[1:NP,1:NT]) <- parm_ifrac_L * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##late latent
initial(SD[1:NP,1:NT]) <- parm_ifrac_SD * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##subclinical disease
initial(CD[1:NP,1:NT]) <- parm_ifrac_CD * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##clinical disease
initial(ATT[1:NP,1:NT]) <- parm_ifrac_ATT * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##AntiTB treatment
initial(epTB[1:NP,1:NT]) <- parm_ifrac_epTB * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##early post-TB
initial(lpTB[1:NP,1:NT]) <- parm_ifrac_lpTB * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##late post-TB

## ## test
## initial(Ntot[,]) <- parm_init_PPD[i] * parm_ifrac_prevTPT[j]

################## dynamics

## === TB states
## uninfected
deriv(U[,]) <-  - infections[i,j] - mort*U[i,j] +
  inflow * frac_U[i,j] + moves_U[i,j]
## early latent
deriv(E[,]) <- infections[i,j] + reinfections[i,j] -
  stabilizations[i,j] - fastprogs[i,j] - mort*E[i,j] +
  inflow * frac_E[i, j] + moves_E[i, j] + tpt_E[i, j]
## late latent
deriv(L[,]) <- - Lreinfections[i,j] + stabilizations[i,j] + selfcures[i,j] -
  slowprogs[i,j] - mort*L[i,j] +
  inflow * frac_L[i, j] + moves_L[i, j] + tpt_L[i, j]
## subclinical disease
deriv(SD[,]) <- fastprogs[i,j] + slowprogs[i,j] - mort*SD[i,j] - worsens[i,j] +
  inflow * frac_SD[i, j] + moves_SD[i, j]
## clinical disease
deriv(CD[,]) <- worsens[i,j]- mort*CD[i,j] + relapses[i,j] +
  inflow * frac_CD[i, j] + moves_CD[i, j] - tbstops[i, j]
## AntiTB treatment
deriv(ATT[,]) <- detects[i,j] - ATT[i,j]/att_time - mort*ATT[i,j] +
  inflow * frac_ATT[i, j] + moves_ATT[i, j]
## early post-TB
deriv(epTB[,]) <- (1-txf)*ATT[i,j]/att_time -
  relapses[i,j] - epTB[i,j]/late_post_time - mHR*mort*epTB[i,j] +
  inflow * frac_epTB[i, j] + moves_epTB[i, j] + tpt_epTB[i, j]
## late post-TB
deriv(lpTB[,]) <- epTB[i,j]/late_post_time - mHR*mort*lpTB[i,j] - Preinfections[i,j] +
  inflow * frac_lpTB[i, j] + moves_lpTB[i, j] + tpt_lpTB[i, j]

## === economic states & counters
## economic parameters
disc_rate <- user(0.035)
LifeExp <- user(40)

## unit costs
uc_screening <- user(1) # LTBI screening at entry
uc_tpt <- user(1)       # TPT following screening
uc_attscreen <- user(1) # ATT for those found via screening at entry
uc_attppd <- user(1)    # ATT for those found passively within the system
uc_attout <- user(1)    # ATT following release

## HRQoL
hrqol <- user(0.3)        # HRQoL decrement while CD
hrqolptb <- user(0.05)     # HRQoL decrement while post TB


## inflow prevalences/fractions:
inflow_LTBI <- (parm_frac_E + parm_frac_L + parm_frac_epTB + parm_frac_lpTB) #LTBI = TPT eligible
inflow_TBD <- (parm_frac_SD + parm_frac_CD)                                  #TB disease = ATT eligible
inflow_rest <- (1-inflow_TBD-inflow_LTBI-parm_frac_ATT)                      #rest, excl ATT

## economic states
initial(CC0) <- 0 # undiscounted cumulative costs
initial(CC) <- 0 # discounted cumulative costs
initial(cATTtp) <- 0 # ATT true positive counter
initial(cATTfp) <- 0 # ATT false positive counter
initial(cTPT) <- 0 # cumulative TPT counter
initial(dLYL) <- 0 # discounted LYL
initial(deaths) <- 0 # TB deaths
initial(qoldec) <- 0 # QoL decrement due to TB

## TODO this will need stratifying by TB & FPs included?
## dynamics
deriv(CC0) <- (((uc_screening + uc_tpt * inflow_LTBI) * inflow_toTPT_L +
                uc_attscreen * inflow_TBD * inflow_toATT_TB) * inflow +
               uc_attppd * sum(detects[1:4,1:NT]) + m * uc_attout * sum(detects[5,1:NT])) *
  if (t > int_time) 1 else 0
deriv(CC) <- (((uc_screening + uc_tpt * inflow_LTBI) * inflow_toTPT_L +
                uc_attscreen * inflow_TBD * inflow_toATT_TB) * inflow +
              uc_attppd * sum(detects[1:4,1:NT]) + m * uc_attout * sum(detects[5,1:NT])) *
  if (t > int_time) exp(-(t - int_time) * disc_rate) else 0

## deriv(CC0) <- ((uc_screening*scr_frac+ uc_tpt*tpt_frac+ uc_attscreen * att_frac)* inflow +
##                uc_attppd * sum(detects[1:4,1:NT]) + m * uc_attout * sum(detects[5,1:NT])) *
##   if (t > int_time) 1 else 0

deriv(cATTtp) <- 0 #ATT true positive counter
deriv(cATTfp) <- 0 #ATT false positive counter
deriv(cTPT) <- 0   #cumulative TPT counter
deriv(dLYL) <- (sum(tbmort)+(m-1)*sum(tbmort[5,1:NT])) * #includes extra outside
  (1 - exp(-disc_rate * LifeExp)) / (disc_rate + 1e-15) *
  if (t > int_time) exp(-(t - int_time) * disc_rate) else 0
deriv(deaths) <- (sum(tbmort)+(m-1)*sum(tbmort[5,1:NT])) * #includes extra outside
  if (t > int_time) 1 else 0
deriv(qoldec) <- (sum(qolrate) + (m-1)*sum(qolrate[5,1:NT])) * #includes extra outside
  if (t > int_time) exp(-(t - int_time) * disc_rate) else 0

## ## test
## deriv(Ntot[, ]) <- inflow * inflow_top[i] * inflow_TPTv[j] + moves_Ntot[i, j]

## === extras
bet <- foi / (parm_ifrac_SD + parm_ifrac_CD) ## compute beta from prevalence and foi
ppop[, ] <- U[i,j] + E[i,j] + L[i,j] + SD[i,j] + CD[i,j] + ATT[i,j] + epTB[i,j] + lpTB[i,j] #population
fi <- if (staticfoi > 0) foi else bet * (sum(SD) + sum(CD)) / (sum(ppop[1:4, 1:NT]) + 1e-2) # FOI
rfv[1:4] <- 1
rfv[5] <- 0
staticfoi <- user(1) #is the model static


## === TB processes
infections[,] <- fi * rfv[i] * U[i,j]
Lreinfections[, ] <- ptn * fi * rfv[i] * L[i, j]
Ereinfections[, ] <- ptn * fi * rfv[i] * epTB[i, j]
Preinfections[, ] <- ptn * fi * rfv[i] * lpTB[i, j]
reinfections[,] <- Ereinfections[i,j] + Lreinfections[i,j] + Preinfections[i,j]
stabilizations[,] <- stb * E[i,j]
fastprogs[,] <- prg * E[i,j] * hrv[j]
slowprogs[,] <- eps * L[i,j] * hrv[j]
worsens[,] <- SD[i,j] / wsn
relapses[,] <- rel * epTB[i,j]
selfcures[,] <- (1-CFR) * CD[i,j] / drn
detects[,] <- (CDR/(1-CDR)) * CD[i,j] / drn
tbstops[,] <- CD[i,j]/(drn*(1-CDR)) #see notes around parameters
tbmort[,] <- CFR * CD[i,j]/drn +
  txf * ATT[i,j]/att_time +
  (mHR-1)*mort*epTB[i,j] +
  (mHR-1)*mort*lpTB[i,j]

## driver of hrqol
qolrate[,] <- hrqol * CD[i,j] + hrqolptb * epTB[i,j] + hrqolptb * lpTB[i,j]

## TB parameters:
## a = 1 / d = rate out with no detection
## ATT : no-ATT = CDR : (1-CDR)
## total rate w/detection = a + b; CDR=b/(a+b) -> b = CDR/(1-CDR)* a -> a+b = 1/(1-CDR) / d
## 1 / ((1-CDR)*d) = rate out with detection
## notes ~ b * D = (CDR/(1-CDR)) * D / d
## -> non-notes ~  D / d
## deaths = cfr * D / d; self-cure = (1-cfr) * D / d
rel <- user() #relapse rate
eps <- user() #slow progn rate
prg <- user() #fast progn rate
stb <- user() #early latent stabilization rate
CDR <- user() #detection
CFR <- user() #CFR untreated TB
drn <- user() #TB duration untreated (clinical)
wsn <- user() #duration subclinical symptoms
mHR <- user() #post-TB mortality HR
txf <- user() #TB mortality on ATT
ptn <- user() #HR protection of reinfection LTBI
foi <- user(0.01)

## other parameters
att_time <- user(0.5)     #duration of ATT
late_post_time <- user(2) #duration defining early post-TB
mort <- user(0.02)        #mortality rate
tptHR <- user(0.3)        #HR for TPT
hrv[1] <- 1
hrv[2] <- tptHR #vector to apply to middle
hrv[3] <- 1

## === intervention
## screening coverages:
##     attributes Uninfected, LTBI, TB, Previous
## outcomes:
##     ATT, TPT, no change
## timing:
int_time <- user() # time intervention is turned on & outputs calculated from
## probabilities within inflow:
## baseline/SOC
inflow_toATT_TB0 <- user(0) # NOTE fp ATT doesn't affect state
inflow_toTPT_L0 <- user(0) # assume perfect spec
## intervention
inflow_toATT_TB1 <- user(0) # NOTE fp ATT doesn't affect state
inflow_toTPT_L1 <- user(0) # assume perfect spec
## switch
inflow_toATT_TB <- if (t > int_time) inflow_toATT_TB1 else inflow_toATT_TB0
inflow_toTPT_L <- if (t > int_time) inflow_toTPT_L1 else inflow_toTPT_L0



## vector form for compactness
inflow_TPTv[1] <- (1-inflow_toTPT_L)
inflow_TPTv[2] <- (inflow_toTPT_L)
inflow_TPTv[3] <- 0
## inflow to remand
inflow_top[1] <- 1
inflow_top[2:NP] <- 0
## the non-TPT layers
TPT_top[1] <- 1
TPT_top[2:3] <- 0


## === TB split on arrival, including TPT
frac_U[,] <- parm_frac_U * inflow_top[i] * TPT_top[j] #uninfected
frac_E[,] <- parm_frac_E * inflow_top[i] * inflow_TPTv[j] #early latent
frac_L[,] <- parm_frac_L * inflow_top[i] * inflow_TPTv[j] #late latent
frac_SD[, ] <- parm_frac_SD * inflow_top[i] * (1 - inflow_toATT_TB) * TPT_top[j] # subclinical disease
frac_CD[, ] <- parm_frac_CD * inflow_top[i] * (1 - inflow_toATT_TB) * TPT_top[j] # clinical disease
frac_ATT[, ] <- (parm_frac_ATT + (parm_frac_SD + parm_frac_CD) * inflow_toATT_TB) *
  inflow_top[i] * TPT_top[j] # ATT
frac_epTB[,] <- parm_frac_epTB * inflow_top[i] * inflow_TPTv[j] #early post-TB
frac_lpTB[,] <- parm_frac_lpTB * inflow_top[i] * inflow_TPTv[j] #late post-TB


## === PPD processes: transitions in 1st index

## 1 = remand
moves_U[1, ] <- -(remand_short + remand_long + remand_release) * U[1, j] +
  previous_remand * U[5, j]
moves_E[1, ] <- -(remand_short + remand_long + remand_release) * E[1, j] +
  previous_remand * E[5, j]
moves_L[1, ] <- -(remand_short + remand_long + remand_release) * L[1, j] +
  previous_remand * L[5, j]
moves_SD[1, ] <- -(remand_short + remand_long + remand_release) * SD[1, j] +
  previous_remand * SD[5, j]
moves_CD[1, ] <- -(remand_short + remand_long + remand_release) * CD[1, j] +
  previous_remand * CD[5, j]
moves_ATT[1, ] <- -(remand_short + remand_long + remand_release) * ATT[1, j] +
  previous_remand * ATT[5, j]
moves_epTB[1, ] <- -(remand_short + remand_long + remand_release) * epTB[1, j] +
  previous_remand * epTB[5, j]
moves_lpTB[1, ] <- -(remand_short + remand_long + remand_release) * lpTB[1, j] +
  previous_remand * lpTB[5, j]

## 2 = short stay
moves_U[2, ] <- -(short_release + short_open) * U[2, j] +
  remand_short * U[1, j] + long_short * U[3, j]
moves_E[2, ] <- -(short_release + short_open) * E[2, j] +
  remand_short * E[1, j] + long_short * E[3, j]
moves_L[2, ] <- -(short_release + short_open) * L[2, j] +
  remand_short * L[1, j] + long_short * L[3, j]
moves_SD[2, ] <- -(short_release + short_open) * SD[2, j] +
  remand_short * SD[1, j] + long_short * SD[3, j]
moves_CD[2, ] <- -(short_release + short_open) * CD[2, j] +
  remand_short * CD[1, j] + long_short * CD[3, j]
moves_ATT[2, ] <- -(short_release + short_open) * ATT[2, j] +
  remand_short * ATT[1, j] + long_short * ATT[3, j]
moves_epTB[2, ] <- -(short_release + short_open) * epTB[2, j] +
  remand_short * epTB[1, j] + long_short * epTB[3, j]
moves_lpTB[2, ] <- -(short_release + short_open) * lpTB[2, j] +
  remand_short * lpTB[1, j] + long_short * lpTB[3, j]

## 3 = long stay
moves_U[3, ] <- -(long_short + long_release) * U[3, j] +
  remand_long * U[1, j]
moves_E[3, ] <- -(long_short + long_release) * E[3, j] +
  remand_long * E[1, j]
moves_L[3, ] <- -(long_short + long_release) * L[3, j] +
  remand_long * L[1, j]
moves_SD[3, ] <- -(long_short + long_release) * SD[3, j] +
  remand_long * SD[1, j]
moves_CD[3, ] <- -(long_short + long_release) * CD[3, j] +
  remand_long * CD[1, j]
moves_ATT[3, ] <- -(long_short + long_release) * ATT[3, j] +
  remand_long * ATT[1, j]
moves_epTB[3, ] <- -(long_short + long_release) * epTB[3, j] +
  remand_long * epTB[1, j]
moves_lpTB[3, ] <- -(long_short + long_release) * lpTB[3, j] +
  remand_long * lpTB[1, j]

## 4 = open
moves_U[4, ] <- -(open_release) * U[4, j] +
  short_open * U[2, j]
moves_E[4, ] <- -(open_release) * E[4, j] +
  short_open * E[2, j]
moves_L[4, ] <- -(open_release) * L[4, j] +
  short_open * L[2, j]
moves_SD[4, ] <- -(open_release) * SD[4, j] +
  short_open * SD[2, j]
moves_CD[4, ] <- -(open_release) * CD[4, j] +
  short_open * CD[2, j]
moves_ATT[4, ] <- -(open_release) * ATT[4, j] +
  short_open * ATT[2, j]
moves_epTB[4, ] <- -(open_release) * epTB[4, j] +
  short_open * epTB[2, j]
moves_lpTB[4, ] <- -(open_release) * lpTB[4, j] +
  short_open * lpTB[2, j]

## 5 = previously detained
moves_U[5, ] <- -(previous_remand) * U[5, j] +
  remand_release * U[1, j] + short_release * U[2, j] +
  long_release * U[3, j] + open_release * U[4, j]
moves_E[5, ] <- -(previous_remand) * E[5, j] +
  remand_release * E[1, j] + short_release * E[2, j] +
  long_release * E[3, j] + open_release * E[4, j]
moves_L[5, ] <- -(previous_remand) * L[5, j] +
  remand_release * L[1, j] + short_release * L[2, j] +
  long_release * L[3, j] + open_release * L[4, j]
moves_SD[5, ] <- -(previous_remand) * SD[5, j] +
  remand_release * SD[1, j] + short_release * SD[2, j] +
  long_release * SD[3, j] + open_release * SD[4, j]
moves_CD[5, ] <- -(previous_remand) * CD[5, j] +
  remand_release * CD[1, j] + short_release * CD[2, j] +
  long_release * CD[3, j] + open_release * CD[4, j]
moves_ATT[5, ] <- -(previous_remand) * ATT[5, j] +
  remand_release * ATT[1, j] + short_release * ATT[2, j] +
  long_release * ATT[3, j] + open_release * ATT[4, j]
moves_epTB[5, ] <- -(previous_remand) * epTB[5, j] +
  remand_release * epTB[1, j] + short_release * epTB[2, j] +
  long_release * epTB[3, j] + open_release * epTB[4, j]
moves_lpTB[5, ] <- -(previous_remand) * lpTB[5, j] +
  remand_release * lpTB[1, j] + short_release * lpTB[2, j] +
  long_release * lpTB[3, j] + open_release * lpTB[4, j]

## ## TEST
## ## 1 = remand
## moves_Ntot[1, ] <- -(remand_short + remand_long + remand_release) * Ntot[1, j] +
##   previous_remand * Ntot[5, j]
## ## 2 = short stay
## moves_Ntot[2, ] <- -(short_release + short_open) * Ntot[2, j] +
##   remand_short * Ntot[1, j] + long_short * Ntot[3, j]
## ## 3 = long stay
## moves_Ntot[3, ] <- -(long_short + long_release) * Ntot[3, j] +
##   remand_long * Ntot[1, j]
## ## 4 = open
## moves_Ntot[4, ] <- -(open_release) * Ntot[4, j] +
##   short_open * Ntot[2, j]
## ## 5 = previously detained
## moves_Ntot[5, ] <- -(previous_remand) * Ntot[5, j] +
##   remand_release * Ntot[1, j] + short_release * Ntot[2, j] +
##   long_release * Ntot[3, j] + open_release * Ntot[4, j]


## ==== TPT transitions
## NOTE only for LTBI-test +ve as set
## off: stay off - handled at inflow
tpt_E[, 1] <- 0
tpt_L[, 1] <- 0
tpt_epTB[, 1] <- 0
tpt_lpTB[, 1] <- 0
## from on -> off
tpt_E[, 2] <- -E[i, 2] / tpt_drn
tpt_L[, 2] <- -L[i, 2] / tpt_drn
tpt_epTB[, 2] <- -epTB[i, 2] / tpt_drn
tpt_lpTB[, 2] <- -lpTB[i, 2] / tpt_drn
## into was from on
tpt_E[, 3] <- E[i, 2] / tpt_drn
tpt_L[, 3] <- L[i, 2] / tpt_drn
tpt_epTB[, 3] <- epTB[i, 2] / tpt_drn
tpt_lpTB[, 3] <- lpTB[i, 2] / tpt_drn

tpt_drn <- user()  #TPT duration


################## dimensions
## TB states
dim(U) <- c(NP,NT) #uninfected
dim(E) <- c(NP,NT) #early latent
dim(L) <- c(NP,NT) #late latent
dim(SD) <- c(NP,NT) #subclinical disease
dim(CD) <- c(NP,NT) #clinical disease
dim(ATT) <- c(NP,NT) #AntiTB treatment
dim(epTB) <- c(NP,NT) #early post-TB
dim(lpTB) <- c(NP,NT) #late post-TB

## extras
dim(ppop) <- c(NP, NT) # total pop
dim(rfv) <- c(NP)      #relative FOI vector

## ## test:
## dim(Ntot) <- c(NP, NT) #test
## dim(moves_Ntot) <- c(NP, NT)

## TB processes
dim(infections) <- c(NP,NT)
dim(Ereinfections) <- c(NP,NT)
dim(Lreinfections) <- c(NP,NT)
dim(Preinfections) <- c(NP,NT)
dim(reinfections) <- c(NP,NT)
dim(selfcures) <- c(NP,NT)
dim(stabilizations) <- c(NP,NT)
dim(fastprogs) <- c(NP,NT)
dim(slowprogs) <- c(NP,NT)
dim(worsens) <- c(NP,NT)
dim(relapses) <- c(NP,NT)
dim(detects) <- c(NP,NT)
dim(tbstops) <- c(NP,NT)
dim(tbmort) <- c(NP,NT)

## TB split on arrival
dim(frac_U) <- c(NP,NT) #uninfected
dim(frac_E) <- c(NP,NT) #early latent
dim(frac_L) <- c(NP,NT) #late latent
dim(frac_SD) <- c(NP,NT) #subclinical disease
dim(frac_CD) <- c(NP,NT) #clinical disease
dim(frac_ATT) <- c(NP,NT) #AntiTB treatment
dim(frac_epTB) <- c(NP,NT) #early post-TB
dim(frac_lpTB) <- c(NP,NT) #late post-TB

dim(inflow_TPTv) <- NT
dim(inflow_top) <- NP

## PPD transitions
dim(moves_U) <- c(NP,NT)
dim(moves_E) <- c(NP,NT)
dim(moves_L) <- c(NP,NT)
dim(moves_SD) <- c(NP,NT)
dim(moves_CD) <- c(NP,NT)
dim(moves_ATT) <- c(NP,NT)
dim(moves_epTB) <- c(NP,NT)
dim(moves_lpTB) <- c(NP,NT)

## TPT transitions & helpers
dim(tpt_E) <- c(NP, NT)
dim(tpt_L) <- c(NP, NT)
dim(tpt_epTB) <- c(NP, NT)
dim(tpt_lpTB) <- c(NP, NT)
dim(hrv) <- 3
dim(TPT_top) <- 3

## input arrays
dim(parm_ifrac_prevTPT) <- c(NT)
dim(parm_init_PPD) <- c(NP)

## other
dim(qolrate) <- c(NP,NT)


################## outputs
## TB notification rate within prison: NOTE does not include screening detects
output(notif100k) <- 1e5*sum(detects[1:4, 1:NT]) / (sum(ppop[1:4, 1:NT]) + 1e-2)
output(ppdpop) <- sum(ppop[1:4, 1:NT])
output(tbincout) <- (sum(fastprogs[5,1:NT]) + sum(slowprogs[5,1:NT]) + sum(relapses[5,1:NT]))*m
output(tbincall) <- (sum(fastprogs[5,1:NT]) + sum(slowprogs[5,1:NT]) + sum(relapses[5,1:NT]))*(m-1) +
  (sum(fastprogs) + sum(slowprogs) + sum(relapses))
