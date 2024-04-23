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

## other parameters
att_time <- user(0.5)     #duration of ATT
late_post_time <- user(2) #duration defining early post-TB
mort <- user(0.02)        #mortality rate

## economic parameters
disc_rate <- user(0.03)



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


## initial PPD split TODO to consider interaction with TB
parm_init_PPD[] <- user()

################## initial state
initial(U[1:NP,1:NT]) <- parm_ifrac_U * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##uninfected
initial(E[1:NP,1:NT]) <- parm_ifrac_E * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##early latent
initial(L[1:NP,1:NT]) <- parm_ifrac_L * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##late latent
initial(SD[1:NP,1:NT]) <- parm_ifrac_SD * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##subclinical disease
initial(CD[1:NP,1:NT]) <- parm_ifrac_CD * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##clinical disease
initial(ATT[1:NP,1:NT]) <- parm_ifrac_ATT * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##AntiTB treatment
initial(epTB[1:NP,1:NT]) <- parm_ifrac_epTB * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##early post-TB
initial(lpTB[1:NP,1:NT]) <- parm_ifrac_lpTB * parm_init_PPD[i] * parm_ifrac_prevTPT[j] ##late post-TB

## economic states
initial(CC0) <- 0 #undiscounted cumulative costs
initial(CC) <- 0 #discounted cumulative costs
initial(cATTtp) <- 0 #ATT true positive counter
initial(cATTfp) <- 0 #ATT false positive counter
initial(cTPT) <- 0   #cumulative TPT counter

################## dynamics
## TODO currently no TPT in

## === TB states
## uninfected
deriv(U[,]) <-  - infections[i,j] - mort*U[i,j] +
  inflow * frac_U[i,j] + moves_U[i,j]
## early latent
deriv(E[,]) <- infections[i,j] + reinfections[i,j] + selfcures[i,j] - stabilizations[i,j] - fastprogs[i,j] - mort*E[i,j] +
  inflow * frac_E[i,j] + moves_E[i,j]
## late latent
deriv(L[,]) <- - Lreinfections[i,j] + stabilizations[i,j] - slowprogs[i,j] - mort*L[i,j] +
  inflow * frac_L[i,j] + moves_L[i,j]
## subclinical disease
deriv(SD[,]) <- fastprogs[i,j] + slowprogs[i,j] - selfcures[i,j]- mort*SD[i,j] - worsens[i,j] +
  inflow * frac_SD[i,j] + moves_SD[i,j]
## clinical disease
deriv(CD[,]) <- worsens[i,j]- mort*CD[i,j] + relapses[i,j] +
  inflow * frac_CD[i,j] + moves_CD[i,j]
## AntiTB treatment
deriv(ATT[,]) <- detects[i,j] * CD[i,j] - ATT[i,j]/att_time - mort*ATT[i,j] +
  inflow * frac_ATT[i,j] + moves_ATT[i,j]
## early post-TB TODO consider treatment outcomes
deriv(epTB[,]) <- ATT[i,j]/att_time - relapses[i,j] - epTB[i,j]/late_post_time - mort*epTB[i,j] +
  inflow * frac_epTB[i,j] + moves_epTB[i,j]
## late post-TB
deriv(lpTB[,]) <- epTB[i,j]/late_post_time - mort*lpTB[i,j] - Preinfections[i,j] +
  inflow * frac_lpTB[i,j] + moves_lpTB[i,j]

## economic states & counters
deriv(CC0) <- 0
deriv(CC) <- 0 - disc_rate * CC
deriv(cATTtp) <- 0 #ATT true positive counter
deriv(cATTfp) <- 0 #ATT false positive counter
deriv(cTPT) <- 0   #cumulative TPT counter


## === TB processes TODO write out & parametrize these based on other models
infections[,] <- 0
Ereinfections[,] <- 0
Lreinfections[,] <- 0
Preinfections[,] <- 0
reinfections[,] <- Ereinfections[i,j] + Lreinfections[i,j] + Preinfections[i,j]
selfcures[,] <- 0
stabilizations[,] <- 0
fastprogs[,] <- 0
slowprogs[,] <- 0
selfcures[,] <- 0
worsens[,] <- 0
relapses[,] <- 0
detects[,] <- 0



## screening rates:
##     attributes Uninfected, LTBI, TB, Previous
## outcomes:
##     ATT, TPT, no change
## probabilities within inflow:
inflow_toATT_TB <- user(0) #NOTE fp ATT doesn't affect state
inflow_toATT_P <- user(0) #NOTE fp ATT doesn't affect state - just for counting
inflow_toATT_rest <- user(0) #NOTE fp ATT doesn't affect state - just for counting
inflow_toTPT_L <- user(0)    #assume perfect spec
## vector form for compactness
inflow_TPTv[1] <- (1-inflow_toTPT_L)
inflow_TPTv[2] <- (inflow_toTPT_L)
inflow_TPTv[3] <- 0
## inflow to remand
inflow_top[1] <- 1
inflow_top[2:NP] <- 0


## === TB split on arrival, including TPT
frac_U[,] <- parm_frac_U * inflow_top[i] #uninfected
frac_E[,] <- parm_frac_E * inflow_top[i] * inflow_TPTv[j] #early latent
frac_L[,] <- parm_frac_L * inflow_top[i] * inflow_TPTv[j] #late latent
frac_SD[,] <- parm_frac_SD * inflow_top[i] * (1-inflow_toATT_TB) #subclinical disease
frac_CD[,] <- parm_frac_CD * inflow_top[i] * (1-inflow_toATT_TB) #clinical disease
frac_ATT[,] <- (parm_frac_ATT + (parm_frac_SD + parm_frac_CD) * inflow_toATT_TB) * inflow_top[i] #ATT
frac_epTB[,] <- parm_frac_epTB * inflow_top[i] * inflow_TPTv[j] #early post-TB
frac_lpTB[,] <- parm_frac_lpTB * inflow_top[i] * inflow_TPTv[j] #late post-TB


## === PPD processes: transitions in 1st index
## TODO check these transitions
moverate[,] <- 0
## remand -> short:   1->2
moverate[1,2] <- -remand_short
moverate[2,1] <- remand_short
## remand -> long:    1->3
moverate[1,3] <- -remand_long
moverate[3,1] <- remand_long
## remand -> release: 1->5
moverate[1,5] <- -remand_release
moverate[5,1] <- remand_release
## long -> short:     3->2
moverate[3,2] <- -long_short
moverate[2,3] <- long_short
## short -> release:  2->5
moverate[2,5] <- -short_release
moverate[5,2] <- short_release
## short -> open:      2->4
moverate[2,4] <- -short_open
moverate[4,2] <- short_open
## long -> release:   3->5
moverate[3,5] <- -long_release
moverate[5,3] <- long_release
## open -> release:   4->5
moverate[4,5] <- -open_release
moverate[5,4] <- open_release
## previous -> remand 5->1
moverate[5,1] <- -previous_remand
moverate[1,5] <- previous_remand

## ## NOTE check syntax has worked
## print("moverate 11 :{moverate[1,1]}")
## print("moverate 12 :{moverate[1,2]}")

## dX_{ij}/dt = \sum_k R_{ik} X_{kj}
## matrix multiplications
moves_Ut[,,] <- moverate[i,j] * U[j,k] #temp
moves_U[,] <- sum(moves_Ut[i,,j])      #matrix mult
moves_Et[,,] <- moverate[i,j] * E[j,k] #temp
moves_E[,] <- sum(moves_Et[i,,j])      #matrix mult
moves_Lt[,,] <- moverate[i,j] * L[j,k] #temp
moves_L[,] <- sum(moves_Lt[i,,j])      #matrix mult
moves_SDt[,,] <- moverate[i,j] * SD[j,k] #temp
moves_SD[,] <- sum(moves_SDt[i,,j])      #matrix mult
moves_CDt[,,] <- moverate[i,j] * CD[j,k] #temp
moves_CD[,] <- sum(moves_CDt[i,,j])      #matrix mult
moves_ATTt[,,] <- moverate[i,j] * ATT[j,k] #temp
moves_ATT[,] <- sum(moves_ATTt[i,,j])      #matrix mult
moves_epTBt[,,] <- moverate[i,j] * epTB[j,k] #temp
moves_epTB[,] <- sum(moves_epTBt[i,,j])      #matrix mult
moves_lpTBt[,,] <- moverate[i,j] * lpTB[j,k] #temp
moves_lpTB[,] <- sum(moves_lpTBt[i,,j])      #matrix mult


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
dim(moverate) <- c(NP,NP)
dim(moves_Ut) <- c(NP,NP,NT)
dim(moves_U) <- c(NP,NT)
dim(moves_Et) <- c(NP,NP,NT)
dim(moves_E) <- c(NP,NT)
dim(moves_Lt) <- c(NP,NP,NT)
dim(moves_L) <- c(NP,NT)
dim(moves_SDt) <- c(NP,NP,NT)
dim(moves_SD) <- c(NP,NT)
dim(moves_CDt) <- c(NP,NP,NT)
dim(moves_CD) <- c(NP,NT)
dim(moves_ATTt) <- c(NP,NP,NT)
dim(moves_ATT) <- c(NP,NT)
dim(moves_epTBt) <- c(NP,NP,NT)
dim(moves_epTB) <- c(NP,NT)
dim(moves_lpTBt) <- c(NP,NP,NT)
dim(moves_lpTB) <- c(NP,NT)

## input arrays
dim(parm_ifrac_prevTPT) <- c(NT)
dim(parm_init_PPD) <- c(NP)

################## outputs
## TODO
