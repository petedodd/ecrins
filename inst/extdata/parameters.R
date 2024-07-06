## === arguments
## TB splits on arrival inflow NOTE must sum to 1
arg_tbsplits <- list(
  parm_frac_U = 1 - 1/20 - 1/10 - 1/1000 - 2/5000 - 2/100, #uninfected
  parm_frac_E = 1/20, #early latent
  parm_frac_L = 1/10, #late latent
  parm_frac_SD = 1/1000, #subclinical disease
  parm_frac_CD =  1/5000, #clinical disnease
  parm_frac_ATT = 1/5000, #AntiTB treatment
  parm_frac_epTB = 1/100, #early post-TB
  parm_frac_lpTB = 1/100 #late post-TB
)

## NOTE must sum to 1
arg_itbsplits <- list(
  parm_ifrac_U = 1 - 1/20 - 1/10 - 1/1000 - 2/5000 - 2/100, #uninfected
  parm_ifrac_E = 1/20, #early latent
  parm_ifrac_L = 1/10, #late latent
  parm_ifrac_SD = 1/1000, #subclinical disease
  parm_ifrac_CD = 1/5000, #clinical disease
  parm_ifrac_ATT = 1/5000, #AntiTB treatment
  parm_ifrac_epTB = 1/100, #early post-TB
  parm_ifrac_lpTB = 1/100 #late post-TB
)
## and array of TPT state fracs
arg_itbsplits$parm_ifrac_prevTPT <- c(1,0.0,0.0) #NOTE must sum to 1


## === PPD transitions
arg_PPDtr <- list(
  inflow = 1,          # inflow rate (less recidivists)
  remand_short = 1,    # remand -> short:   1->2
  remand_long = 1,     # remand -> long:    1->3
  remand_release = 1,  # remand -> release: 1->5
  long_short= 1,       # long -> short:     3->2
  short_release = 1,   # short -> release:  2->5
  short_open = 1,      # short -> open:     2->4
  long_release = 1,    # long -> release:   3->5
  open_release = 1,    # open -> release:   4->5
  previous_remand = 1  # previous -> remand 5->1
)

## === PPD initial stocks
arg_PPDinit <- list(parm_init_PPD=c(5:1))

## === TB hyperparms
hyperparms <- list(
  staticfoi = 1, #static or dynamic based on >0 or not
  ## --------------------------------------------------- transmission
  ## bet=list(meanlog=log(10),sdlog=0.75),         #bet,      #beta
  ptn=list(shape1=20.7,shape2=77.9),            #psi:protn Andrews
  foi=list(meanlog=log(1e-2),sdlog=0.5),    #ari0
  ## --------------------------------------------------- progression
  stb=list(meanlog=0.62, sdlog=0.068),       #kappa:arig Ragonnet
  prg=list(meanlog=-2.837,sdlog=0.32),         #eps: pp Ragonnet
  eps=list(meanlog=-6.89,sdlog=0.58),         #nu: Ragonnet
  rel=list(meanlog=-3.95,sdlog=0.27),         #omega: relapse Crampin NOTE x-ref
  ## --------------------------------------------------- detection
  CDR=list(shape1=41.80,shape2=5.22),        #UK: CDR based on WHO data & Laura Anderson data
  ## --------------------------------------------------- timescales
  drn=list(meanlog=1.1,sdlog=0.2),               #durnX log(3)
  ## --------------------------------------------------- CFRs
  txf=list(shape1=157+5,shape2= 3426-157-5), #UK data:
  ## https://www.gov.uk/government/publications/tuberculosis-in-england-2023-report-data-up-to-end-of-2022/tb-treatment-and-outcomes-england-2022
  CFR=list(shape1=25.48, shape2= 33.78),
  ## --------------------------------------------------- other
  m=list(meanlog=log(0.94),sdlog=0.1), #multiplier for population transmission = R
  tptHR = list(meanlog=-1.772,sdlog=0.089), #HR TPT protection in TBI+, Martinez et al 0.17 (0.14-0.2)
  tpt_drn = 20,       #durn of TPT: Salazar-Austin
  wsn = list(meanlog=-0.693,sdlog=0.97),  #durn AS TB D from getLNparms(0.5,1) = Frascella half of TB SC + 1ydrn
  mHR = list(meanlog=0.131,sdlog=0.071),             #post-TB mortality HR
  att_time = 0.5,      #duration of ATT
  late_post_time=2,    #duration defining early post-TB
  mort=0.02,            #mortality rate
  hrqolptb=list(meanlog=-3.324,sdlog=0.486), # HRQoL decrement while post TB
  hrqol = list(shape1 = 21.15177, shape2 = 42.36706) # GBD decrement 0.333 (0.224–0.454)
  )


qfun <- function(u,L){
  x <- NULL
  if(names(L)[1]=='meanlog') x <- qlnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape1') x <- qbeta(u,L[[1]],L[[2]])
  if(names(L)[1]=='mean') x <- qnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape') x <- qgamma(u,L[[1]],scale=L[[2]])
  if(is.null(x)) x <- u
  x
}

uv2ps <- function(u,HP,returnlist=TRUE){
  for(i in 1:length(HP)){
    if(is.list(HP[[i]])){
      u[i] <- qfun(u[i],HP[[i]])
    } else { #fixed value
      u[i] <- HP[[i]]
    }
  }
  if(returnlist){
    u <- as.list(u)
    names(u) <- names(HP)
  }
  u
}

## === tb parms
arg_tb <- uv2ps(rep(0.5, length(hyperparms)),hyperparms) # natural history

## === tb interventions + HE parms
arg_int <- list(int_time=20,                #time for SOC/BL -> INT
                disc_rate=3.5/100,            #discount rate
                LifeExp = 40,               # life expectancy
                ## === baseline/SOC
                inflow_toATT_TB0=0, # NOTE fp ATT doesn't affect state
                inflow_toATT_L0=0, #
                inflow_toATT_no0=0, #
                inflow_toTPT_TB0=0, #
                inflow_toTPT_L0=0, #
                inflow_toTPT_no0=0, #
                ## intervention
                inflow_toATT_TB1=1, # NOTE fp ATT doesn't affect state
                inflow_toATT_L1=0, #
                inflow_toATT_no1=0, #
                inflow_toTPT_TB1=0, #
                inflow_toTPT_L1=1, #
                inflow_toTPT_no1=0, #
                ## === unit costs
                uc_attppd=1, # ATT for those found passively within the system
                uc_attout=1, # ATT following release
                uc_entry_tpt_TB=1,
                uc_entry_tpt_L=1,
                uc_entry_tpt_no=1,
                uc_entry_att_TB=1,
                uc_entry_att_L=1,
                uc_entry_att_no=1,
                uc_entry_notx_TB=1,
                uc_entry_notx_L=1,
                uc_entry_notx_no=1
                ## === HRQoL
                ## hrqol = 0.333 # HRQoL decrement while CD 0.333 (0.224–0.454) NOTE now hyperparms
                ## hrqolptb=0.04, # HRQoL decrement while post TB NOTE now a hyperparm above
                ## m=1.0         #multiplier for TB events outside prison NOTE now a hyperparm above
                )                            #interventions

## === join all parm types
parms <- c(arg_tbsplits, arg_itbsplits, arg_PPDtr, arg_PPDinit, arg_tb, arg_int)
