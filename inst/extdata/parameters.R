## === arguments
## TB splits on arrival inflow NOTE must sum to 1
arg_tbsplits <- list(
  parm_frac_U = 1 - 1/20 - 1/10 - 1/1000, #uninfected
  parm_frac_E = 1/20, #early latent
  parm_frac_L = 1/10, #late latent
  parm_frac_SD = 1/1000, #subclinical disease
  parm_frac_CD = 0, #clinical disease
  parm_frac_ATT = 0, #AntiTB treatment
  parm_frac_epTB = 0, #early post-TB
  parm_frac_lpTB = 0 #late post-TB
)

## TB splits in stock population NOTE TODO consider RRs by PPD type/i
## NOTE must sum to 1
arg_itbsplits <- list(
  parm_ifrac_U = 1 - 1/20 - 1/10 - 1/1000, #uninfected
  parm_ifrac_E = 1/20, #early latent
  parm_ifrac_L = 1/10, #late latent
  parm_ifrac_SD = 1/1000, #subclinical disease
  parm_ifrac_CD = 0, #clinical disease
  parm_ifrac_ATT = 0, #AntiTB treatment
  parm_ifrac_epTB = 0, #early post-TB
  parm_ifrac_lpTB = 0 #late post-TB
)
## and array of TPT state fracs
arg_itbsplits$parm_ifrac_prevTPT <- c(1,0,0) #NOTE must sum to 1


## === PPD transitions
arg_PPDtr <- list(
  inflow = 1,          # inflow rate (less recidivists)
  remand_short = 1,    # remand -> short:   1->2
  remand_long = 1,     # remand -> long:    1->3
  remand_relapse = 1,  # remand -> release: 1->5
  short_long = 1,      # short -> long:     2->3
  short_relapse = 1,   # short -> release:  2->5
  long_open = 1,       # long -> open:      3->4
  long_relapse = 1,    # long -> release:   3->5
  open_relapse = 1,    # open -> release:   4->5
  previous_remand = 1  # previous -> remand 5->1
)

## === PPD initial stocks
arg_PPDinit <- list(parm_init_PPD=c(5:1))

## === join all parm types
parms <- c(arg_tbsplits,arg_itbsplits,arg_PPDtr,arg_PPDinit)

## === TB hyperparms
hyperparms <- list(
  ## --------------------------------------------------- transmission
  bet=list(meanlog=log(10),sdlog=0.75),         #bet,      #beta
  v=list(shape1=20.7,shape2=77.9),            #psi:protn Andrews
  ari0=list(meanlog=log(3e-2),sdlog=0.75),    #ari0
  ## --------------------------------------------------- progression
  arig=list(meanlog=0.62, sdlog=0.068),       #kappa:arig Ragonnet
  pp=list(meanlog=-2.837,sdlog=0.32),         #eps: pp Ragonnet
  eps=list(meanlog=-6.89,sdlog=0.58),         #nu: Ragonnet
  rel=list(meanlog=-3.95,sdlog=0.27),         #omega: relapse Crampin NOTE x-ref
  ## --------------------------------------------------- detection
  ## CDR=list(shape1=2,shape2=2),                    #K: CDR baseline
  CDR=list(meanlogit=0,sd=0.3),                    #K: CDR baseline
  ## --------------------------------------------------- timescales
  drnX=list(meanlog=1.1,sdlog=0.2),               #durnX log(3)
  ## --------------------------------------------------- CFRs
  txf=list(shape1=2.71,shape2= 87.55),
  cfrn=list(shape1=25.48, shape2= 33.78),
  )

qfun <- function(u,L){
  if(names(L)[1]=='meanlog') x <- qlnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape1') x <- qbeta(u,L[[1]],L[[2]])
  if(names(L)[1]=='mean') x <- qnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape') x <- qgamma(u,L[[1]],scale=L[[2]])
  x
}

uv2ps <- function(u){
  for(i in 1:length(hyperparms)){
    if(is.list(hyperparms[[i]])){
      u[i] <- qfun(u[i],hyperparms[[i]])
    } else { #fixed value
      u[i] <- hyperparms[[i]]
    }
  }
  u
}
