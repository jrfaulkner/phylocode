## Trajectories used for phylodynamic simulations in paper

setwd('~/git_repos/jfaulkner_phylocode/Rev/tests/setup/')

source('phylofuncs.r')

library(phylodyn)


####-------------------------------------------
##     1) Bottleneck  (piecewise constant)
####-------------------------------------------

nsamp <- 500   #total number of samples
nstart <- 50   #number of samples at time zero
ngridB <- 101  #number of grid boundaries
samp.end <- 8  #last sample time
gend <- 8.5    #start of final (infinite) grid cell

set.seed(3)

# generate sample times
samptv_bn <- c(0, sort(runif(nsamp-nstart, 0, samp.end)) )
nsampv_bn <- c(nstart, rep(1, nsamp-nstart))

# generate coalescent times
gene_bn <- coalsim(samp_times = samptv_bn, n_sampled = nsampv_bn, 
		traj = bottleneck_trajB, lower_bound = 0.1, ne.max=1, ne.min=0.1, bstart=6, bend=4  )

# make tree object and write it
tree_bn <- generate_newick(gene_bn)
write.tree(tree_bn$newick,"../tree_only/data/bottleneck.tre",digits=10)

# write grid end
cat(gend,"\n",sep="",file="grid_ends.txt")

# write MLE constant-pop theta for prior specification
mle_bn <- mean(skyLine(summarize_phylo(tree_bn$newick))$theta)
cat(mle_bn,"\n",sep="",file="const_pop_mle.txt")

####-------------------------------------------
##     2) Boom-Bust  (Mexican hat)
####-------------------------------------------

# Function to generate trajectory
mex_hat_traj <- function(x, mf, sf,soff, sscale, snht, bloc, bscale, bht, tmx=NULL){
	if (!is.null(tmx)) {
	  if (length(x)==1) {
			  if(x > tmx) x <- tmx
		}
		if (length(x)>1) x[x>=tmx] <- tmx
	}	 
	   yg <- snht*sin((soff-x)/sscale) + bht*exp(-((x-bloc)^2)/bscale)
		 out <- mf + sf*yg
		 out
}


nsamp <- 2000
nstart <- 50
ngridB <- 101
samp.end <- 11.8
gend <- 12

set.seed(5)

# generate sample times
samptv_mx <- c(0, sort(runif(nsamp-nstart, 0, samp.end)) )
nsampv_mx <- c(nstart, rep(1, nsamp-nstart))

# generate coalescent times
gene_mx <- coalsim(samp_times = samptv_mx, n_sampled = nsampv_mx, 
	traj= mex_hat_traj, mf=0.4, sf=1, soff=5.5, sscale=3,  snht=0.25, bloc=5, bscale=.4, bht=.75  , lower_bound = 0.1 )

# make tree object and write it
tree_mx <- generate_newick(gene_mx)
write.tree(tree_mx$newick,"../tree_only/data/mex_hat.tre",digits=10)

# write grid end
cat(gend,"\n",sep="",file="grid_ends.txt",append=TRUE)

# write MLE constant-pop theta for prior specification
mle_mx <- mean(skyLine(summarize_phylo(tree_mx$newick))$theta)
cat(mle_mx,"\n",sep="",file="const_pop_mle.txt",append=TRUE)




####-------------------------------------------
##     3) Broken Exponential  (piecewise expoential)
####-------------------------------------------



nsamp <- 1000
nstart <- 100
ngridB <- 101
samp.end <- 7.8
gend <- 8

set.seed(1)

# generate sample times
samptv_be <- c(0, sort(runif(nsamp-nstart, 0, samp.end)) )
nsampv_be <- c(nstart, rep(1, nsamp-nstart))

# generate coalescent times
gene_be <- coalsim(samp_times = samptv_be, n_sampled = nsampv_be, 
			traj = piecewise_exp_traj, lower_bound = 0.1,  tbreaks=c(0, 4.5, 5, 10), lnvals=c(log(0.3),log(0.45),log(.15),log(.8))  )

# make tree object and write it
tree_be <- generate_newick(gene_be)
write.tree(tree_be$newick,"../tree_only/data/broken_exponential.tre",digits=10)

# write grid end
cat(gend,"\n",sep="",file="grid_ends.txt",append=TRUE)

# write MLE constant-pop theta for prior specification
mle_be <- mean(skyLine(summarize_phylo(tree_be$newick))$theta)
cat(mle_be,"\n",sep="",file="const_pop_mle.txt",append=TRUE)




####-------------------------------------------
##     4) Non-stationary Gaussian Process
####-------------------------------------------


# function for varying length scales
afunc <- function(xv, asmin, asmax, asig){
  mux <- mean(xv)
  d1 <- dnorm(xv, mean=mux, sd=asig)
  d1s <- d1/max(d1)
  asmin + asmax*d1s	
}

# Nonstationary covariance function from Paciorek and Schervish
cvf.matNSps <- function(x1, x2, as1, as2, nu, sig2){
  t1 <- as1^(1/4)*as2^(1/4)*(1/sqrt(0.5*(as1 + as2)) )
  hh <- abs(x1-x2)
  if (hh==0) hh <- 1e-10
  Q <- hh^2/(0.5*(as1 + as2))
  Kv <- besselK(2*sqrt(nu*Q), nu=nu)
  av <- (2*sqrt(nu*Q))^nu
  con <- (2^(1-nu))/gamma(nu)
  out <- sig2*con*t1*av*Kv
  return(out)
}

## Generate GP function on a grid of 400 
# Set up covariance function and mean
nx <- 400
xseq <- seq(0,12, length=nx)
tasvec <- afunc(xseq, 1, 15, .75) 
cmat <- matrix(0, nx, nx)
tnu <- 4
tsig2 <- 0.1 #c(0.15, 0.2)  #110 and 282
tmu <- rep(0.55, nx)  #c(1, 0.8)

for (i in 1:nx){
  for (j in i:nx){
    cmat[i,j] <- cvf.matNSps(xseq[i], xseq[j], as1=tasvec[i], as2=tasvec[j], nu=tnu, sig2=tsig2)
    if (j!=i) {
      cmat[j,i] <- cmat[i,j]
    } 
  }
}


# generate function
library(MASS)
set.seed(282)  ## this seed will define the function shape
trajGP <- rev(mvrnorm(n=1, mu=tmu, Sigma=cmat))


nsamp <- 2000
nstart <- 200
ngridB <- 101
samp.end <- 11.8
gend <- 12

samptv <- c(0, sort(runif(nsamp-nstart, 0, samp.end)) )
nsampv <- c(nstart, rep(1, nsamp-nstart))
gene <- coalsimGP(samp_times = samptv, n_sampled = nsampv, trajvec=trajGP, tvec=xseq, lower_bound=0.04)

# make tree object and write it
tree_gp <- generate_newick(gene)
write.tree(tree_gp$newick,file="../tree_only/data/GP.tre",digits=10)

# write grid end
cat(gend,"\n",sep="",file="grid_ends.txt",append=TRUE)

# write MLE constant-pop theta for prior specification
mle_gp <- mean(skyLine(summarize_phylo(tree_gp$newick))$theta)
cat(mle_gp,"\n",sep="",file="const_pop_mle.txt",append=TRUE)

