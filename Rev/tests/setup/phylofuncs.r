#####


###  need classic skyline methods for hyperparameter estimation for prior selection

## isochronous classic skyline

# skyLine <- function(ctimes) {
#  tk <- ctimes
#  uk <- diff(tk)
#  nk <- length(tk)
#  out <- numeric(nk)
#  nkv <- nk:2
#  nkv*(nkv-1)*uk/2
# }

## heterochronous classic skyline
# -- this should take a phylo.summary object as produced in phylodyn
# theta.hat.k = sum{j=0,j=j_k} n_{k,j}*(n_{k,j}-1)*w_{k,j}/2
# j_k in {0,...,n-1} ; w_k = (w_k0, ... , w_{kj_k})
# sum{j=0,j_k} w_{k,j} = u_k - so w_{k,j} are subinterval lengths

skyLine <- function(phy){
	# classic skyline for iso or heterochronous data.
	# phy is a list from summarize_phylo() with 
	# coal_times, samp_times, and n_sampled
	# outputs vector of Ne estimates sorted from most recent.
	st <- phy$samp_times 
	ns <- phy$n_sampled 
	ct <- phy$coal_times 
	sm <- data.frame(type='s', time=st, nadd=ns)
	cm <- data.frame(type='c', time=ct, nadd=-1)
	zm <- merge(sm, cm, all=T)
	zm <- zm[order(zm$time), ]
	zm$ncount <- cumsum(zm$nadd)
	zn <- nrow(zm)
	zuid <- numeric(zn)
	zuid[zm$type=='c'] <- 1:nrow(cm)
	for (j in (zn-1):1) {
		if (zuid[j]==0) zuid[j] <- zuid[j+1]
	}
	zm$zuid <- zuid
	wk <- diff(zm$time)
	dm <- data.frame(uid=zm$zuid[-1], ncount=zm$ncount[-zn], wk=round(wk,8) )
	dm$thsub <- dm$ncount*(dm$ncount-1)*dm$wk/2
	adm <- aggregate(dm[,"thsub"], by=list(uid=dm$uid), sum)
	zc1 <- c(0,cm$time[-nrow(cm)])
	zc2 <- cm$time
	ctmid <- zc1 + (zc2-zc1)/2
	
	return(list(k=adm$uid, mid_ctime=ctmid, theta=adm$x))
}


wtMean <- function(yvec, dvec){
  wv <- 1/dvec
  sw <- sum(wv)
  rwv <- wv/sw
  sum(rwv*yvec)
}

wtVar <- function(yvec, dvec){
  wv <- 1/dvec
  sw <- sum(wv)
  rwv <- wv/sw
  mbar <- sum(rwv*yvec)
  vhat <- (1/length(yvec))*sum( (yvec - mbar)^2/dvec )
  vhat
}


varRef <- function(nn, kap, omg2, order=1){
	#nn is number of thetas, kap = 1/gam^2, omg2=omega^2
	#returns vector of marginal variances
	nnv <- 1:nn
	if (order==1) out <- omg2 + (nnv-1)*(1/kap)
	if (order==2) out <- omg2 + (1/kap)*nnv*(nnv-1)*(2*nnv-1)/6
	return(out)
}


set_zeta <- function(phylo, ncell, alpha=0.05, order=1){
	zsky <- skyLine(phylo)
	vld <- var(log(zsky$theta))
	if (order==1)	cmr <- varRef(ncell, 1,  vld, order=1)
	if (order==2)	cmr <- varRef(ncell, 1,  vld, order=2)
	sref <- exp(mean(0.5*log(cmr)))
  uu <- sqrt(vld)
  zz <- uu/(sref*(tan((pi/2)*(1-alpha))))
  zz
}



make_grid <- function(coal_times, samp_times, Ngrid){
	gbds <- range(c(coal_times, samp_times))
	grd <- seq(gbds[1],gbds[2],length=Ngrid)
	grd
}


bottleneck_trajB <- function(t, ne.max=100, ne.min=10, bstart=2, bend=1,...) {
  result = rep(0,length(t))
  result[t <= bend] <- ne.max
  result[t > bend & t < bstart] <- ne.min
  result[t >= bstart] <- ne.max
  return(result)
}

boombust_trajB <- function(t, bust=1, scale=1000, rate=1,...) {
  result = rep(0, length(t))
  result[t <= bust] = scale*exp(rate*(t[t <= bust]-bust))
  result[t >  bust] = scale*exp(rate*(bust-t[t >  bust]))
  return(result)
}

logistic_trajB <- function(t, ne.max=100, ne.min=10, offset=0, a=2,...) {
  t = t + offset
  result = rep(0, length(t))
  result[(t %% 12) <= 6] = ne.min + (ne.max-ne.min)/(1+exp((3-(t[(t %% 12) <= 6] %% 12)) * a))
  result[(t %% 12) >  6] = ne.min + (ne.max-ne.min)/(1+exp(((t[(t %% 12) >  6] %% 12) - 12 + 3) * a))
  return(result)
}


## !! NEED to make these work for scalar or vector t's, and need to account for constant traj for t > tmax

piecewise_exp_traj <- function(t, tbreaks, lnvals, reverse=FALSE, tmx=NULL){
	# tbreaks must include first and last times of t
	# lnvals are log Ne values at sbreaks, tmx is max value of t, above which traj is constant
	ntb <- length(tbreaks)
	lnsv <- numeric(length(t))
	lnsv[1] <- lnvals[1]
	if (!is.null(tmx)) {
	  if (length(t)==1) {
			  if(t > tmx) t <- tmx
		}
		if (length(t)>1) t[t>=tmx] <- tmx
	}	 
	for (j in 1:(ntb-1)) {
 	 b1 <- (lnvals[j+1]-lnvals[j])/(tbreaks[j+1]-tbreaks[j])
 	 b0 <- lnvals[j]-b1*tbreaks[j]
 	 segid <- t > tbreaks[j] & t <= tbreaks[j+1]
 	 lnsv[segid] <- b0 + b1*t[segid]
	}
	out <- exp(lnsv)
	if (reverse==TRUE) out <- rev(out)
	return(out)
}

#' Simulate from inhomogeneous, heterochronous coalescent.
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param upper numeric upper limit on \code{traj} function on its support.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim(0:2, 3:1, unif_traj, lower_bound=10, level=10)
coalsim <- function(samp_times, n_sampled, traj, lower_bound, ...)
{
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower_bound/traj(time, ...))
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
         intercoal_times = c(coal_times[1], diff(coal_times)),
         samp_times = samp_times, n_sampled = n_sampled))
}



#' Simulate from inhomogeneous, heterochronous coalescent.
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param upper numeric upper limit on \code{traj} function on its support.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim(0:2, 3:1, unif_traj, lower_bound=10, level=10)
coalsimGP <- function(samp_times, n_sampled, trajvec, tvec, lower_bound, ...)
{
	# trajvec has gp trajectory (should be dense)
	# tvec is grid boundaries for trajvec (length(trajvec) + 1 )
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
	
	tup <- tvec[-1]
	tlo <- tvec[-length(tvec)]
	tint <- tvec[2]-tvec[1]
	tmids <- tvec[-1]-tint/2
	tadj <- diff(range(tvec))/(10*(length(tvec)-1))

  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
		# convert time to index of trajvec
		dtime <- time
 	 if (sum(time==tlo)>0) dtime <- time + tadj  #adjujst for cell boundaries
	 if (time >= max(tvec)) dtime <- max(tvec) - tadj	 # set long times equal to max
 	 tind <- which(dtime > tlo & dtime < tup)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower_bound/trajvec[tind] )
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
         intercoal_times = c(coal_times[1], diff(coal_times)),
         samp_times = samp_times, n_sampled = n_sampled))
}







## log likelihood for a constant Ne
const_coal_loglike <- function(theta, y, C, D){
	 -theta*sum(y) - exp(-theta)*sum(C*D)
} 

const_coal_loglikeB <- function(theta, y, C, D){
	 -theta*y - exp(-theta)*C*D
} 


## This should create data list for present to stan
make_coal_data <- function(samp_times, n_sampled, coal_times, grid)
{
  ns <- length(samp_times)
  nc <- length(coal_times)
  ng <- length(grid)-1
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
  
  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")
  
  t <- sort(unique(c(samp_times, coal_times, grid)))  #combined times
  alin <- rep(0, length(t))   #number of active lineages
  
  for (i in 1:ns)
    alin[t >= samp_times[i]] <- alin[t >= samp_times[i]] + n_sampled[i]
  
  for (i in 1:nc)
    alin[t >= coal_times[i]] <- alin[t >= coal_times[i]] - 1
  
  #print(l)
  
  if (sum((alin < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")
  
  mask <- alin > 0
  t <- t[mask]  # drop times with zero lineages
  alin <- head(alin[mask], -1) #drops first
  
  gridrep <- rep(0, ng)  #numbers of subintervals within a grid cell
  for (i in 1:ng)
    gridrep[i] <- sum(t > grid[i] & t <= grid[i+1])
  
  C <- 0.5 * alin * (alin-1)  #binomial coefficient for active lineages
  D <- diff(t)  #time step width
  
  y <- rep(0, length(D))  #indicator for coal event within sub interval
  y[t[-1] %in% coal_times] <- 1
  ny <- length(y)
	ncoal <- sum(y)
  
  rep_idx <- cumsum(gridrep)
  rep_idx <- cbind(start=rep_idx-gridrep+1,end=rep_idx)
	
	coalind <- integer(ny)
	cnt <- 1
	for (j in 1:ny){
		coalind[j] <- cnt
		if (y[j]==1) cnt <- cnt+1
	}
	cstart <- integer(ncoal)
	cend <- integer(ncoal)
	for (k in 1:ncoal){
		tmpi <- which(coalind==k)
		cstart[k] <- min(tmpi)
		cend[k] <- max(tmpi)
	}
	grend <- cumsum(gridrep) 
	grstart <- c(1, grend[-ny])
  mle <- optimize(f=const_coal_loglike, interval=c(log(.00000001), log(1E6)), y=y, 
		C=C, D=D, maximum=TRUE)$max
	# for (j in 1:ng){
  # 		  if (j==1) idx <- rep(j,gridrep[j])
  # 			idx <- c(idx, rep(j,gridrep[j]))
  # 	}
  ncoalv <- numeric(ng)
	for (jj in 1:ng){
		ncoalv[jj] <- sum(y[rep_idx[jj,1]:rep_idx[jj,2]])
	}
  
  return(list(J=ng, N=ny, y=y, gridrep=gridrep, Aik=C, coalind=coalind,ncoal=ncoal,ncoalv=ncoalv,cstart=cstart,
		cend=cend,gridstart=grstart,gridend=grend, dalpha=D, rep.idx=rep_idx, log_mu=mle))
}


make_coal_samp_data <- function(samp_times, n_sampled, coal_times, grid)
{
  ns <- length(samp_times)
  nc <- length(coal_times)
  ng <- length(grid)-1
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
  
  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")
  
  t <- sort(unique(c(samp_times, coal_times, grid)))
  alin <- rep(0, length(t))
  
  for (i in 1:ns)
    alin[t >= samp_times[i]] <- alin[t >= samp_times[i]] + n_sampled[i]
  
  for (i in 1:nc)
    alin[t >= coal_times[i]] <- alin[t >= coal_times[i]] - 1
  
  #print(l)
  
  if (sum((alin < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")
  
  mask <- alin > 0
  t <- t[mask]
  alin <- head(alin[mask], -1)
  
  gridrep <- rep(0, ng)
  for (i in 1:ng)
    gridrep[i] <- sum(t > grid[i] & t <= grid[i+1])
  
  C <- 0.5 * alin * (alin-1)
  D <- diff(t)
  
  y <- rep(0, length(D))
  y[t[-1] %in% coal_times] <- 1
  ny <- length(y)
	ncoal <- sum(y)
  
  buckets <- cut(x = samp_times, breaks = t,
                include.lowest = TRUE)
  tab <- aggregate(n_sampled ~ buckets, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$buckets)] <- tab$n_sampled
  count[head(t, -1) >= max(samp_times)] <- NA
  scount <- count[!is.na(count)]
  #max(which(!is.na(count)))
  #can there be na's other than at the end?
  smax <- length(scount)

  gmat <- cbind(grid[-length(grid)], grid[-1])
  b1len <- which(max(samp_times) > gmat[,1] & max(samp_times) <= gmat[,2]) #gives length of b1
  nb1sub <- sum(gridrep[1:b1len])

  rep_idx <- cumsum(gridrep)
  rep_idx <- cbind(rep_idx-gridrep+1,rep_idx)
	
	coalind <- integer(ny)
	cnt <- 1
	for (j in 1:ny){
		coalind[j] <- cnt
		if (y[j]==1) cnt <- cnt+1
	}
	cstart <- integer(ncoal)
	cend <- integer(ncoal)
	for (k in 1:ncoal){
		tmpi <- which(coalind==k)
		cstart[k] <- min(tmpi)
		cend[k] <- max(tmpi)
	}
	grend <- cumsum(gridrep) 
	grstart <- c(1, grend[-ny])

  mle <- optimize(f=const_coal_loglike, interval=c(log(1), log(1E6)), y=y, 
		C=C, D=D, maximum=TRUE)$max
  
 # return(list(t=t, alin=alin, C=C, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, 
 # 				ng=ng, rep_idx=rep_idx, samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid))
 return(list(J=ng, N=ny, y=y, gridrep=gridrep, Aik=C, dalpha=D, nsamps=ns, scount=scount, coalind=coalind, ncoal=ncoal,
	   cstart=cstart, cend=cend, gridstart=grstart,gridend=grend, smax=smax, Jb1 = b1len, Nb1 = nb1sub, log_mu=mle ))

}




###  THIS WILL BE PART OF STAN MODEL
####  maybe as incr log prob function

zcoal_loglik <- function(init, f, grad=FALSE)
{
  if (init$J != length(f))
    stop(paste("Incorrect length for f; should be", init$J))
  fg = rep(f, init$gridrep)
  llnocoal = init$dalpha * init$Aik * exp(-fg)
  lls = - init$y * fg - llnocoal
  ll = sum(lls[!is.nan(lls)])
    
}

zcoal_samp_loglik <- function(init, f, beta0, beta1)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  f = rep(f, init$gridrep)
  llnocoal = init$D * init$C * exp(-f)
  lls = - init$y * f - llnocoal
  #print(lls)
  #print(init$count)
  llsampevents = beta1 * init$count * f
  #print(llsampevents[!is.na(init$count)])
  llsampnoevents = init$D * beta0 * exp(f)^beta1
  #print(llsampnoevents[!is.na(init$count)])
  llsamp = init$ns * log(beta0) + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  llcoal = sum(lls[!is.nan(lls)])
  return(llcoal + llsamp)
}


## squared exponential cov fun for ihPois_GP
expCovih <- function(ssv, pvec){
  # plist is list of parms {sigmaf2, ll}
  sigmaf2 <- pvec[1]
  ll <- pvec[2]
  nn <- length(ssv)
  cmat <- matrix(0, nn, nn)
  for (i in 1:nn){
    for (j in 1:nn){
      cmat[i,j] <- sigmaf2*exp(-0.5*((i-j)/ll)^2)
    }
  }
  cmat
}

muGPconst <- function(ssv, cmu=0){
  rep(cmu, length(ssv))
}


### Functions for generating samples
###  from inhomogeneous Poisson process via thinning.
###  Follows Adams et al (2009)

ihPois_GP <- function(tvec, lamb.max, mu.fun, cov.fun, theta.mu, theta.cov){
   ## GP is generated within this function via GP{mu.fun(sv, theta.mu), cov.fun(sv, theta.cov)}
   ## tvec is the vector (min:max) of space of interest
   ## lamb.max is max Poisson rate
   ## returns locations of thinned (inhomogeneous) events
   library(MASS)
   
   mu.t <- max(tvec) - min(tvec) #measure on T
   #cat("mu.t =", mu.t, "\n")
   J <- rpois(1, mu.t*lamb.max)  # number of events for homogeneous
   #cat("J =", J, "\n")
   sv <- sort(runif(J, min(tvec), max(tvec))) # uniformly distributed and ordered event locations
   #cat("sv =", sv, "\n")
   gmu <- mu.fun(sv, theta.mu)  #get mean of GP on S
   gcov <- cov.fun(sv, theta.cov)  # get cov of GP on S
   gs <- as.vector(mvrnorm(n=1, mu=gmu, Sigma=gcov)) # get GP at sample locations
   #cat("gs = ", gs)
   epi <- NULL  #initialize set of accepted events
   for (j in 1:J){
     rj <- runif(1)  # acceptance prob
	 #cat("j=",j, "rj = ", rj, "\n")
	 if (rj < plogis(gs[j])) epi <- c(epi, sv[j]) # select events
   }
 return(list(times = epi, gpfun = gs))
} 


# testih <- ihPois_GP(tvec = 0:100, lamb.max=5, mu.fun=muGPconst, cov.fun= expCovih, 
					# theta.mu=0, theta.cov=c(0.1, 50) )

					
# pdf("output/plots/test_plots_poissamp.pdf")
  # hist(testih$gpfun, nclass=50)
  # plot(1:length(testih$gpfun), testih$gpfun, type='l') 
# dev.off()

psamp_traj <- function(trng, b0, b1, traj, ...){
   ## The ... are for the traj function.
   ## Traj is a function for Ne 
   ## lamb.max is max Poisson rate
   
   tts <- seq(trng[1], trng[2], length=500)
   lamb.max <- max(b0*traj(tts, ...)^b1)
   mu.t <- trng[2] - trng[1] #measure on T
   #cat("mu.t =", mu.t, "\n")
   J <- rpois(1, mu.t*lamb.max)  # number of events for homogeneous
   #cat("J =", J, "\n")
   sv <- sort(runif(J, trng[1], trng[2])) # uniformly distributed and ordered event locations
   #cat("sv =", sv, "\n")
   lamb.t <- b0*traj(sv, ...)^b1
   lratio <- lamb.t/lamb.max
   epi <- sv[runif(J) < lratio]

 return(list(times = epi, lamb.t = lamb.t))

} 

piecewise_const_traj <- function(tvec, levs, breaks,...){
  ## assumes num breaks = num levs - 1
  ## assumes first lev lies on interval tvec[1]:breaks[1]
  ## and last lev is on interval breaks[nb]:tvec[nt]
  nt <- length(tvec)
  nl <- length(levs)
  outv <- numeric(nt)
  for (i in 1:nl){
    if (i==1) bint <- tvec < breaks[1]
	if (i==nl) bint <- tvec >= breaks[i-1]
	if (i > 1 & i < nl) bint <- tvec >= breaks[i-1] & tvec < breaks[i]
	outv[bint] <- levs[i]
  }
  return(outv)
}


b1_logistic_trajB <- function(t, b1.max=1.5, b1.min=0.5, offset=0, a=2,...) {
  t = t + offset
  result = rep(0, length(t))
  result[(t %% 12) <= 6] = b1.min + (b1.max-b1.min)/(1+exp((3-(t[(t %% 12) <= 6] %% 12)) * a))
  result[(t %% 12) >  6] = b1.min + (b1.max-b1.min)/(1+exp(((t[(t %% 12) >  6] %% 12) - 12 + 3) * a))
  return(result)
}

b1_linear_trajB <- function(t, a0t, a1t, ...){
	a0t + a1t*t
}

b1_sine_trajB <- function(t, a0t, a1t, kt, ...){
	a0t + a1t*sin(t/kt)
}

#piecewise_const_traj(tvec=1:20, levs=c(0,5,1,5), breaks=c(5, 10, 15) )

### NEEDS EDITING BELOW

psamp_fulltraj <- function(trng, b0, b1traj, traj, ...){
   ## The ... are for both traj functions.
   ## Traj is a function for Ne 
   ## lamb.max is max Poisson rate
   
   tts <- seq(trng[1], trng[2], length=2000)
   b1s <- b1traj(tts, ...)
   lamb.max <- max(b0*traj(tts, ...)^b1s)
   mu.t <- trng[2] - trng[1] #measure on T
   #cat("mu.t =", mu.t, "\n")
   J <- rpois(1, mu.t*lamb.max)  # number of events for homogeneous
   #cat("J =", J, "\n")
   sv <- sort(runif(J, trng[1], trng[2])) # uniformly distributed and ordered event locations
   #cat("sv =", sv, "\n")
   lamb.t <- b0*traj(sv, ...)^b1traj(sv, ...)
   lratio <- lamb.t/lamb.max
   epi <- sv[runif(J) < lratio]

 return(list(times = epi, lamb.t = lamb.t))

} 


est_alpha <- function(xv, b1r){
	a1 <- (b1r[2]-b1r[1])/(max(xv)-min(xv))
	a0 <- b1r[2] - a1*max(xv)
	return(list(a0=a0, a1=a1))
}


est_logb0 <- function(xv, yv, b1r){
	av <- est_alpha(xv, b1r)
	bx <- av$a0 + av$a1*xv
	lb0 <- mean(log(yv) - bx*log(xv))
	return(lb0) 
}

est_b1 <- function(xv, yv, b1r){
	av <- est_alpha(xv, b1r)
	bx <- av$a0 + av$a1*xv
	return(bx)
}

quan_ssq <- function(parmv, qset){
	qq1 <- qnorm(0.025, mean=parmv[1], sd=exp(parmv[2]))
	qq2 <- qnorm(0.975, mean=parmv[1], sd=exp(parmv[2]) )
	ssq <- (qq1-qset[1])^2 + (qq2-qset[2])^2
	return(ssq) 
}
set_logb0_prior <- function(xv, yv, blim){
	q1 <- est_logb0(xv, yv, blim)
	q2 <- est_logb0(xv, yv, c(0.95,1.05))
	qset <- c(q1, q2)
	pint <- c(-2, 0.5)
	bf <- optim(par=pint, fn=quan_ssq, qset=qset)
	return(list(logb0mean=bf$par[1], logb0sd=exp(bf$par[2]) ))
}




extract_theta <- function(mfit, obstype="normal",  alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.th <- rstan::extract(mfit, "theta")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.th1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.th1 <- rbind(tmp.th1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.th1)[[2]], pattern="theta")
	  zth <-  grep(x=dimnames(tmp.th1)[[2]], pattern="ztheta")
	  thind <- setdiff(ath, zth)
	  tmp.th <- tmp.th1[ , thind]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="theta")
  	zth <-  grep(x=colnames(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- mfit[ , thind]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="theta")
  	zth <-  grep(x=names(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- as.matrix(mfit[ , thind])
  }
  plow <- alpha/2
  phigh <- 1 - alpha/2

  if (obstype=="normal"){
   tmp.md <- apply(tmp.th, 2, median)
   tmp.l <- apply(tmp.th, 2, quantile, probs=plow)
   tmp.u <- apply(tmp.th, 2, quantile, probs=phigh)
  }
  if (obstype=="poisson"){
   tmp.md <- exp(apply(tmp.th, 2, median))
   tmp.l <- exp(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  }
  if (obstype=="binomial"){
   tmp.md <- plogis(apply(tmp.th, 2, median))
   tmp.l <- plogis(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- plogis(apply(tmp.th, 2, quantile, probs=phigh))
  }
   out <- list(postmed=tmp.md, bci.lower=tmp.l, bci.upper=tmp.u)
   out

}



extract_theta_A <- function(mfit, obstype="normal",  alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.th <- rstan::extract(mfit, "theta_A")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.th1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.th1 <- rbind(tmp.th1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.th1)[[2]], pattern="theta_A")
	  zth <-  grep(x=dimnames(tmp.th1)[[2]], pattern="ztheta")
	  thind <- setdiff(ath, zth)
	  tmp.th <- tmp.th1[ , thind]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="theta_A")
  	zth <-  grep(x=colnames(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- mfit[ , thind]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="theta_A")
  	zth <-  grep(x=names(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- as.matrix(mfit[ , thind])
  }
  plow <- alpha/2
  phigh <- 1 - alpha/2

  if (obstype=="normal"){
   tmp.md <- apply(tmp.th, 2, median)
   tmp.l <- apply(tmp.th, 2, quantile, probs=plow)
   tmp.u <- apply(tmp.th, 2, quantile, probs=phigh)
  }
  if (obstype=="poisson"){
   tmp.md <- exp(apply(tmp.th, 2, median))
   tmp.l <- exp(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  }
  if (obstype=="binomial"){
   tmp.md <- plogis(apply(tmp.th, 2, median))
   tmp.l <- plogis(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- plogis(apply(tmp.th, 2, quantile, probs=phigh))
  }
   out <- list(postmed=tmp.md, bci.lower=tmp.l, bci.upper=tmp.u)
   out

}





extract_theta_B <- function(mfit, obstype="normal",  alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.th <- rstan::extract(mfit, "theta_B")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.th1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.th1 <- rbind(tmp.th1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.th1)[[2]], pattern="theta_B")
	  zth <-  grep(x=dimnames(tmp.th1)[[2]], pattern="ztheta")
	  thind <- setdiff(ath, zth)
	  tmp.th <- tmp.th1[ , thind]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="theta_B")
  	zth <-  grep(x=colnames(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- mfit[ , thind]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="theta_B")
  	zth <-  grep(x=names(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- as.matrix(mfit[ , thind])
  }
  plow <- alpha/2
  phigh <- 1 - alpha/2

  if (obstype=="normal"){
   tmp.md <- apply(tmp.th, 2, median)
   tmp.l <- apply(tmp.th, 2, quantile, probs=plow)
   tmp.u <- apply(tmp.th, 2, quantile, probs=phigh)
  }
  if (obstype=="poisson"){
   tmp.md <- exp(apply(tmp.th, 2, median))
   tmp.l <- exp(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  }
  if (obstype=="binomial"){
   tmp.md <- plogis(apply(tmp.th, 2, median))
   tmp.l <- plogis(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- plogis(apply(tmp.th, 2, quantile, probs=phigh))
  }
   out <- list(postmed=tmp.md, bci.lower=tmp.l, bci.upper=tmp.u)
   out

}







extract_beta1 <- function(mfit, obstype="normal",  alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.th <- rstan::extract(mfit, "beta1")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.th1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.th1 <- rbind(tmp.th1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.th1)[[2]], pattern="beta1")
	  zth <-  grep(x=dimnames(tmp.th1)[[2]], pattern="zbeta1")
	  thind <- setdiff(ath, zth)
	  tmp.th <- tmp.th1[ , thind]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="beta1")
  	zth <-  grep(x=colnames(mfit), pattern="zbeta1")
  	thind <- setdiff(ath, zth)
  	tmp.th <- mfit[ , thind]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="beta1")
  	zth <-  grep(x=names(mfit), pattern="zbeta1")
  	thind <- setdiff(ath, zth)
  	tmp.th <- as.matrix(mfit[ , thind])
  }
  plow <- alpha/2
  phigh <- 1 - alpha/2

  if (obstype=="normal"){
   tmp.md <- apply(tmp.th, 2, median)
   tmp.l <- apply(tmp.th, 2, quantile, probs=plow)
   tmp.u <- apply(tmp.th, 2, quantile, probs=phigh)
  }
  if (obstype=="poisson"){
   tmp.md <- exp(apply(tmp.th, 2, median))
   tmp.l <- exp(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  }
  if (obstype=="binomial"){
   tmp.md <- plogis(apply(tmp.th, 2, median))
   tmp.l <- plogis(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- plogis(apply(tmp.th, 2, quantile, probs=phigh))
  }
   out <- list(postmed=tmp.md, bci.lower=tmp.l, bci.upper=tmp.u)
   out

}



extract_lambda <- function(mfit, alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.p <- rstan::extract(mfit, "LLambda")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.p1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.p1 <- rbind(tmp.p1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.p1)[[2]], pattern="LLambda")
	  tmp.p <- tmp.p1[ , ath]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="LLambda")
 	  tmp.p <- tmp.p1[ , ath]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="LLambda")
	  tmp.p <- as.matrix(mfit[ , ath])
  }
  plow <- alpha/2
  phigh <- 1 - alpha/2

   if (is.array(tmp.p)==TRUE){
     tmp.md <- apply(tmp.p, 2, median)
     tmp.l <- apply(tmp.p, 2, quantile, probs=plow)
     tmp.u <- apply(tmp.p, 2, quantile, probs=phigh)
   }
   if (is.array(tmp.p)==FALSE){
     tmp.md <- median(tmp.p)
     tmp.l <- quantile(tmp.p, probs=plow)
     tmp.u <- quantile(tmp.p, probs=phigh)   
   }
   out <- list(postmed=exp(tmp.md), bci.lower=exp(tmp.l), bci.upper=exp(tmp.u))
   out

}


plot_trace <- function(postob, vname, pscale="original", stack=FALSE, colset="color"){

   # Error checks
   if (missing(postob)) stop("Need to specify object containing posterior draws.")
   if (class(postob)!="array") stop("Posterior draws must be stored in an array.  See 'extract' or 'as.array' functions in rstan.")
   ## add check for vname specification - must be a character and be in dimnames of postob
   if (missing(vname)) stop("Need to specify variable name of parameter to plot.")
   if (!(vname %in% dimnames(postob)$par)) stop("Variable name must be a parameter name in postob.")
   if (!(pscale %in% c("original", "log", "inv"))) stop("Check specification of 'pscale'. Must be either 'original', 'log', or 'inv'.")
   if (!(colset %in% c("color", "gray", "black"))) stop("Check specification of 'colset'. Must be either 'color', 'gray', or 'black'.")

    # set up plots
	nc <- dim(postob)[2] #number of chains
	nit <- dim(postob)[1] #number of iterations per chain
	tnit <- nc*nit  #total iterations
	exv <- postob[, , vname]
	vrng <- range(exv)
	#chtck <-
	if (colset=="color") vcv <- rainbow(nc)
	if (colset=="black") vcv <- rep("black", nc)
	if (colset=="gray") {
		glev <- seq(0, 1 - 1/nc, by=1/nc)
		vcv <- gray(glev)
	}
  if (stack==FALSE){
  	if (pscale=="original") {
  		plot(1:tnit, 1:tnit, type="n", ylim=vrng,
  			main=paste("trace of",vname, "by chain") , xlab="iteration", ylab=vname)
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, exv[ ,ii], lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  	if (pscale=="log") {
  		plot(1:tnit, 1:tnit, type="n", ylim=log(vrng),
  			main=paste("trace of log",vname, "by chain") , xlab="iteration", ylab=paste("log",vname) )
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, log(exv[ ,ii]), lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  	if (pscale=="inv") {
  		plot(1:tnit, 1:tnit, type="n", ylim=rev(1/vrng),
  			main=paste("trace of inverse",vname, "by chain") , xlab="iteration", ylab=paste("1/",vname,sep=""))
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, 1/exv[ ,ii], lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  } # end stack false
	if (stack==TRUE){
		par(mfrow=c(4,1), oma=c(2, 2, 2, 0), mar=c(3,2,2,1) )
		if (pscale=="original") {
			for (ii in 1:nc){
			  plot(1:nit, exv[ ,ii], ylim=vrng, main="", xlab="", ylab="",
			  	type="l", lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=vname)
			mtext(side=3, outer=T, line=1, text=paste("trace of",vname, "by chain") )

		}
		if (pscale=="log") {
			for (ii in 1:nc){
			  plot(1:nit, log(exv[ ,ii]), ylim=log(vrng), main="", xlab="", ylab="",
			  	type="l",lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=paste("log",vname))
			mtext(side=3, outer=T, line=1, text=paste("trace of log",vname, "by chain") )
		}
		if (pscale=="inv") {
			for (ii in 1:nc){
			  plot(1:nit, 1/exv[ ,ii], ylim=rev(1/vrng), main="", xlab="", ylab="",
			  	type="l", lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=paste("1/",vname,sep=""))
			mtext(side=3, outer=T, line=1, text=paste("trace of inverse",vname, "by chain")  )
		}
  }  #end stack true
  par(mfrow=c(1,1), oma=c(0, 0, 0, 0), mar=c(5,4,4,2) )

}
