#' Simulation of the Likelihood Ratio Statistic

library(MASS)
eLRT <- function(

	Z = NULL, # not include the intercept term// 
	E = E, # GReX
	G = G, #  cis-SNP matrix for GWAS data
	muxi = NULL,
	nsim = 10000,
	log_grid_hi = 8,
	log_grid_lo = -10,
	gridlength  = 100,#length for grid points of LRT//
	lambda0 = 0,
	seed = NA,
	parallel = c("no", "multicore", "snow"),#multicore can substantially improve the computation
	ncpus = 1L,
	cl = NULL) {

	REML = FALSE
	parallel <- match.arg(parallel)
	have_mc <- have_snow <- FALSE
	if (parallel != "no" && ncpus > 1L) {
		if (parallel == "multicore")
		have_mc <- .Platform$OS.type != "windows"
		else if (parallel == "snow")
		have_snow <- TRUE
	if (!have_mc && !have_snow)
	ncpus <- 1L
	}

	#print(have_mc)
	q <- 1
	m <- NCOL(G)     # no. of cis-SNPs snp
	if (is.null(Z)) {X <- as.matrix(cbind(1, E))}
	else            {X <- as.matrix(cbind(1, Z, E))}

	#X <- as.matrix(cbind(1, Z, E))
	n <- k <- NROW(X)     # no. of obs
	#k <- NROW(K)     # no. of random effects
	p <- NCOL(X)     # no of fixed effects

	#compute eigenvalues
	#mu <- (svd(sqrtSigma %*% t(qr.resid(qr(X), Z)), nu = 0, nv = 0)$d)^2
	#mu <- (svd(sqrtSigma %*% t(qr.resid(qr(X), Z)), nu = 0, nv = 0)$d)^2
	#xi <- (svd(sqrtSigma %*% t(Z), nu = 0, nv = 0)$d)^2
	if (is.null(muxi))
	{
	#print("compute mu and xi")
	if (m >= n) {
		K <- G %*% t(G)
		eig <- eigenK(K, sqrtk = TRUE)
		rm(K)
		xi  <- sort(eig$value, decreasing = TRUE)
		Px <- diag(n) - X %*% ginv(t(X) %*% X) %*% t(X)
		Px <- eig$sqrtK %*% Px %*% eig$sqrtK
		mu <- sort(eigenK(Px, sqrtk = FALSE)$value, decreasing = TRUE)
		rm(Px)
	}

	else 
	{
		K <- t(G) %*% G
		xi <- sort(eigenK(K, sqrtk = FALSE)$value, decreasing = TRUE)
		#Px <- diag(n) - X %*% ginv(t(X) %*% X) %*% t(X)
		#Px <- t(G) %*% Px %*% G

		GtX <- t(G) %*% X
		XtX <- t(X) %*% X
		Px <- K - GtX %*% ginv(XtX) %*% t(GtX)

		mu <- sort(eigenK(Px, sqrtk = FALSE)$value, decreasing = TRUE)
		#index        <- which(mu<1e-8)
		#if (length(index)>0) {mu[index] = 0}
		rm(K)
		rm(Px)
		rm(GtX)
		rm(XtX)
	}
	}

	else 
	{
	#print("not compute mu and xi")
	mu <- muxi[,1]
	xi <- muxi[,2]
	}

	#norm eigenvalues
	mu <- mu / max(mu,xi)
	xi <- xi / max(mu,xi)
	k  <- min(k, m)
	#print(k)
	lambda_grid <-c(0, exp(seq(log_grid_lo, log_grid_hi, length = gridlength - 1)))

	if (!is.na(seed))
	set.seed(seed)

	res <- if (ncpus > 1L && (have_mc || have_snow)) {
	nsim <- as.integer(ceiling(nsim/ncpus))
	if (have_mc) {#parallel = "multicore"
		#print("simulation is fitting")
		tmp <- parallel::mclapply(seq_len(ncpus), function(i){
		RLRsimCpp(p = as.integer(p), k = as.integer(k),
		n = as.integer(n), nsim = as.integer(nsim),
		g = as.integer(gridlength),
		q = as.integer(q), mu = as.double(mu),
		lambda = as.double(lambda_grid),
		lambda0 = as.double(lambda0), xi = as.double(xi),
		REML = as.logical(FALSE))
		}, mc.cores = ncpus)
	do.call(mapply, c(tmp, FUN=c))
	}
	else {
		if (have_snow) {
		if (is.null(cl)) {
		cl <- parallel::makePSOCKcluster(rep("localhost",ncpus))
		if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
		{parallel::clusterSetRNGStream(cl)}
		tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
		print("simulation is fitting")
		RLRsimCpp(p = as.integer(p), k = as.integer(k),
		n = as.integer(n), nsim = as.integer(nsim),
		g = as.integer(gridlength),
		q = as.integer(q), mu = as.double(mu),
		lambda = as.double(lambda_grid),
		lambda0 = as.double(lambda0), xi = as.double(xi),
		REML = as.logical(FALSE))
		})
		parallel::stopCluster(cl)
		do.call(mapply, c(tmp, FUN=c))
		}
	else {
		tmp <- parallel::parLapply(cl, seq_len(ncpus), function(i){
		print("simulation is fitting")
		RLRsimCpp(p = as.integer(p), k = as.integer(k),
		n = as.integer(n), nsim = as.integer(nsim),
		g = as.integer(gridlength),
		q = as.integer(q), mu = as.double(mu),
		lambda = as.double(lambda_grid),
		lambda0 = as.double(lambda0), xi = as.double(xi),
		REML = as.logical(FALSE))
		})
		do.call(mapply, c(tmp, FUN=c))
	}
	}
	}
	}
	else {#parallel = "no"
		print("simulation is fitting")
		RLRsimCpp(p = as.integer(p), k = as.integer(k),
		n = as.integer(n), nsim = as.integer(nsim),
		g = as.integer(gridlength),
		q = as.integer(q), mu = as.double(mu),
		lambda = as.double(lambda_grid),
		lambda0 = as.double(lambda0), xi = as.double(xi),
		REML = as.logical(FALSE))
	}
	#print("simulation finished")
	lambda <- lambda_grid[res$lambdaind]
	ret <- res$res
	attr(ret, "lambda") <- lambda_grid[res$lambdaind]
	return(ret)
}

aLRT = function(simLike, nmom) {
pchibarsqAUD=function(p,df=1,mix=0.5,a=1,lower.tail=TRUE,log.p=FALSE) 
{
  p=p/a
  df=rep(df,length.out=length(p))
  mix=rep(mix,length.out=length(p))
  c1=ifelse(df==1, if (lower.tail) 1
    else 0, pchisq(p,df-1,lower.tail=lower.tail))
  c2=pchisq(p,df,lower.tail=lower.tail)
  r=mix*c1+(1-mix)*c2
  if (log.p) 
    log(r)
  else r
}
pmix = matrix(NA,length(nmom),4)
for (j in 1:length(nmom))
{
  x=simLike[1:nmom[j]]
  nzero=which(x==0)
  if (length(nzero)>0)
   {
   pQ1=mean(x==0)
   x1=x[-nzero]
   k=var(x1)/(2*mean(x1))
   v=2*(mean(x1)^2)/var(x1)
   p_pQ=1-pchibarsqAUD(obsLike,df=v,mix=pQ1,a=k)
  }
  else 
  {
  pQ1=mean(x==0)
  k=var(x)/(2*mean(x))
  v=2*(mean(x)^2)/var(x)
  p_pQ=pchisq(obsLike/k,v,lower.tail = FALSE)
  }
  pmix[j,] = c(pvalue=p_pQ,k=k,v=v,p=pQ1)
}
 colnames(pmix)=c("pvalue","k","v","p")
 return(pmix)
}
