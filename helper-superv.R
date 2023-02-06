# function: convert from likert to NRI
likert2nri = function(likert,ref,nrilist,imputemissing=NA) {
	# impute
	likert = as.matrix(likert)
	likert = replace(likert,is.na(likert),imputemissing)
	ref = as.matrix(ref)
	ref = replace(ref,is.na(ref),imputemissing)
	# convert to NRI
	nri = as.data.frame(sapply(nrilist,function(f) {
		unname(f(likert,ref))
	}))
	return(nri)
}

# function: SCUMP estimates
scump = function(ytr,xtr,xte) {
	# assert
	stopifnot(length(setdiff(ytr,0:1))==0|length(setdiff(ytr,
	c("bot","human")))==0)
	stopifnot(length(ytr)==nrow(xtr))
	stopifnot(ncol(xtr)==ncol(xte))
	# class label format
	if(!is.numeric(ytr)) {
		ytr = ifelse(ytr=="bot",1,0)
	}
	# estimate class feature distributions
	mu0 = colMeans(xtr[ytr==0,,drop=FALSE])
	mu1 = colMeans(xtr[ytr==1,,drop=FALSE])
	sig0 = cov(xtr[ytr==0,,drop=FALSE])
	sig1 = cov(xtr[ytr==1,,drop=FALSE])
	# estimate mixing proportions
	dens0 = mvtnorm::dmvnorm(xte,mu0,sig0)
	dens1 = mvtnorm::dmvnorm(xte,mu1,sig1)
	dll = function(lam) {
		sum((dens1-dens0)/(lam*dens1+(1-lam)*dens0))
	} # derivative of log likelihood
	if(sign(dll(0))==sign(dll(1))) {
		loglik = function(lam) {
			sum(log(lam*dens1+(1-lam)*dens0))
		}
		lam = ifelse(loglik(1)>loglik(0),1,0)
	} else {
		lam = uniroot(dll,0:1)$root
	} # bisection method
	# put together
	return(list(lam=lam,mu0=mu0,sig0=sig0,mu1=mu1,sig1=sig1))
}

# function: predicted probability of bot
sniffbayes = function(xte,fit) {
	dens0 = mvtnorm::dmvnorm(xte,fit$mu0,fit$sig0)
	dens1 = mvtnorm::dmvnorm(xte,fit$mu1,fit$sig1)
	nume = fit$lam*dens1
	deno = fit$lam*dens1+(1-fit$lam)*dens0
	return(nume/deno)
}

# function: prediction probability of null class
sniffnhst = function(xte,ref) {
	nulldist = mahalanobis(ref,colMeans(ref),cov(ref))
	obs = mahalanobis(xte,colMeans(ref),cov(ref))
	pval = sapply(obs,function(v) {
		mean(v>nulldist)
	})
	return(pval)
}
