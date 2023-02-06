# function: generate bots
rbot = function(n,pointscale,numitems,prob) {
	samp = t(replicate(n,{
		sample(1:pointscale,size=numitems,prob=prob,replace=TRUE)
	})) # each row is a participant
	return(samp)
}

# function: sample from data
samplerows = function(n,data,replace=TRUE) {
	sampindex = sample(1:nrow(data),size=n,replace=replace)
	return(as.matrix(data[sampindex,,drop=FALSE]))
}

# function: principal component analysis loadings and scores
pcaloadscore = function(x,numpc) {
	r = cor(x) # correlation matrix
	loading = with(eigen(r),vectors%*%
	sqrt(diag(values)))[,1:numpc,drop=FALSE] # loadings, signs not adjusted
	signtot = ifelse(colSums(loading)==0,+1,sign(colSums(loading)))
	loadingsigned = loading%*%diag(signtot,nrow=numpc,
	ncol=numpc) # loadings, signs adjusted
	score = scale(x)%*%solve(r,loadingsigned) # scores
	return(list(loadings=loadingsigned,scores=score))
}

# function: FMT matrix of item characteristics
fmtmic = function(x,numpc,numiter=30) {
	# function to normalize
	eucnormalize = function(v) {
		v/sqrt(sum(v^2))
	}
	# initial representation
	itemchar = pcaloadscore(x,numpc=numpc)$loadings
	# iterative PCA
	for(iter in 1:numiter) {
		itemchar = pcaloadscore(itemchar,
		numpc=numpc)$scores # orthogonalize
		itemchar = t(apply(itemchar,1,eucnormalize)) # normalize
	}
	return(itemchar)
}

# function: response coherence
fmtrcoh = function(x,ref,numpc,numiter=30,undefined=0) {
	# FMT matrix of item characteristics
	itemchar = fmtmic(ref,numpc=numpc,
	numiter=numiter) # rows are items; columns are factors
	# strategy and length
	respstrat = t(cor(itemchar,t(as.matrix(x)))) # response strategy
	coher = sqrt(rowSums(respstrat^2)) # euclidean length
	narm = ifelse(is.na(coher),undefined,coher) # impute over undefined
	return(narm)
}

# function: get item pairs
getitempairs = function(m) {
	# assert
	stopifnot(isSymmetric(m))
	# initialize pairs
	paired = data.frame(a=numeric(0),b=numeric(0)) # empty data frame
	# loop appending pairs
	m[lower.tri(m,diag=TRUE)] = NA # upper triangle only
	while(TRUE) {
		newpair = which(m==max(m,na.rm=TRUE),arr.ind=TRUE)[1,]
		paired = rbind(paired,newpair)
		m[newpair,] = NA
		m[,newpair] = NA
		if(all(is.na(m))) {
			break
		}
	}
	return(setNames(paired,c("a","b")))
}

# function: response reliability
fmtrrel = function(x,ref,numpc,numiter=30,undefined=-1) {
	# function to normalize
	eucnormalize = function(v) {
		v/sqrt(sum(v^2))
	}
	# FMT matrix of item characteristics
	itemchar = fmtmic(ref,numpc=numpc,
	numiter=numiter) # rows are items; columns are factors
	# split half strategies
	paired = getitempairs(itemchar%*%t(itemchar)) # get item pairs
	itemchar_a = itemchar[paired$a,,drop=FALSE]
	itemchar_b = itemchar[paired$b,,drop=FALSE]
	x_a = as.matrix(x[,paired$a,drop=FALSE])
	x_b = as.matrix(x[,paired$b,drop=FALSE])
	respstrat_a = t(cor(itemchar_a,t(x_a))) # strategy A
	respstrat_b = t(cor(itemchar_b,t(x_b))) # strategy B
	# inner product and correction
	rel = sapply(1:nrow(x),function(i) {
		normeda = eucnormalize(respstrat_a[i,])
		normedb = eucnormalize(respstrat_b[i,])
		sum(normeda*normedb)
	}) # inner product
	relsb = sapply(rel,function(r) {
		max(c(2*r/(1+r),-1))
	}) # spearman brown corrected
	narm = ifelse(is.na(relsb),undefined,relsb) # impute over undefined
	return(narm)
}

# function: mahalanobis distance
mahal = function(x,ref) {
	unname(sqrt(mahalanobis(x,center=colMeans(ref),cov=cov(ref))))
}

# function: person total correlation
ptc = function(x,ref,undefined=-1) {
	thecor = apply(x,1,function(v) {
		cor(v,colMeans(ref))
	}) # compute correlation
	narm = ifelse(is.na(thecor),undefined,thecor) # impute over undefined
	return(narm)
}

# function: plot NRI space
nriplot = function(xte,xtr=NULL,yhat=NULL,yte=NULL,ytr=NULL,plottrain=TRUE,
getpc=FALSE,trplotlim=FALSE) {
	# principal component scores
	if(getpc) {
		if(!is.null(xtr)) {
			stdtr = scale(xtr)
			stdte = scale(xte,center=colMeans(xtr),
			scale=sapply(xtr,sd))
			eig = eigen(cor(stdtr))
			loctr = stdtr%*%eig$vectors[,1:2]
			locte = stdte%*%eig$vectors[,1:2]
		} else {
			eig = eigen(cor(xte))
			locte = scale(xte)%*%eig$vectors[,1:2]
		}
		thexlab = "PC1"
		theylab = "PC2"
	} else {
		if(!is.null(xtr)) {
			loctr = xtr
		}
		locte = xte
		thexlab = colnames(xte)[1]
		theylab = colnames(xte)[2]
	} # convert to PC scores or not
	# draw plot
	if(trplotlim) {
		thexlim = range(loctr[,1])
		theylim = range(loctr[,2])
	} else {
		thexlim = range(locte[,1])
		theylim = range(locte[,2])
	}
	plot(locte[,1],locte[,2],type="n",xlab=thexlab,ylab=theylab,
	xlim=thexlim,ylim=theylim)
	if(plottrain) {
		if(!is.null(ytr)) {
			points(loctr[,1],loctr[,2],pch=ifelse(ytr=="bot",4,1),
			col="green")
		} else {
			points(loctr[,1],loctr[,2],col="green",pch=15)
		} # training labels available or not
	} # plot training set
	if(!is.null(yte)&!is.null(yhat)) {
		points(locte[,1],locte[,2],pch=ifelse(yte=="bot",4,1),
		col=ifelse(yhat=="kill","red","blue"),lwd=2)
	} else if(is.null(yte)&!is.null(yhat)) {
		points(locte[,1],locte[,2],col=ifelse(yhat=="kill","red",
		"blue"),pch=15)
	} else if(!is.null(yte)&is.null(yhat)) {
		points(locte[,1],locte[,2],pch=ifelse(yte=="bot",4,1),lwd=2)
	} else if(is.null(yte)&is.null(yhat)) {
		points(locte[,1],locte[,2],pch=15)
	}
	# return test set coordinates
	return(locte[,1:2])
}

# function: ROC curve and AUC
rocauc = function(yb,yc,add=FALSE,plot=TRUE,...) {
	# assert
	stopifnot(length(yb)==length(yc))
	# unique values
	sortuni = sort(unique(yc))
	# get ROC coordinates
	fp = sapply(sortuni,function(v) {
		mean(yc[yb==0]>v)
	}) # false positive rate
	tp = sapply(sortuni,function(v) {
		mean(yc[yb==1]>v)
	}) # true positive rate
	coords = data.frame(x=fp,y=tp)
	# draw plot
	if(plot) {
		if(!add) {
			plot(NA,NA,type="n",xlab="false positive",
			ylab="true positive",xlim=0:1,ylim=0:1)
		}
		with(coords,points(x,y,type="l",...))
	}
	# calculate AUC with trapezoids
	return(with(coords,-sum(diff(x)*(y[-1]+y[-nrow(coords)])/2)))
}

# function: performance metrics
perfmet = function(y,yhat,yhatp=NULL) {
	# assert
	stopifnot(length(setdiff(y,0:1))==0|length(setdiff(y,
	c("bot","human")))==0)
	stopifnot(length(setdiff(yhat,0:1))==0|length(setdiff(yhat,
	c("kill","spare")))==0)
	# confusion table
	confusion = table(y,yhat)
	# class label format
	if(!is.numeric(y)) {
		y = ifelse(y=="bot",1,0)
	}
	if(!is.numeric(yhat)) {
		yhat = ifelse(yhat=="kill",1,0)
	}
	# compute metrics
	acc = mean(y==yhat) # accuracy
	spec = mean(yhat[y==0]==0) # specificity
	sens = mean(yhat[y==1]==1) # sensitivity
	ppv = mean(y[yhat==1]==1) # positive predictive value
	npv = mean(y[yhat==0]==0) # negative predictive value
	killrate = mean(yhat==1) # kill rate
	if(!is.null(yhatp)) {
		auc = rocauc(yb=y,yc=yhatp,add=FALSE,plot=FALSE)
	} else {
		auc = NA
	}
	# put together
	outcomemeasures = list(acc=acc,spec=spec,sens=sens,ppv=ppv,npv=npv,
	killrate=killrate,auc=auc)
	return(list(confusion=confusion,
	outcomemeasures=unlist(outcomemeasures)))
}
