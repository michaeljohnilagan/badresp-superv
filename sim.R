# required packages and helper functions
library("furrr")
library("future")
library("mvtnorm")
source("helper-general.R")
source("helper-superv.R")

# run cell replicate
runsimrepl = function(seed,humandata,pointscale,nrilist,n1tr,respstyletr,
n,contam,n0tr,respstylete,classifier) {
	# set seed
	set.seed(seed)
	# sample sizes
	n0 = round(n*(1-contam))
	n1 = n-n0
	# generate test data likert
	zte0 = samplerows(n0,data=humandata,replace=TRUE) # human
	zte1 = rbot(n1,pointscale=pointscale,numitems=ncol(humandata),
	prob=respstylete) # bot
	yte = c(rep("human",nrow(zte0)),rep("bot",nrow(zte1))) # labels
	# generate training data likert
	ztr0 = samplerows(n0tr,data=humandata,replace=TRUE) # human
	ztr1 = rbot(n1tr,pointscale=pointscale,numitems=ncol(humandata),
	prob=respstyletr) # bot
	ytr = c(rep("human",nrow(ztr0)),rep("bot",nrow(ztr1))) # labels
	# convert to NRIs
	xtr = likert2nri(rbind(ztr0,ztr1),ref=ztr0,nrilist=nrilist,
	imputemissing=median(1:pointscale)) # train
	xte = likert2nri(rbind(zte0,zte1),ref=ztr0,nrilist=nrilist,
	imputemissing=median(1:pointscale)) # test
	# SCUMP
	fit = scump(ytr=ytr,xtr=xtr,xte=xte) # estimates
	scumppred = sniffbayes(xte=xte,fit=fit) # predictions
	# binary predictions
	if(classifier=="scumpbayes") {
		nominalreject = NA
		yhatp = scumppred
		yhat = ifelse(yhatp>0.5,"kill","spare") # SCUMP bayes
	} else if(classifier=="spec99") {
		nominalreject = 0.99
		yhatp = sniffnhst(xte=xte,ref=xtr[ytr=="human",,drop=FALSE])
		yhat = ifelse(yhatp>nominalreject,"kill",
		"spare") # specificity 99
	} else {
		stop("not implemented")
	}
	# evaluate classifier
	met = perfmet(y=yte,yhat=yhat,yhatp=yhatp)
	# put together
	return(list(zte0=zte0,zte1=zte1,yte=yte,ztr0=ztr0,ztr1=ztr1,ytr=ytr,
	xtr=xtr,xte=xte,fit=fit,yhat=yhat,yhatp=yhatp,met=met,
	scumppred=scumppred,nominalreject=nominalreject))
}

# run cell replicate, just the outcome measures
runsimreplom = function(seed,humandata,pointscale,nrilist,n1tr,respstyletr,
n,contam,n0tr,respstylete,classifier) {
	runsimrepl(seed,humandata,pointscale,nrilist,n1tr,respstyletr,
	n,contam,n0tr,respstylete,classifier)$met$outcomemeasures
}

# simulation study constants
numrepl = 1e3 # number of replicates per cell
humandata = with(new.env(),{
	datoriginal = read.csv("./data.csv",header=TRUE)[,1:32] # load data
	humandata = replace(datoriginal,datoriginal==-1,NA) # code missing
	humandata[apply(humandata,1,function(v) {
	!all(is.na(v))}),] # remove rows all missing
}) # human data
pointscale = 5 # number of likert categories
nrilist = list(mahal=function(x,ref) {
	mahal(x=x,ref=ref) # mahalanobis distance
},ptc=function(x,ref) {
	ptc(x=x,ref=ref) # person total correlation
},fmtrcoh=function(x,ref) {
	fmtrcoh(x=x,ref=ref,numpc=4,numiter=30) # FMT response coherence
},fmtrrel=function(x,ref) {
	fmtrrel(x=x,ref=ref,numpc=4,numiter=30) # FMT response reliability
})[c("mahal","ptc","fmtrcoh","fmtrrel")] # NRIs
n1tr = 2e3 # how many bots in training
respstyletr = c(1,1,1,1,1) # training bots response style

# simulation study factors
simfactors = list(n=c(500,1e3),contam=c(0.05,0.25,0.5,0.75,0.95),
n0tr=c(50,100,200,400),respstylete=setNames(list(c(1,1,1,1,1),c(1,2,4,2,1)),
c("unif","mrs")),classifier=c("scumpbayes","spec99"))

# try single replicate
with(new.env(),{
	foo = runsimrepl(2147,humandata=humandata,pointscale=5,
	nrilist=nrilist,n1tr=2e3,respstyletr=c(1,1,1,1,1),n=100,contam=0.5,
	n0tr=100,respstylete=c(1,1,1,1,1),classifier="scumpbayes") # execute
	foo$met # confusion table
})

# parallel processing
future::plan(future::multisession)

# simulation run
simresu = array(list(),sapply(simfactors,length)) # metrics replicates
simresumean = simresu # metrics mean
Sys.time(); for(i1 in 1:length(simfactors$n)) 
for(i2 in 1:length(simfactors$contam)) 
for(i3 in 1:length(simfactors$n0tr)) 
for(i4 in 1:length(simfactors$respstylete)) 
for(i5 in 1:length(simfactors$classifier)) {
	# scenario
	cellid = i1*100+i2*10+i3*1
	message(Sys.time())
	print(simfactors$n[i1])
	print(simfactors$contam[i2])
	print(simfactors$n0tr[i3])
	print(simfactors$respstylete[[i4]])
	print(simfactors$classifier[i5])
	# work
	curr = furrr::future_map(cellid+1:numrepl,runsimreplom,
	humandata=humandata,pointscale=pointscale,nrilist=nrilist,n1tr=n1tr,
	respstyletr=respstyletr,n=simfactors$n[i1],
	contam=simfactors$contam[i2],n0tr=simfactors$n0tr[i3],
	respstylete=simfactors$respstylete[[i4]],
	classifier=simfactors$classifier[i5],
	.progress=FALSE,.options=furrr::furrr_options(seed=NULL))
	# consolidate
	simresu[[i1,i2,i3,i4,i5]] = do.call(rbind,curr) # metrics
	simresumean[[i1,i2,i3,i4,i5]] = colMeans(simresu[[i1,i2,i3,i4,
	i5]]) # means
}; Sys.time()
warns = warnings() # save warnings
print(warns)

# table of results
simtabpar = do.call(expand.grid,{
	with(simfactors,list(n=n,contam=contam,n0tr=n0tr,
	respstyle=names(respstylete),classifier=classifier))
}) # scenarios
simtabom = with(new.env(),{
	simtabid = sapply(simtabpar,function(v) {
		as.integer(factor(v))
	}) # convert scenarios table to indices
	t(apply(simtabid,1,function(r) {
		simresumean[[r[1],r[2],r[3],r[4],r[5]]] 
	})) # loop over table rows
}) # outcome measures
simtab = cbind(simtabpar,simtabom) # full table

# visualize outcome measures
vis = function(i1,i4,i5,xlab,ylab,legendpos=NULL,main=NULL) {
	# assert
	stopifnot(i1%in%1:length(simfactors$n))
	stopifnot(i4%in%1:length(simfactors$respstyle))
	stopifnot(i5%in%1:length(simfactors$classifier))
	# coordinates
	coords = subset(simtab,n==simfactors$n[i1]&
	respstyle==names(simfactors$respstylete)[i4]&
	classifier==simfactors$classifier[i5])
	# plotting
	plot(simtab[[xlab]],simtab[[ylab]],type="n",xlab=xlab,ylab=ylab,
	main=main) # initialize
	sapply(1:length(simfactors$contam),function(i) {
		currline = subset(coords,contam==simfactors$contam[i])
		points(currline[[xlab]],currline[[ylab]],col=i+1,type="b",
		pch=c(2,rep(1,length(simfactors$n0tr)-2),3))
	}) # trend lines per contamination rate
	if(!is.null(legendpos)) {
		legend(legendpos,legend=paste(100*simfactors$contam,"%",sep=""),
		lwd=2,lty=1,col=1:length(simfactors$contam)+1,bg="white")
	} # legend
	return(coords)
}

# end session
save.image("badrespsuperv-v20220505t2140.RData")
devtools::session_info()
