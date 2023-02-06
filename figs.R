# load workspace
load("badrespsuperv-v20220505t2140.RData")

# figure 1
foo = runsimrepl(650,humandata=humandata,pointscale=5,
nrilist=nrilist,n1tr=2e3,respstyletr=c(1,1,1,1,1),n=50,contam=0.5,n0tr=100,
respstylete=c(1,1,1,1,1),classifier="scumpbayes") # execute
pdf("intro.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
fooloc = with(foo,{
	nriplot(xte=setNames(xte[,c(1,2)],c("mahalanobis distance",
	"person total correlation")),
	xtr=NULL,yhat=NULL,yte=yte,ytr=ytr,
	plottrain=FALSE,getpc=FALSE,trplotlim=FALSE)
}) # make plot
with(new.env(),{
	y = ifelse(foo$yte=="bot",1,0)
	x = foo$xte[,c("mahal","ptc")]
	mod = glm(y~.,data=x,family=binomial(link="logit"))
	modcoef = coef(mod)
	abline(-modcoef[1:2]/modcoef[3],lty=1)
	legend("topright",c("human","bot"),pch=c(1,4),cex=1.35,bg="white")
}) # linear boundary by logistic regression
dev.off()

# figure 2a
pdf("hist50.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
curve(0.5*dnorm(x,2.5,1),from=-3,to=+7,lwd=2,lty=2,xlab="nonresponsivity",
ylab="density",main=NA,ylim=c(0,dnorm(0)))
curve(0.5*dnorm(x,0,1),from=-3,to=+7,lwd=2,add=TRUE)
abline(v=qnorm(0.95,0,1))
legend("topright",c("human","bot"),lwd=2,lty=1:2,cex=1.35,bg="white")
dev.off()

# figure 2b
pdf("hist90.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
curve(0.9*dnorm(x,2.5,1),from=-3,to=+7,lwd=2,lty=2,xlab="nonresponsivity",
ylab="density",main=NA,ylim=c(0,dnorm(0)))
curve(0.1*dnorm(x,0,1),from=-3,to=+7,lwd=2,add=TRUE)
abline(v=qnorm(0.95,0,1))
legend("topleft",c("human","bot"),lwd=2,lty=1:2,cex=1.35,bg="white")
dev.off()

# figure 3
pdf("rocnorm.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
with(new.env(),{
	mu0 = 0
	sig0 = 1
	mu1 = 2.5
	sig1 = 1
	thegrid = seq(+7,-3,length.out=200)
	coords = data.frame(fp=1-pnorm(thegrid,mu0,sig0),
	tp=1-pnorm(thegrid,mu1,sig1))
	plot(NA,NA,type="l",xlim=0:1,ylim=0:1,xlab="false positive rate",
	ylab="true positive rate",lwd=2)
	abline(0:1,lwd=1,col="gray")
	with(coords,points(fp,tp,type="l",lwd=2))
	cutoffspec = qnorm(0.95,mu0,sig0)
	abline(v=0.05,lty=2)
	cutoff50 = uniroot(function(x) {
		0.5*dnorm(x,mu1,sig1)-0.5*dnorm(x,mu0,sig0)
	},lower=-3,upper=+7)$root
	points(1-pnorm(cutoff50,mu0,sig0),1-pnorm(cutoff50,mu1,sig1))
	text(1-pnorm(cutoff50,mu0,sig0),1-pnorm(cutoff50,mu1,sig1)+0.07,"50%",
	cex=1.2)
	cutoff90 = uniroot(function(x) {
		0.9*dnorm(x,mu1,sig1)-0.1*dnorm(x,mu0,sig0)
	},lower=-3,upper=+7)$root
	points(1-pnorm(cutoff90,mu0,sig0),1-pnorm(cutoff90,mu1,sig1))
	text(1-pnorm(cutoff90,mu0,sig0),1-pnorm(cutoff90,mu1,sig1)-0.07,"90%",
	cex=1.2)
})
dev.off()
integrate(function(x) {
	1-pnorm(qnorm(1-x,0,1),2.5,1)
},lower=0,upper=1)$value # AUC
set.seed(224)
mean(replicate(10e3,{
	rnorm(1,0,1)<rnorm(1,2.5,1)
})) # also AUC

# figure 4a
foo1 = runsimrepl(650,humandata=humandata,pointscale=5,
nrilist=nrilist,n1tr=2e3,respstyletr=c(1,1,1,1,1),n=100,contam=0.5,n0tr=100,
respstylete=c(1,1,1,1,1),classifier="scumpbayes") # execute
pdf("demoscatscump.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
with(foo1,{
	nriplot(xte=setNames(xte[,1:2],c("mahalanobis distance",
	"person total correlation")),
	xtr=xtr[,1:2],yhat=yhat,yte=yte,ytr=ytr,
	plottrain=TRUE,getpc=FALSE,trplotlim=TRUE)
	legend("topright",c("flagged","spared","training"),pch=15,
	col=c("red","blue","green"),cex=1.35,bg="white")
	legend("bottomleft",c("human","bot"),pch=c(1,4),cex=1.35,bg="white")
}) # make plot
dev.off()
foo1$met$confusion

# figure 4b
foo2 = runsimrepl(650,humandata=humandata,pointscale=5,
nrilist=nrilist,n1tr=2e3,respstyletr=c(1,1,1,1,1),n=100,contam=0.5,n0tr=100,
respstylete=c(1,1,1,1,1),classifier="spec99") # execute
pdf("demoscatspec.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
with(foo2,{
	nriplot(xte=setNames(xte[,1:2],c("mahalanobis distance",
	"person total correlation")),
	xtr=xtr[,1:2],yhat=yhat,yte=yte,ytr=ytr,
	plottrain=TRUE,getpc=FALSE,trplotlim=TRUE)
	legend("topright",c("flagged","spared","training"),pch=15,
	col=c("red","blue","green"),cex=1.35,bg="white")
	legend("bottomleft",c("human","bot"),pch=c(1,4),cex=1.35,bg="white")
}) # make plot
dev.off()
foo2$met$confusion

boxplot(spec~classifier,data=subset(simtab,n==1e3&n0tr==100&respstyle!="unif"))
abline(h=0.99,lty=2)

# rename columns in simtab
colnames(simtab)[6:8] = c("accuracy","specificity","sensitivity")

# figure 6b
pdf("specaccunifscump.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
vis(2,1,1,"specificity","accuracy",main=NA)
abline(v=c(0.95,0.99),lty=3)
legend(x=0.5,y=0.52,c("5%","25%","50%","75%","95%"),lty=1,col=1:5+1,
title="contamination",cex=0.85)
legend(x=0.69,y=0.49,c("50","100","200","400"),pch=c(2,1,1,3),
title="training humans",cex=0.85)
dev.off()

# figure 6a
pdf("specaccunifnonscump.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
vis(2,1,2,"specificity","accuracy",main=NA)
abline(v=c(0.95,0.99),lty=3)
legend(x=0.5,y=0.52,c("5%","25%","50%","75%","95%"),lty=1,col=1:5+1,
title="contamination",cex=0.85)
legend(x=0.69,y=0.49,c("50","100","200","400"),pch=c(2,1,1,3),
title="training humans",cex=0.85)
dev.off()

# figure 6d
pdf("specaccmrsscump.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
vis(2,2,1,"specificity","accuracy",main=NA)
abline(v=c(0.95,0.99),lty=3)
legend(x=0.5,y=0.52,c("5%","25%","50%","75%","95%"),lty=1,col=1:5+1,
title="contamination",cex=0.85)
legend(x=0.69,y=0.49,c("50","100","200","400"),pch=c(2,1,1,3),
title="training humans",cex=0.85)

dev.off()

# figure 6c
pdf("specaccmrsnonscump.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
vis(2,2,2,"specificity","accuracy",main=NA)
abline(v=c(0.95,0.99),lty=3)
legend(x=0.5,y=0.52,c("5%","25%","50%","75%","95%"),lty=1,col=1:5+1,
title="contamination",cex=0.85)
legend(x=0.69,y=0.49,c("50","100","200","400"),pch=c(2,1,1,3),
title="training humans",cex=0.85)
dev.off()

# figure 5
pdf("demoroc.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
with(new.env(),{
	plot(NA,NA,xlim=0:1,ylim=0:1,xlab="false positive rate",
	ylab="true positive rate")
	abline(0:1,col="gray")
	col1 = "purple"
	col2 = "green4"
	lty1 = 1
	lty2 = 1
	auc1 = rocauc(ifelse(foo1$yte=="bot",1,0),foo1$scumppred,add=TRUE,
	plot=TRUE,lty=lty1,lwd=2,col=col1)
	auc2 = rocauc(ifelse(foo2$yte=="bot",1,0),foo2$yhatp,add=TRUE,
	plot=TRUE,lty=lty2,lwd=2,col=col2)
	legend("bottomright",c("SCUMP Bayes","99% specificity"),
	lty=c(lty1,lty2),lwd=2,col=c(col1,col2),cex=1.35,bg="white")
	with(as.list(foo1$met$outcomemeasures),points(1-spec,sens,lwd=2,
	pch=15,col=col1))
	with(as.list(foo2$met$outcomemeasures),points(1-spec,sens,lwd=2,
	pch=15,col=col2))
	message("SCUMP bayes AUC: ",auc1)
	message("99% nominal specificity AUC: ",auc2)
})
dev.off()

# figure 8
pdf("aucagg.pdf")
par(cex.lab=1.35, cex.axis=1.35, cex=1.35)
with(new.env(),{
	ref = subset(simtab,classifier%in%c("scumpbayes","spec99"))
	ref$group = paste(ref$classifier,ref$respstyle)
	m = aggregate(ref$auc,ref[c("group","n0tr")],mean)
	dotchart(m$x,labels=m$n0tr,groups=factor(m$group,
	levels=c("scumpbayes unif","scumpbayes mrs","spec99 unif",
	"spec99 mrs")),
	xlab="AUC")
})
dev.off()