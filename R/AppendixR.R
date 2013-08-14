################################################################################
#  AppendixR.R
################################################################################
# Distorted-distance models for directional dispersal: 
# a general framework with application to a wind-dispersed tree
# Appendix R, Example R-script for maximum likelihood fitting of anisotropic dispersal kernels,  
# specifically elliptic distorted-distance kernels, through inverse modelling
# Bram van Putten, Marco D. Visser, Helene C. Muller-Landau and Patrick A. Jansen
#################################################################################  
# paramlist gives all the parameters
# paramlist=list(alpha,lambda,psi,beta,gamma,deltax,deltay,k) where 
# alpha  = seed production parameter (must be fit)
# lambda = rate parameter from the exponential distribution (must be fit)
# beta    = coherency parameter
# gamma = drift parameter (see table 1) 
# psi     = rotation angle parameter
# deltax,deltay = Shift parameters
# k      = negative binomial clumping parameter (must be fit)

# whichelliptic specifies which elliptic model to use:
# "full" for full elliptic model: need to specify all parameters: alpha, lambda, gamma, beta, psi, deltax,deltay,k
# "constrained-shift" for constrained-shift elliptic model: alpha, lambda, gamma, beta, deltax,deltay,k (psi is set equal to the shift direction)
# "no-shift" for no-shift elliptic model: parameters alpha, lambda, gamma, beta, psi, k (deltax and deltay are zero)
# "isotropic" for isotropic model: only parameters alpha, lambda, and k are relevant

# the following function simulates a dataset, fits parameters to the simulated dataset, and makes graphs
# to run this code, copy and paste the commands within this function into the R console after sourcing the entire code file
copyandpaste=function() {
    simwhichelliptic="no-shift"
    simparamlist=list(alpha=2,lambda=20,psi=-0.5,beta=1,gamma=1,k=1)
    testdata=simulatedata(whichelliptic=simwhichelliptic, paramlist=simparamlist, ntrap=200,ntree=200) 
    fitresults=fitmodel (whichelliptic="no-shift", startparamlist=list(alpha=2,lambda=50,psi=0,beta=1,gamma=1,k=1),fitmethod='Nelder-Mead',
            trapdata=testdata$trapdata,treedata=testdata$treedata)
    par(mfrow=c(1,2),las=1,mar=c(4,4,1,1),bty='o')
	graphseedshadow(whichelliptic=simwhichelliptic,paramlist=simparamlist,graphtitle="Simulated model")
	graphseedshadow(whichelliptic=fitresults$whichelliptic,paramlist=fitresults$paramlist,graphtitle="Fitted model")
} #end copyandpaste

####################### Define functions ###############################

# distance distribution Fx; here we only use the exponential distribution as an example (see Table 2)
Fx=function(r,lambda=100){dexp(r,1/lambda)}

# fecundity model phi; as in the main text.
phi=function(z,alpha=1000){exp(alpha)*z^2}

# Full eliptic distorted-distance function
rdist=function(x,y,Xp,Yp,beta=1,gamma,psi,deltax,deltay){
# x,y     = Cartesian coordinates of trap
# Xp, Yp  = Cartesian coordinates of tree
    co=cos(psi)
    si=sin(psi)
    xx=(x-Xp)+deltax
    yy=(y-Yp)+deltay
    u=((xx*co)+(yy*si))/beta
    v=-(xx*si)+(yy*co)
    return(dist1par(u,v,gamma))
}

# function for use in dist5par above
dist1par=function(u,v,gamma){
# gamma = Drift parameter (see table 1)
    a=-2*gamma*sqrt(gamma^2+1)*(u+gamma)*sqrt(u^2+2*gamma*u+v^2+gamma^2+1)
    b=u^2*(2*gamma^2+1)
    c=2*gamma*u*(2*gamma^2+1)
    d=(gamma^2+1)*(v^2+2*gamma^2)
    return(sqrt(a+b+c+d))
}

# Function to create  a distorted distance matrix between traps and trees
create.distmat=function(traplocations,treelocations,gamma,psi,beta=1,deltax=0,deltay=0){
DistMat=matrix(ncol=dim(treelocations)[1],nrow=dim(traplocations)[1])	
    for(i in 1:dim(treelocations)[1]){
        		DistMat[,i]=rdist(traplocations$x,traplocations$y,treelocations$x[i],
       			treelocations$y[i],gamma=gamma,psi=psi,beta=beta,deltax=deltax,deltay=deltay)	
    }
    return(DistMat)
}

# Seed shadow using polar coordinates
SS=function(r,dbh,alpha,lambda,quadsize=0.5,gamma=0,beta=1){
# r = distorted distance to source
# dbh = size of adult (diameter at breast height).
# alpha,lambda = parameters of the dispersal kernel (fx) and fecundity model (phi)
# quadsize = area of seedtrap in m2 [0.5 m2 in this study]
    return(quadsize*phi(z=dbh,alpha=alpha)*((1/(2*pi*beta*sqrt(gamma^2+1)))/r)*Fx(r=r,lambda=lambda))
}

#Function to calculate expected seed density (Mu)
Mu=function(traplocations,treelocations,paramlist){
# traplocations = dataframe with trap coordinates (x,y) and seed counts (n)
# treelocations = dataframe with tree coordinates (x,y) and size (dbh)
# paramlist        = full (translated) set of parameters, e.g. start parameters

	# build distorted distance matrix
	R=create.distmat(traplocations,treelocations,beta=paramlist$beta,gamma=paramlist$gamma,
	psi=paramlist$psi,deltax=paramlist$deltax,deltay=paramlist$deltay)
	
	#create expected seed density matrix P
	P=R
	for(i in 1:length(treelocations[,1])){
		P[,i]=SS(R[,i],dbh=treelocations$dbh[i],alpha=paramlist$alpha,lambda=paramlist$lambda,
			gamma=paramlist$gamma,beta=paramlist$beta)
	}
	
    #calculate expected seed densities per trap, and return values.
    return(apply(P,1,sum))
} # end Mu



###################### Fit models to data  ###################

fitmodel=function(whichelliptic="no-shift", startparamlist=c(alpha=2,lambda=0.02,psi=0,beta=1,gamma=1,k=1),
        fitmethod='Nelder-Mead',trapdata,treedata) {
# whichelliptic sets type of model to fit (see top of code file)
# startparamlist sets start parameters c(alpha,lambda,psi,beta,gamma,deltax,deltay,k) (see top of code file)
# fitmethod sets optimization algorithm: Default is 'Nelder-Mead'; use 'L-BFGS-B' to constrain search dimensions (see ?optim)
    startparamlist=checkparlist(whichelliptic,startparamlist)
    if (is.na(startparamlist$alpha))     return()
    startparvector=parvectorfromparlist(whichelliptic,startparamlist)
    if (is.na(startparvector[1])) return()

    # fit model to data, this may take a while depending on your system.
    estimation=optim(par=startparvector,fn=Likfunc,method=fitmethod,whichelliptic=whichelliptic,
            trapdata=trapdata,treedata=treedata)

	results=list(whichelliptic=whichelliptic,loglike=-estimation$value,
            AIC=2*estimation$value+2*length(estimation$par),
            paramlist=parlistfromparvector(whichelliptic,estimation$par))
	print(results)
	return(results)
}

# define the likelihood function - calculates the negative of the log likelihood
Likfunc=function(optimpar,whichelliptic,trapdata,treedata){
	paramlist=parlistfromparvector(whichelliptic,optimpar)
    ExpDen=Mu(trapdata,treedata,paramlist)
    sum(-dnbinom(trapdata$n,size=paramlist$k,mu=ExpDen,log=T))
}


parvectorfromparlist=function(whichelliptic,parlist) {
# Translates from a full named list of all model parameters
# to a vector of just the free parameters to be optimized over
    if (whichelliptic=="isotropic") {
        parvector=c(parlist$alpha,parlist$lambda,parlist$k)
        return(parvector)
    } else if (whichelliptic=="no-shift") {
        parvector=c(parlist$alpha,parlist$lambda,parlist$beta,parlist$gamma,
                parlist$psi,parlist$k)
        return(parvector)
    } else if (whichelliptic=="constrained-shift") {
        parvector=c(parlist$alpha,parlist$lambda,parlist$beta,parlist$gamma,
                parlist$deltax,parlist$deltay,parlist$k)
        return(parvector)
    } else if (whichelliptic=="full") {
        parvector=c(parlist$alpha,parlist$lambda,parlist$beta,parlist$gamma,
                parlist$psi,parlist$deltax,parlist$deltay,parlist$k)
        return(parvector)
    } else {
        print(paste("ERROR: whichelliptic=",whichelliptic,"not recognized!"))
        print("Must be one of the following: isotropic, no-shift, constrained-shift, full")
        return(NA)
    }
} 


parlistfromparvector=function(whichelliptic,parvector) {
# Translates from a vector list of the free parameters used in the optimization
# to a full named list of all model parameters used in the calls to Mu, etc.
    if (whichelliptic=="isotropic") {
        parlist=list(alpha=parvector[1],lambda=parvector[2],beta=1,gamma=0,
                psi=0,deltax=0,deltay=0,k=parvector[3])
        return(parlist)
    } else if (whichelliptic=="no-shift") {
        parlist=list(alpha=parvector[1],lambda=parvector[2],beta=parvector[3],gamma=parvector[4],
                psi=parvector[5],deltax=0,deltay=0,k=parvector[6])
        return(parlist)
    } else if (whichelliptic=="constrained-shift") {
        parlist=list(alpha=parvector[1],lambda=parvector[2],beta=parvector[3],gamma=parvector[4],
                psi=atan(parvector[6]/parvector[5]),deltax=parvector[5],deltay=parvector[6],k=parvector[7])
        return(parlist)
    } else if (whichelliptic=="full") {
        parlist=list(alpha=parvector[1],lambda=parvector[2],beta=parvector[3],gamma=parvector[4],
                psi=parvector[5],deltax=parvector[6],deltay=parvector[7],k=parvector[8])
        return(parlist)
    } else {
        print(paste("ERROR: whichelliptic=",whichelliptic,"not recognized!"))
        print("Must be one of the following: isotropic, no-shift, constrained-shift, full")
        return(list(alpha=NA))
    }
}

checkparlist=function(whichelliptic,parlist) {
# checks that parameters have the right values given the type of fit
# Translates from a full named list of all model parameters
# to a vector of just the free parameters to be optimized over
    print("Checking parameters specified in parlist.")
    print("If any required parameters are missing, an error will occur.")
    if (whichelliptic=="isotropic") {
        print("Since model is isotropic, setting beta=1, gamma=0, psi=0, deltax=0, deltay=0.")
        print("For the isotropic model, must specify alpha, lambda and k.")
        parlist=list(alpha=parlist$alpha,lambda=parlist$lambda,beta=1,gamma=0,
            psi=0,deltax=0,deltay=0,k=parlist$k)
    } else if (whichelliptic=="no-shift") {
        print("Since model is no-shift, setting deltax=0,deltay=0")
        print("For the no-shift model, must specify alpha, lambda, beta, gamma, psi, and k.")
        parlist=list(alpha=parlist$alpha,lambda=parlist$lambda,beta=parlist$beta,gamma=parlist$gamma,
            psi=parlist$psi,deltax=0,deltay=0,k=parlist$k)
    } else if (whichelliptic=="constrained-shift") {
        print("Since model is constrained-shift, setting psi=arctan(deltay/deltax)")
        print("For the constrained model, must specify alpha, lambda, beta, gamma, deltax, deltay, and k.")
        parlist=list(alpha=parlist$alpha,lambda=parlist$lambda,beta=parlist$beta,gamma=parlist$gamma,
            psi=parlist$deltay/parlist$deltax,deltax=parlist$deltax,deltay=parlist$deltay,k=parlist$k)
    } else if (whichelliptic=="full") {
        print("For the full model, must specify alpha, lambda, beta, gamma, psi deltax, deltay, and k.")
        parlist=list(alpha=parlist$alpha,lambda=parlist$lambda,beta=parlist$beta,gamma=parlist$gamma,
            psi=parlist$psi,deltax=parlist$deltax,deltay=parlist$deltay,k=parlist$k)
    } else {
        print(paste("ERROR: whichelliptic=",whichelliptic,"not recognized!"))
        print("Must be one of the following: isotropic, no-shift, constrained-shift, full")
        return(list(alpha=NA))
    }
    print("Parameter check complete.")
    return(parlist)
} 

######################### Build simulated datasets to test fitting ################################

# Seed trap creator: choose random trap locations
traplocator=function(t=1000,plotlimits=c(0,1000,0,500)){
# t traps are created in a square defined by limits
    x=runif(t,plotlimits[1],plotlimits[2])
    y=runif(t,plotlimits[3],plotlimits[4])
    return(data.frame(x,y))
}

# Tree creator: choose random tree locations
treelocator=function(a=1000,plotlimits=c(0,1000,0,500),dbhmin=300,dbhrate=0.005){
# a trees are created in a square defined by limits
# with exponentially distributed size (dbh in mm) defined by dbhrate
    x=runif(a,plotlimits[1],plotlimits[2])
    y=runif(a,plotlimits[3],plotlimits[4])
    dbh=dbhmin+rexp(10,dbhrate)
    return(data.frame(x,y,dbh))
}

simulatedata=function(whichelliptic="no-shift",paramlist=list(alpha=2,lambda=0.05,psi=6,beta=0.5,gamma=1,k=1),
        ntrap=200,ntree=200,plotlimits=c(0,1000,0,500),dbhmin=300,dbhrate=0.005) {
    paramlist=checkparlist(whichelliptic,paramlist)
    trapdata=traplocator(ntrap,plotlimits) # create ntrap seed traps
    treedata=treelocator(ntree,plotlimits,dbhmin,dbhrate) # create ntree trees

    # "sample seeds" in newly created fake forest
    expectedseeds=Mu(trapdata,treedata,paramlist)
    # realized trap counts
    trapdata$n=rnbinom(dim(trapdata)[1],size=paramlist$k,mu=expectedseeds)
    return(list(treedata=treedata,trapdata=trapdata))
}

######################## Graph a seed shadow #####################
# Graph a seed shadow (as in Figure 7 A)
graphseedshadow=function(whichelliptic,paramlist,xlim=c(-40,80),ylim=c(-80,40),dbh=500,graphtitle="",levels=c(1000,500,300,200,100,50,10)) {
    # plot details
    paramlist=checkparlist(whichelliptic,paramlist)
    x=seq(xlim[1],xlim[2],length=500)
    y=seq(ylim[1],ylim[2],length=500)
    
    # calculate seed density map for generating model, for a tree of size 500 mm dbh
    rdistmatrix=outer(x,y,rdist,Xp=0,Yp=0,psi=paramlist$psi,beta=paramlist$beta,
            gamma=paramlist$gamma,deltax=paramlist$deltax,deltay=paramlist$deltay)
    mumatrix=SS(rdistmatrix,dbh=dbh,alpha=paramlist$alpha,lambda=paramlist$lambda,
            gamma=paramlist$gamma,beta=paramlist$beta)
    contour(x,y,mumatrix,levels=c(1000,500,300,200,100,50,10),drawlabels=T, col='grey40',lwd=1.2,labcex=1,ylab='y',xlab='x',main=graphtitle)
    abline(v=0,h=0,lty='dashed')
    vec=(xlim[2]-xlim[1])/5 #vec = arrow length
    arrows(0,0,cos(paramlist$psi)*vec,sin(paramlist$psi)*vec,col='black',lwd=2,length=0.15)
} # end graphseedshadow

################################ End of code ###################################
