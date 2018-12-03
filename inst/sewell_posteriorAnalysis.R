# rm(list = setdiff(ls(),"LSMDN_result"))
attach(LSMDN_result)

#----------------------------------------------------
jpeg(file="Bin.jpg")
par(mar=c(4, 6, 2, 2) + 0.1)
plot(Bin[1:N],type="l",ylab=expression(beta[IN]),
     xlab="",cex.lab=2.5,cex.axis=1.5)
dev.off()
jpeg(file="Bout.jpg")
par(mar=c(4, 6, 2, 2) + 0.1)
plot(Bout[1:N],type="l",ylab=expression(beta[OUT]),
     xlab="",cex.lab=2.5,cex.axis=1.5)
dev.off()
jpeg(file="s2.jpg")
par(mar=c(4, 6, 2, 2) + 0.1)
plot(s2[1:N],type="l",ylab=expression(sigma^2),
     xlab="",cex.lab=2.5,cex.axis=1.5)
dev.off()
jpeg(file="t2.jpg")
par(mar=c(4, 6, 2, 2) + 0.1)
plot(t2[1:N],type="l",ylab=expression(tau^2),
     xlab="",cex.lab=2.5,cex.axis=1.5)
dev.off()

N1 <- N-burnin


#----------------------------------------------------
EBin <- mean(Bin[-c(1:burnin)])
EBout <- mean(Bout[-c(1:burnin)])
EBin;EBout
print(sum(Bin[-c(1:burnin)]>Bout[-c(1:burnin)])/N1)
Ew <- apply(w[,-c(1:burnin)],1,mean)
EX <- array(0,dim=c(p,TT,n))
for(it in c(1:N1)){
  EX <- EX + X[[burnin+it]]
}
EX <- EX/N1


xl <- c(1.15*min(EX[1,,]),1.15*max(EX[1,,]))
yl <- c(1.15*min(EX[2,,]),1.15*max(EX[2,,]))
lims <- range(c(xl,yl))

jpeg(file="latentPositions.jpg")
par(mar=c(2,2,4,2))
for(tt in 1:TT){
  plot(t(EX[,tt,]),xlim=lims,ylim=lims,xlab="",ylab="",
       pch=16,cex=0.75,xaxt="n",yaxt="n",
       main="Posterior mean latent positions showing temporal dynamics")
  if(tt>1) arrows(EX[1,tt-1,],EX[2,tt-1,],EX[1,tt,],EX[2,tt,],length=0.1)
#  if(tt==1) textxy(EX[1,tt,],EX[2,tt,],labs=c(1:26)[-21]) #load "calibrate" package
  if(tt<TT) par(new=TRUE)
}
dev.off()

#----------------------------------------------------
require("ROCR") # install.packages("ROCR")
predY <- array(0,dim=dim(Y))
for(tt in 1:TT){
  for(i in 1:n){
    for(j in c(1:n)[-i]){
      dijt <- sqrt(t(EX[,tt,i]-EX[,tt,j])%*%(EX[,tt,i]-EX[,tt,j]))
      expon <- EBin*(1-dijt/Ew[j])+EBout*(1-dijt/Ew[i])
      predY[i,j,tt] <- 1/(1+exp(-expon))
    }
  }
}

if(MissData){
  for(tt in 1:TT){
    ROCR_pred <- c(ROCR_pred,upperTriangle(predY[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_pred <- c(ROCR_pred,lowerTriangle(predY[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_Y <- c(ROCR_Y,upperTriangle(Y[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_Y <- c(ROCR_Y,lowerTriangle(Y[-Missing[[tt]],-Missing[[tt]],tt]))
  }
  pred <- prediction( predictions=ROC_pred,labels=ROC_Y )
  performance( pred, "auc")
  AUCPMean <- slot(performance( pred, "auc"),"y.values")[[1]]
}else{
  ROC_pred <- ROC_Y <- numeric(0)
  for(tt in 1:TT){
    ROC_pred <- c(ROC_pred,c(predY[,,tt][upper.tri(predY[,,tt])]))
    ROC_pred <- c(ROC_pred,c(predY[,,tt][lower.tri(predY[,,tt])]))
    ROC_Y <- c(ROC_Y,c(Y[,,tt][upper.tri(Y[,,tt])]))
    ROC_Y <- c(ROC_Y,c(Y[,,tt][lower.tri(Y[,,tt])]))
  }
  pred <- prediction( predictions=round(ROC_pred),labels=as.numeric(ROC_Y>0) )
  performance( pred, "auc")
  AUCPMean <- slot(performance( pred, "auc"),"y.values")[[1]]  
}
print(AUCPMean)

#----------------------------------------------------
require("animation") # install.packages("animation")
L <- seq(from=1,to=0,length=51)[-51]
oopt = ani.options(interval = 0.01, nmax = TT*length(L))
Xtemp2 <- t(EX[,1,])
for(tt in 2:TT)Xtemp2 <- rbind(Xtemp2,t(EX[,tt,]))
Xtemp2 <- Xtemp2 -
  kronecker(rep(1,n*TT),matrix(apply(Xtemp2,2,mean),1,p))

## use a loop to create images one by one
saveHTML(
{
  par(mar=c(2,2,2,2)+0.1)  
  for (tt in 1:(TT-1)) {
    for(l in 1:length(L)){
      plot(L[l]*t(EX[,tt,])+(1-L[l])*t(EX[,tt+1,]),xlab="",ylab="",xlim=lims,ylim=lims,
           pch=16,xaxt="n",yaxt="n", cex=2)
      ani.pause() ## pause for a while ('interval')
    }
  }
  
},img.name="LSMDNPlot",htmlfile="LSMDN.html",autoplay=FALSE,
interval=0.01,title="Latent Space Movements",
imgdir="LSMDNimages",ani.height=800,ani.width=800,
outdir="~/LSMDNMovie",
description="Posterior means of the latent space positions, with the movements interpolated")


#----------------------------------------------------
###
###Nodal Influence
###
require("circular") # install.packages("circular")
lambda <- 0.05
p0 <- 0.5


kappaij <- matrix(0,n,n)
temp <- numeric(TT-1)
X1 <- X2 <- array(0,c(n,TT,N-burnin+1))
for(it in burnin:N){
  X1[,,it-burnin+1] <- t(X[[it]][1,,])
  X2[,,it-burnin+1] <- t(X[[it]][2,,])
}
count=0
pb <- txtProgressBar(min=1,max=n*(n-1)/2,style=3)
for(i in 1:(n-1)){
  for( j in (i+1):n){ 
    kappaij[i,j] <- c_postzeroprob(X1[i,,],X2[i,,],X1[j,,],X2[j,,],s2[burnin:N],
                                   lambda,p0)
    kappaij[j,i] <- c_postzeroprob(X1[j,,],X2[j,,],X1[i,,],X2[i,,],s2[burnin:N],
                                   lambda,p0)
    count = count+1
    setTxtProgressBar(pb,count)
  }
}
close(pb)
Prob0 <- matrix(0,n,n)
for(i in 1:n){
  for(j in c(1:n)[-i]){
    Prob0[i,j] <- 1/(1+kappaij[i,j])
  }
}
NegInfl <- matrix(0,n,n)
for( i in 1:n){
  for(j in c(1:n)[-i]){
    temp <- numeric(TT-1)
    for(tt in 2:TT){
      temp[tt-1] <- atan2(EX[2,tt,j]-EX[2,tt-1,i],EX[1,tt,j]-EX[1,tt-1,i])-
        atan2(EX[2,tt,i]-EX[2,tt-1,i],EX[1,tt,i]-EX[1,tt-1,i])
    }
    if(mean(abs(temp))>pi/2){ Prob0[i,j] <- 1;NegInfl[i,j] <- 1;}#
  }
}
index <- which(apply(Prob0,1,function(x) sum(x < 0.5)-1)>0)
Results <- apply(Prob0[index,],1,function(x) return(which(x < 0.5)))
for(i in 1:length(index)){
  Results[[i]] <- Results[[i]][which(Results[[i]] != index[i])]
}
Refined <- list()
for(i in 1:length(index)){
  temp <- NULL
  for(j in 1:length(Results[[i]])){
    for(tt in 1:TT){
      dijt <- sqrt(t(EX[,tt,index[i]]-EX[,tt,Results[[i]][j]])%*%(EX[,tt,index[i]]-EX[,tt,Results[[i]][j]]))
      if(dijt < Ew[index[i]] | dijt < Ew[Results[[i]][j]]){
        temp <- c(temp,Results[[i]][j])
      }
      Refined[[i]]<-unique(temp)
    }
  }
}
print(index)
print(Refined)

conc <- list()
temp <- numeric(TT-1)
for( i in 1:length(index)){
  conc[[i]] <- matrix(0,TT-1,3)
  muX <- numeric(length(Refined[[i]]))
  Angle <- matrix(0,TT-1,length(muX))
  for(ind in 1:length(Refined[[i]])){
    for(tt in 2:TT){
      Angle[tt-1,ind] <- -atan2(EX[2,tt,Refined[[i]][ind]]-EX[2,tt-1,index[i]],
                                EX[1,tt,Refined[[i]][ind]]-EX[1,tt-1,index[i]])
      Rotate <- matrix(c(cos(Angle[tt-1,ind]),sin(Angle[tt-1,ind]),-sin(Angle[tt-1,ind]),cos(Angle[tt-1,ind])),2,2)
      temp[tt-1] <- (Rotate%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))[1]
    }
    muX[ind] <- max(mean(temp),0)
  }
  ind <- which.max(muX)
  conc[[i]][,3] <- rep(Refined[[i]][ind],TT-1)
  muX <- max(muX)
  for(tt in 2:TT){
    temp[tt-1] <- sqrt(t(EX[,tt,index[i]]-EX[,tt-1,index[i]])%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))
    conc[[i]][tt-1,1] <- muX*temp[tt-1]
  }
  conc[[i]][,1] <- conc[[i]][,1]/mean(s2[burnin:N])
  conc[[i]][,2] <- -Angle[,ind]
}

yl <- xl <- matrix(0,length(index),2)
for(i in 1:length(index)){
  tempX <- mean(range(EX[1,,c(Refined[[i]],index[i])]))
  tempY <- mean(range(EX[2,,c(Refined[[i]],index[i])]))
  rangeX <- range(EX[1,,c(Refined[[i]],index[i])])
  rangeY <- range(EX[2,,c(Refined[[i]],index[i])])
  Size <- max(diff(rangeX),diff(rangeY))
  xl[i,] <- tempX + c(-1,1)*Size/2
  yl[i,] <- tempY + c(-1,1)*Size/2
  xl[i,] <- xl[i,] + c(-1,1)*0.05*Size
  yl[i,] <- yl[i,] + c(-1,1)*0.05*Size
}

Theta <- seq(from=0,to=2*pi,length=2001)[-2001]
MAX <- NULL
for(i in 1:length(index)){
  for(tt in 1:(TT-1)){
    MAX <- c(MAX,dvonmises(Theta,mu=circular(conc[[i]][tt,2]),kappa=conc[[i]][tt,1]))
  }
}
MAX <- max(MAX)
RSc <- 0.05*apply(xl,1,diff)

for(i in 1:length(index)){
  cols <- "black"; 
  cols <- c(cols,rep("blue",length(Refined[[i]])));
  PCH <- 16;
  PCH <- c(PCH, rep(17,length(Refined[[i]])))
  
  temp <- NULL
  for(tt in 1:TT) temp <- rbind(temp,t(EX[,tt,c(index[i],Refined[[i]])]))
  jpeg(file=paste("NodalInfluence_Actor",index[i],".jpg",sep=""))
  par(mar=c(1,1,1,1)+0.1)
  plot(temp,col="white",ylab="",xlab="",xlim=xl[i,],ylim=yl[i,],xaxt="n",yaxt="n")
  for(tt in 1:TT){ 
    points(t(EX[,tt,c(index[i],Refined[[i]])]),pch=PCH,col=cols,cex=2)
    if(tt>1) arrows(EX[1,tt-1,c(index[i],Refined[[i]])],EX[2,tt-1,c(index[i],Refined[[i]])],
                    EX[1,tt,c(index[i],Refined[[i]])],EX[2,tt,c(index[i],Refined[[i]])],
                    length=0.2,col=cols)
    if(tt < TT){
      dTheta <- dvonmises(Theta,mu=circular(conc[[i]][tt,2]),kappa=conc[[i]][tt,1])
      ind <- 1000*(MAX-dTheta)/MAX; 
      if(min(ind)<1) ind[which.min(ind)]=1; 
      if(max(ind)>1000) ind[which.max(ind)]=1000
      #MAX=max(dTheta)
      #Cols <- gray(1:1000/1000)[1000*(max(dTheta)-dTheta)/max(dTheta)]
      Cols <- gray(1:1000/1000)[ind]
      tempSeq <- seq(from=1.75,to=0.000001,length=1000)
      #CEX <- tempSeq[1000*(max(dTheta)-dTheta)/max(dTheta)]
      CEX <- tempSeq[ind]
      par(new=TRUE)
      plot(x=EX[1,tt,index[i]]+RSc[i]*cos(Theta),y=EX[2,tt,index[i]]+RSc[i]*sin(Theta),
           col=Cols,type="b",pch=16,ylab="",xlab="",cex=CEX,
           ylim=yl[i,],xlim=xl[i,],xaxt="n",yaxt="n",bty="n")
    }
  }
  dev.off()
}

#---

Theta <- seq(from=-pi,to=pi,length=500)
conc <- list()
temp <- numeric(TT-1)
for( i in 1:length(index)){
  conc[[i]] <- numeric(length(Refined[[i]]))
  for(j in 1:length(conc[[i]])){
    ind <- Refined[[i]][j]
    for(tt in 2:TT){
      Angle <- -atan2(EX[2,tt,ind]-EX[2,tt-1,index[i]],EX[1,tt,ind]-EX[1,tt-1,index[i]])
      Rotate <- matrix(c(cos(Angle),sin(Angle),-sin(Angle),cos(Angle)),2,2)
      temp[tt-1] <- (Rotate%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))[1]
    }
    muX <- max(mean(temp),0)
    for(tt in 2:TT){
      temp[tt-1] <- sqrt((EX[,tt,index[i]]-EX[,tt-1,index[i]])%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))
    }
    conc[[i]][j] <- muX*mean(temp)/mean(s2[burnin:N])
  }
}
for(i in 1:length(index)){
  MAX <- dvonmises(0,mu=circular(0),kappa=max(conc[[i]]))
  
  jpeg(filename=
    paste("NodInflAggregateWrapped_Actor",index[i],".jpg",sep=""))
  par(mar=c(2, 2, 2, 2))
  for(j in 1:length(Refined[[i]])){
    dTheta <- dvonmises(Theta,mu=circular(0),kappa=conc[[i]][order(conc[[i]])][j])
    ind <- 1000*(MAX-dTheta)/MAX 
    if(min(ind)<1) ind[which.min(ind)]=1; 
    if(max(ind)>1000) ind[which.max(ind)]=1000
    tempSeq <- seq(from=3,to=0.0001,length=1000)
    CEX <- tempSeq[ind]
    #     Cols <- gray(1:2000/2000)[2000*(max(dTheta)-dTheta)/max(dTheta)]
    Cols <- gray(1:1000/1000)[ind]
    plot(j/length(Refined[[i]])*cos(Theta),j/length(Refined[[i]])*sin(Theta),
         col=Cols,type="b",pch=16,ylab="",xlab="",cex=CEX,
         ylim=c(-1.2,1.2),xlim=c(-1.2,1.2),xaxt="n",yaxt="n")
    if(j < length(Refined[[i]])) par(new=TRUE)    
  }
  text(x=c(1.2,1.1*cos(pi/4),1.1*cos(3*pi/4),1.1*cos(pi/4),1.1*cos(3*pi/4)),
       y=c(0,1.3*sin(pi/4),1.3*sin(3*pi/4),1.3*sin(-pi/4),1.3*sin(-3*pi/4)),
       labels=c(0,expression(frac(pi,4)),expression(frac(3*pi,4)),
                expression(-frac(pi,4)),expression(-frac(3*pi,4))))
  dev.off()
  jpeg(file=
    paste("NodInflAggregate_Actor",index[i],".jpg",sep=""))
  par(mar=c(4, 4, 2, 2))
  for(j in 1:length(Refined[[i]])){
    curve(dvonmises(x,mu=0,kappa=conc[[i]][order(conc[[i]],decreasing=TRUE)][j]),
          xlim=c(-pi,pi),col=1,
          xlab="",ylab="",xaxt="n",
          ylim=c(0.9*dvonmises(pi,mu=0,kappa=range(conc[[i]])[1]),
                 1.1*dvonmises(0,mu=0,kappa=range(conc[[i]])[2])))
    abline(h=0,v=0,lty=3)
    if(j < length(Refined[[i]])) par(new=TRUE)
  }
  axis(1,at=c(-pi,-pi/2,0,pi/2,pi),labels=c(expression(-pi),expression(-pi/2),
                                            0,expression(pi/2),expression(pi)))
  dev.off()
}


#----------------------------------------------------
###
###Prediction
###
X2T <- X1T <- matrix(0,n,N1)
for(it in (burnin+1):N){
  X1T[,it-burnin] <- X[[it]][1,TT,]
  X2T[,it-burnin] <- X[[it]][2,TT,]
}
yhat <- c.prediction(t(EX[,TT,]),s2[-c(1:burnin)],X1T,X2T,Bin[-c(1:burnin)],Bout[-c(1:burnin)],w[,-c(1:burnin)])


#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------


