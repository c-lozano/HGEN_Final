
# Part A ####
makePlots <- c(c=T,e=T,f=T) # to reduce RMD render time
runSimulations <- c(c=T,e=T,f=T)

standardize <- function(x){ # for standardizing G later
  if (is.matrix(x) && !any(dim(x)==1)){
    xStd <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    for (col in 1:ncol(x)){
      xStd[,col] <-(x[,col]-mean(x[,col]))
      if(sd(x[,col])!=0){
        xStd[,col] <- xStd[,col]/sd(x[,col])
      }
      # xStd[,col] <- x[,col]/sd(x[,col])
    }
  }
  else{
    xStd <- (x-mean(x))
    if(sd(x)!=0){
      xStd <- xStd/sd(x)
    }
    # xStd <- x/sd(x)
  }
  return(xStd)
}
unStandardize <- function(xStd,x){
  if (is.matrix(x) && !any(dim(x)==1)){
    xUnStd <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    for (col in 1:ncol(x)){
      xUnStd[,col] <-xStd[,col]*sd(x[,col]) + mean(x[,col])
      if(is.integer(x)){
        xUnStd[,col] <- as.integer(xUnStd[,col])
      }
    }
  }
  else{
    xUnStd <- xStd*sd(x) + mean(x)
    if(is.integer(x)){
      xUnStd <- as.integer(xUnStd)
    }
  }
  return(xUnStd)
}
gSim <- function(N,S,numCausal=500){ # simulate genotypes
  alFreq <- rbeta(S,4,4) # simulate allele frequencies
  G <- matrix(NA,nrow=N,ncol=S)
  for (j in 1:S){
    G[,j] <- rbinom(N,2,alFreq[j])
  }
  if(numCausal==0){varBeta <- 0}
  else{varBeta <- h2/numCausal}
  beta <- matrix(rep(0,S))
  causals <- sample(S,numCausal)
  beta[causals,] <- causalBetas <- rnorm(numCausal,0,sqrt(varBeta))
  
  return(list(G=G,b=beta,cB=causalBetas))
}

pSim <- function(G, h2,b,cB,numCausal=500){
  N <- nrow(G)
  S <- ncol(G)
  varE <- 1-h2
  E <- matrix(rnorm(N,0,sqrt(varE)))
  stdG <- standardize(G)
  geneticComponent <-stdG%*%b
  y <- geneticComponent + E
  h2Est <- var(geneticComponent)/var(y)
  return(list(stdG=stdG,b=b,cB=cB,genComp=geneticComponent,y=y,h2Est=h2Est))
}


genSim <- function(N,S,h2,numCausal=500){
  genetics <- gSim(N,S,numCausal)
  output <- pSim(genetics$G,h2,genetics$b,genetics$cB,numCausal)
  return(list(G=genetics$G,stdG=output$stdG,b=genetics$b,cB=genetics$cB,y=output$y,h2Est=output$h2Est))
}


# Part B ####

offspringSim <- function(G0,h2,b,cB,numCausal=500){
  N <- nrow(G0)
  S <- ncol(G0)
  G1 <- matrix(0,nrow=N,ncol=S)
  
  for (i in 1:N){
    parents <- sample(N,2)
    for (j in 1:S){
      for (p in 1:2){
        if (G0[parents,j][p]==1){
          G1[i,j] <- G1[i,j]+rbinom(1,1,0.5)
        }
        else if (G0[parents,j][p]==2){
          G1[i,j] <- G1[i,j]+1
        }
      }
    }
  }
  output <- pSim(G1,h2,b,cB,numCausal)
  return(list(G=G1,stdG=output$stdG,y=output$y,h2Est=output$h2Est))
}


# Part C ####

HEReg <- function(y,G){
  pSim <- y %*% t(y)
  K <- cor(t(G))
  Y <- pSim[upper.tri(pSim)] 
  X <- K[upper.tri(K)] 
  reg <- lm(Y ~ X) 
  return(reg$coefficients[2]) 
}


N <- 1e3
S <- 500
h2 <- 0.3


gens <- 10
sims <- 50
if(runSimulations['c']){
  h2Ests <- matrix(NA,nrow=gens,ncol=sims)
  HERegEsts <- matrix(NA,nrow=gens,ncol=sims)
  for (s in 1:sims){
    gen0 <- genSim(N,S,h2)
    h2Ests[1,s] <- gen0$h2Est
    HERegEsts[1,s] <- HEReg(gen0$y,gen0$stdG)
    prevGen <- gen0
    for (g in 2:gens){
      nextGen <- offspringSim(prevGen$G,h2,gen0$b,gen0$cB)
      h2Ests[g,s] <- nextGen$h2Est
      HERegEsts[g,s] <- HEReg(nextGen$y,nextGen$stdG)
      prevGen <- nextGen
      cat('\014')
      cat('progress: ',round(((g+(gens*(s-1)))/(gens*sims))*100),'%',sep='')
    }
  }
  h2EstMeans <- rowMeans(h2Ests)
  HERegEstMeans <- rowMeans(HERegEsts)
  write.csv(data.frame('True h2'=h2EstMeans,'HE Regression h2 Means'=HERegEstMeans), file='Non-AM h2 means.csv')
}

if(makePlots['c']){
  png('FigC.png',width=5,height=5,units='in',res=300)
  h2Means <- read.csv('Non-AM h2 means.csv',row.names=1)
  h2EstMeans <- h2Means$True.h2
  HERegEstMeans <- h2Means$HE.Regression.h2.Means
  plot(0,xlim=c(1,gens),ylim=c(0,1),type='n',xlab='Generations',ylab='',main=bquote(h^2~'without assortive mating'))
  mtext(bquote(h^2),side=2,line=2.5,crt=90,las=1)
  abline(h=h2,lty=2,lwd=2)
  lines(1:gens,h2EstMeans,lwd=2,col='blue')
  lines(1:gens,HERegEstMeans,col='red',lwd=2)
  legend('topright',legend=c(bquote('True'~h^2),bquote('HE regr.'~h^2~'approx.'),bquote('initial'~h^2)),col=c('blue','red','black'),lty=c(1,1,2),lwd=2)
  dev.off()
}


# Part D ####

require(MASS)
correlatedParents <- function(y,r,tol){
  if(dim(y)[2]!=1) y <- t(y)
  N <- nrow(y)
  parent1Inds <- sample(N,N/2)
  parent2Inds <- (1:N)[-parent1Inds]
  corrMat <- matrix(c(1,r,r,1),nrow=2)
  mvDat <- mvrnorm(N/2,c(0,0),Sigma=corrMat,empirical=TRUE)
  desiredP1Ranks <- rank(mvDat[,1],ties='first')
  desiredP2Ranks <- rank(mvDat[,2],ties='first')
  newP1Inds <- parent1Inds[order(y[parent1Inds,])[desiredP1Ranks]]
  newP2Inds <- parent2Inds[order(y[parent2Inds,])[desiredP2Ranks]]
  rho <- cor(y[newP1Inds,],y[newP2Inds,])
  iter <- 0
  while(rho<(r-tol) || rho>(r+tol)){
    parent1Inds <- sample(N,N/2)
    parent2Inds <- (1:N)[-parent1Inds]
    mvDat <- mvrnorm(N/2,c(0,0),Sigma=corrMat,empirical=TRUE)
    desiredP1Ranks <- rank(mvDat[,1],ties='first')
    desiredP2Ranks <- rank(mvDat[,2],ties='first')
    newP1Inds <- parent1Inds[order(y[parent1Inds,])[desiredP1Ranks]]
    newP2Inds <- parent2Inds[order(y[parent2Inds,])[desiredP2Ranks]]
    rho <- cor(y[newP1Inds,],y[newP2Inds,])
  }
  parentInds <- matrix(c(newP1Inds,newP2Inds),ncol=2)
  colnames(parentInds) <- c('parent1','parent2')
  return(parentInds)
}



offspringSimAM <- function(G0,y0,h2,b,cB,r,tol,numCausal=500){
  N <- nrow(G0)
  S <- ncol(G0)
  G1 <- matrix(0,nrow=N,ncol=S)
  
  parentSets <- replicate(2,correlatedParents(y0,r,tol))
  parentInds <- matrix(NA,nrow=N,ncol=2)
  colnames(parentInds) <- c('parent1', 'parent2')
  parentInds[1:(N/2),] <- parentSets[,,1]
  parentInds[((N/2)+1):N,] <- parentSets[,,2]
  
  for (i in 1:N){
    parents <- parentInds[i,]
    for (j in 1:S){
      for (p in 1:2){
        if (G0[parents,j][p]==1){
          G1[i,j] <- G1[i,j]+rbinom(1,1,0.5)
        }
        else if (G0[parents,j][p]==2){
          G1[i,j] <- G1[i,j]+1
        }
      }
    }
  }
  output <- pSim(G1,h2,b,cB,numCausal)
  return(list(G=G1,stdG=output$stdG,y=output$y,b=output$b,cB=output$cB,h2Est=output$h2Est))
}


# Part E ####

N <- 1e3
S <- 500
h2 <- 0.3
AMCorr <- 0.3
AMTolerance <- 0.01

gens <- 10
sims <- 50
if(runSimulations['e']){
  h2Ests <- matrix(NA,nrow=gens,ncol=sims)
  HERegEsts <- matrix(NA,nrow=gens,ncol=sims)
  for (s in 1:sims){
    gen0 <- genSim(N,S,h2)
    h2Ests[1,s] <- gen0$h2Est
    HERegEsts[1,s] <- HEReg(gen0$y,gen0$stdG)
    prevGen <- gen0
    for (g in 2:gens){
      nextGen <- offspringSimAM(prevGen$G,prevGen$y,h2,prevGen$b,prevGen$cB,AMCorr,AMTolerance)
      h2Ests[g,s] <- nextGen$h2Est
      HERegEsts[g,s] <- HEReg(nextGen$y,nextGen$stdG)
      prevGen <- nextGen
      cat('\014')
      cat('progress: ',round(((g+(gens*(s-1)))/(gens*sims))*100),'%',sep='')
    }
  }
  
  h2EstMeans <- rowMeans(h2Ests)
  HERegEstMeans <- rowMeans(HERegEsts)
  write.csv(data.frame('True h2'=h2EstMeans,'HE Regression h2 Means'=HERegEstMeans), file='AM h2 means.csv')
}

if(makePlots['e']){
  png('FigE.png',width=5,height=5,units='in',res=300)
  h2Means <- read.csv('AM h2 means.csv',row.names=1)
  h2EstMeans <- h2Means$True.h2
  HERegEstMeans <- h2Means$HE.Regression.h2.Means
  plot(0,xlim=c(1,gens),ylim=c(0,1),type='n',xlab='Generations',ylab='',main=bquote(h^2~'under assortive mating'))
  mtext(bquote(h^2),side=2,line=2.5,crt=90,las=1)
  abline(h=h2,lty=2,lwd=2)
  lines(1:gens,h2EstMeans,lwd=2,col='blue')
  lines(1:gens,HERegEstMeans,col='red',lwd=2)
  legend('topright',legend=c(bquote('True'~h^2),bquote('HE regr.'~h^2~'approx.'),bquote('initial'~h^2)),col=c('blue','red','black'),lty=c(1,1,2),lwd=2)
  dev.off()
}


# Part F ####

esses <- c(5e2,5e3)
propCausals <- c(1e-2,1e-1,1)
h2 <- 0.3
AMCorr <- 0.3
AMTolerance <- 0.01

gens <- 10
sims <- 10

if(runSimulations['f']){
  ticker <- 0
  prevPercent <- 0
  startTime <- prevTime <- Sys.time()
  h2EstMeansLowS <- matrix(NA,nrow=gens,ncol=length(propCausals))
  HERegEstMeansLowS <- matrix(NA,nrow=gens,ncol=length(propCausals))
  colnames(h2EstMeansLowS) <- colnames(HERegEstMeansLowS) <- c('low pC', 'mid pC', 'high pC')
  h2EstMeansHighS <- matrix(NA,nrow=gens,ncol=length(propCausals))
  HERegEstMeansHighS <- matrix(NA,nrow=gens,ncol=length(propCausals))
  colnames(h2EstMeansHighS) <- colnames(HERegEstMeansHighS) <- c('low pC', 'mid pC', 'high pC')
  h2EstMeans <- list(lowS=h2EstMeansLowS, highS=h2EstMeansHighS)
  HERegEstMeans <- list(lowS=HERegEstMeansLowS, highS=HERegEstMeansHighS)
  
  for (pC in 1:length(propCausals)){
    for (s in 1:length(esses)){
      numCausal <- as.integer(esses[s]*propCausals[pC])
      h2Ests <- matrix(NA,nrow=gens,ncol=sims)
      HERegEsts <- matrix(NA,nrow=gens,ncol=sims)
      for (sim in 1:sims){
        gen0 <- genSim(N,esses[s],h2,numCausal)
        h2Ests[1,sim] <- gen0$h2Est
        HERegEsts[1,sim] <- HEReg(gen0$y,gen0$stdG)
        prevGen <- gen0
        ticker <- ticker+1
        newPercent <- round(ticker/(sims*gens*length(esses)*length(propCausals))*100)
        newTime <- Sys.time()
        if (newPercent>prevPercent){
          span <- difftime(newTime,startTime,units='mins')
          eta <- (span/newPercent)*(100-newPercent)
          prevPercent <- newPercent
          cat('\014')
          cat('\n\nprogress: ',(newPercent),'%',sep='')
          cat('\nelapsed: ',as.numeric(span), ' mins',sep='')
          cat('\n','ETA: ',as.numeric(eta),' mins',sep='')
        }
        
        for (g in 2:gens){
          nextGen <- offspringSimAM(prevGen$G,prevGen$y,h2,prevGen$b,prevGen$cB,AMCorr,AMTolerance,numCausal)
          h2Ests[g,sim] <- nextGen$h2Est
          HERegEsts[g,sim] <- HEReg(nextGen$y,nextGen$stdG)
          prevGen <- nextGen
          
          ticker <- ticker+1
          newPercent <- round(ticker/(sims*gens*length(esses)*length(propCausals))*100)
          newTime <- Sys.time()
          if (newPercent>prevPercent){
            span <- difftime(newTime,startTime,units='mins')
            eta <- (span/newPercent)*(100-newPercent)
            prevPercent <- newPercent
            cat('\014')
            cat('\n\nprogress: ',(newPercent),'%',sep='')
            cat('\nelapsed: ',as.numeric(span), ' mins',sep='')
            cat('\n','ETA: ',as.numeric(eta),' mins',sep='')
          }
        }
      }
      h2EstMeans[[s]][,pC] <- rowMeans(h2Ests)
      HERegEstMeans[[s]][,pC] <- rowMeans(HERegEsts)
    }
  }
  write.csv(data.frame(h2EstMeans),file='AM f1 true h2.csv')
  write.csv(data.frame(HERegEstMeans),file='AM f1 HE Reg.csv')
}

if(makePlots['f']){
  png('FigF.png',width=12,height=6,units='in',res=300)
  par(mfrow=(c(1,2)))
  for (s in 1:length(esses)){
    cols <- palette.colors(3,'Dark2')
    h2EstMeans <- read.csv('AM f1 true h2.csv',row.names=1)
    HERegEstMeans <- read.csv('AM f1 HE Reg.csv',row.names=1)
    plot(0,xlim=c(1,gens),ylim=c(0,1),type='n',xlab='Generations',ylab='',main=paste0('S = ',esses[s]))
    mtext(bquote(h^2),side=2,line=2.5,crt=90,las=1)
    abline(h=h2,lty=3,lwd=3)
    for (pC in 1:length(propCausals)){
      lines(1:gens,h2EstMeans[,3*(s-1)+pC],col=cols[pC],lwd=3,lty=5)
      lines(1:gens,HERegEstMeans[,3*(s-1)+pC],col=cols[pC],lwd=3,lty=1)
    }
    if(s==1){
      legend('topright',legend=c(bquote('initial'~h^2),
                                 bquote('True'~h^2),
                                 bquote('HE regr.'~h^2~'approx.')),
             col='black',
             lty=c(3,5,1),
             lwd=c(3,2,2))
      
      legend('topleft',legend=c( propCausals[1],
                                 propCausals[2],
                                 propCausals[3]),
             fill=cols, border=cols,
             # lty=c(3,5,1,rep(1,3)),
             # lwd=c(3,rep(2,5)),
             title='Proportion SNPs causal:')
    }
  }
  dev.off()
}

# Wizard Frog ####

#                              .-----.
#                             /7  .  (
#                            /   .-.  \
#              :)           /   /   \  \
#                          / `  )   (   )
#                         / `   )   ).  \
#                       .'  _.   \_/  . |
#      .--.           .' _.' )`.        |
#     (    `---...._.'   `---.'_)    ..  \
#      \            `----....___    `. \  |
#       `.           _ ----- _   `._  )/  |
#         `.       /"  \   /"  \`.  `._   |
#           `.    ((O)` ) ((O)` ) `.   `._\
#             `-- '`---'   `---' )  `.    `-.
#                /                  ` \      `-.
#              .'                      `.       `.
#             /                     `  ` `.       `-.
#      .--.   \ ===._____.======. `    `   `. .___.--`     .''''.
#     ' .` `-. `.                )`. `   ` ` \          .' . '  8)
#    (8  .  ` `-.`.               ( .  ` `  .`\      .'  '    ' /
#     \  `. `    `-.               ) ` .   ` ` \  .'   ' .  '  /
#      \ ` `.  ` . \`.    .--.     |  ` ) `   .``/   '  // .  /
#       `.  ``. .   \ \   .-- `.  (  ` /_   ` . / ' .  '/   .'
#         `. ` \  `  \ \  '-.   `-'  .'  `-.  `   .  .'/  .'
#           \ `.`.  ` \ \    ) /`._.`       `.  ` .  .'  /
#            |  `.`. . \ \  (.'               `.   .'  .'
#         __/  .. \ \ ` ) \                     \.' .. \__
#  .-._.-'     '"  ) .-'   `.                   (  '"     `-._.--.
# (_________.-====' / .' /\_)`--..__________..-- `====-. _________)
#                  (.'(.'

