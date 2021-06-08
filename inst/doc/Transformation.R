## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load libraries, echo = T-------------------------------------------------
library(ecotraj)

## -----------------------------------------------------------------------------
sites = rep(1,4)
surveys=1:4
spdata = rbind(c(35,30,20,15),
               c(50,25,15,10),
               c(65,20,10,5),
               c(80,15,5,0))

## -----------------------------------------------------------------------------
D = vegan::vegdist(spdata, "bray")
is.metric(D)
D

## ---- echo=TRUE, fig = TRUE, fig.height=4, fig.width=6, fig.align="center"----
trajectoryPCoA(D,sites,surveys, survey.labels = T)

## -----------------------------------------------------------------------------
trajectoryLengths(D,sites,surveys)
trajectoryAngles(D,sites,surveys)
trajectoryDirectionality(D,sites,surveys)

## -----------------------------------------------------------------------------
sqrtD = sqrt(D)
sqrtD


## ---- echo=TRUE, fig = TRUE, fig.height=4, fig.width=6, fig.align="center"----
trajectoryPCoA(sqrtD,sites,surveys, survey.labels = T)

## -----------------------------------------------------------------------------
trajectoryLengths(sqrtD,sites,surveys)
trajectoryAngles(sqrtD,sites,surveys)
trajectoryAngles(sqrtD,sites,surveys, all=TRUE)
trajectoryDirectionality(sqrtD,sites,surveys)

## -----------------------------------------------------------------------------
Nsteps = 50
CC = 50
Nreplace <- CC*0.05

## -----------------------------------------------------------------------------
x <- c(0, 1, 0, 67, 1, 3, 0, 2, 2, 2, 1, 6, 2, 0, 0, 2, 5, 1, 6, 0)
poffspring <- c(0, 0, 0.002, 0.661 ,0 ,0, 0.037, 0.281, 0, 0, 0, 0.008, 0, 0, 0.005, 0.003, 0, 0, 0, 0)

## -----------------------------------------------------------------------------
m <- matrix(0, nrow=Nsteps+1, ncol=length(x))
m[1, ] = x
for(k in 1:Nsteps) {
  pdeath <-x/sum(x) #Equal probability of dying
  deaths<-rmultinom(1,Nreplace, pdeath)
  x <- pmax(x - deaths,0)
  offspring = rmultinom(1,Nreplace, as.vector(poffspring))
  x <- x + offspring
  m[k+1, ]<-x
}

## -----------------------------------------------------------------------------
Sj <- seq(1,Nsteps+1, by=4) #Sample every four steps
mj <- m[Sj,]
surveys = 1:length(Sj)
sites = rep(1,length(Sj))

## -----------------------------------------------------------------------------
D <- vegan::vegdist(mj,"bray")

## -----------------------------------------------------------------------------
is.metric(D, tol=0.0000001)

## ---- echo=TRUE, fig = TRUE, fig.height=5, fig.width=6, fig.align="center"----
pcoa<-trajectoryPCoA(D, sites, surveys, selection=1,length=0.1, axes=c(1,2), survey.labels = T)
pcoaD = dist(pcoa$points)

## ---- echo=TRUE, fig = TRUE, fig.height=5, fig.width=6, fig.align="center"----
sqrtD = sqrt(D)
pcoaSqrt = trajectoryPCoA(sqrtD, sites, surveys, selection=1,length=0.1, axes=c(1,2), survey.labels = T)

## ---- echo=TRUE, fig = TRUE, fig.height=5, fig.width=6, fig.align="center"----
res <- smacof::mds(D, ndim = length(Sj)-1, type = "interval")
mmdsD <- dist(res$conf)
trajectoryPlot(res$conf, sites, surveys, selection=1,length=0.1, axes=c(1,2), survey.labels = T)

## -----------------------------------------------------------------------------
smacof::stress0(D,pcoaSqrt$points, type="interval")
smacof::stress0(D,pcoa$points, type="interval")
smacof::stress0(D,res$conf, type="interval")

## -----------------------------------------------------------------------------
anglesD = trajectoryAngles(D,sites,surveys)
anglesSqrtD = trajectoryAngles(sqrtD,sites,surveys)
anglesPcoaD = trajectoryAngles(pcoaD,sites,surveys)
anglesmmdsD = trajectoryAngles(mmdsD,sites,surveys)

df<-as.data.frame(rbind(anglesD, anglesSqrtD, anglesPcoaD, anglesmmdsD))
row.names(df)<-c("local", "global.sqrt", "global.pcoa", "global.mmds")
round(df,2)

## -----------------------------------------------------------------------------
trajectoryDirectionality(D,sites,surveys)
trajectoryDirectionality(sqrtD,sites,surveys)
trajectoryDirectionality(pcoaD,sites,surveys)
trajectoryDirectionality(mmdsD,sites,surveys)

