# DepletionSimulation
code from "Effects of resource spatial distribution and tow overlap in the accuracy and precision of common methods used to estimate dredge efficiency for bottom trawl fisheries"

# Visualization of the depletion experiment
This animation shows how the sequence of tows with a dredge efficiency of 0.6 depletes the scallop bed and the associated relationship between the capture of each survey haul and the cumulative capture of the experiment.
![](animation1.gif)

# loading libraries
```
library(Rcpp)
library(raster)
library(rgeos)
library(raptr)
library(PBSmapping)
library(FSA)
```

# Compile C++ code to create a virtual landscape with varying degrees of spatial aggregation of suitable sites
Code modified from Hiebeler (2000) POPULATIONS ON FRAGMENTED LANDSCAPES WITH SPATIALLY STRUCTURED HETEROGENEITIES: LANDSCAPE GENERATION AND LOCAL DISPERSAL
https://doi.org/10.1890/0012-9658(2000)081[1629:POFLWS]2.0.CO;2

```
sourceCpp("Landscape.cpp")
```
# Functions to obtain the proportion of swept area 1, 2, .... n times
```
PropBarridosRaster = function(Inter){
R=raster::raster(res=c(0.1,0.1), ext=extent(Inter), crs=proj4string(Inter))
Q<-fasterize::fasterize(sf=sf::st_as_sf(Inter), raster=R, field="Value",  fun='sum')
return(table(Q[])/length(which(Q[]>0)))
}
```

# Function to create the landscape
```
getPaisaje = function(d=300, p0=0.3, q00=0.3, niters=10000000, maxErr= 0.000001,
                      proj=utm, ext=extP, n=1000000){
  RR=raster::raster(nrow=d, ncol=d, crs=utm, ext=extP)
  L=CreateLandscape(xsize=d, ysize=d, p0=p0, q00=q00, iterates=niters, maxErr=maxErr)
  R=matrix(L, d, d)
  RR[]=R[]
  indi=which(RR[]==1)
  individuals = xyFromCell(RR,  sample(indi, size=n, replace=T))
  individuals[] = individuals[] + rnorm(2*n, 0, 1)[]
  individuals=cbind(1:nrow(individuals), individuals)
  colnames(individuals)=c("EID", "X", "Y")
  individuals = as.EventData(individuals)
  return(individuals)
}
```
# Function to generate survey hauls
```
genLances = function(proj=utm, extP=extP, sover=0.5, posError=5){
  s=150 # vertical gap to the border of the box
  posx0=NULL
  posx1=NULL
  posy0=NULL
  posy1=NULL
  lal=700 # mean length of tows
  
  for(barrido in 1:4){
    posx0=c(posx0, rnorm(40,mean=seq(extP@xmin+s, extP@xmax-s,length=40), sd=sover))
    posx1=c(posx1, rnorm(40,mean=seq(extP@xmin+s, extP@xmax-s,length=40), sd=sover))
    posy0=c(posy0, rnorm(40, mean= extP@ymin+s, sd=sover))
    posy1=c(posy1, posy0+rnorm(40,mean=lal, sd=sover))
  }
  
  j=1
  po = SpatialPoints(data.frame(xo=posx0[j], yo=posy0[j]))
  r=runif(1,0, posError)
  a=runif(1,-pi,pi)
  po2= SpatialPoints(data.frame(xo=posx0[j]+ r*cos(a), yo=posy0[j]+r*sin(a)))
  
  pf = SpatialPoints(data.frame(xf=posx1[j], yf=posy1[j]))
  r=runif(1, 0, posError)
  a=runif(1,-pi,pi)
  pf2 = SpatialPoints(data.frame(xf=posx1[j]+ r*cos(a), yf=posy1[j]+r*sin(a)))
  
  Lineas = as(rbind(po,pf), "SpatialLines")
  Lineas2 = as(rbind(po2,pf2), "SpatialLines")
  
  for(j in 2:length(posx0)){
    po = SpatialPoints(data.frame(xo=posx0[j], yo=posy0[j]))
    r=runif(1, 0, posError)
    a=runif(1,-pi,pi)
    po2= SpatialPoints(data.frame(xo=posx0[j]+ r*cos(a), yo=posy0[j]+r*sin(a)))
    
    pf = SpatialPoints(data.frame(xf=posx1[j], yf=posy1[j]))
    r=runif(1, 0, posError)
    a=runif(1,-pi,pi)
    pf2 = SpatialPoints(data.frame(xf=posx1[j]+ r*cos(a), yf=posy1[j]+r*sin(a)))
    Lineas = rbind(Lineas, as(rbind(po,pf), "SpatialLines"))
    Lineas2 = rbind(Lineas2, as(rbind(po2,pf2), "SpatialLines"))
  }
  proj4string(Lineas) = proj
  proj4string(Lineas2) = proj
  PolyLances = buffer(Lineas,  width=2.5/2, dissolve=F)
  PolyLances2 = buffer(Lineas2,  width=2.5/2, dissolve=F)
  ##############################################################################
  # shuffle tows to break sequence
  ##############################################################################
  # il = sample(1:length(Lineas), size=length(Lineas))
  # PolyLances=PolyLances[il,]
  # PolyLances2=PolyLances2[il]
  PolyLances$ID= 1:length(Lineas)
  PolyLances2$ID= 1:length(Lineas2)
  return(list(PolyLances, PolyLances2))
}
```
# Function to get captures from survey hauls paths
```
getCapturas=function(Lances, vieiras, Efi){
pbLances = convert2PolySet(Lances[[1]], n=length(Lances[[1]]))
Capturas = rep(NA, length(Lances[[1]]))
TablaR=list()

for(l in 1:length(Lances[[2]])){
    Inter = raster::intersect(Lances[[2]][l,], Lances[[2]][1:l,])
    Inter$Value=rep(1, length(Inter))
    TablaR[[l]] = PropBarridosRaster(Inter)
    vie=findPolys(vieiras, pbLances[pbLances$PID==l,])
    CR=rbinom(n=length(vie$EID), Efi, size=1)
    CF=vie$EID[which(CR==1)]
    Capturas[l]=length(CF)
    
    if(length(match(CF, vieiras$EID))>0)
      vieiras=vieiras[-match(CF, vieiras$EID),]
  }
 TablaP=matrix(0, length(TablaR), length(TablaR))
  
  for(k in 1:length(TablaR)){
    cols=as.numeric(names(TablaR[[k]]))
    TablaP[k,cols]=TablaR[[k]]
  }
 return(list(Capturas, TablaP))
}
```

# Function to estimate efficiency for the patch model
```
Ajuste=function(par){
    d=par[1]
    e=par[2]
    espe=rep(NA, length(Capturas))
    
    for(j in 1:length(Capturas))
      espe[j]= e*sum(a[j]*d*Tabla[j,]*(1-e)^(expo))
    
    s=sqrt(sum((log1p(Capturas)-log1p(espe))^2)/(length(Capturas)-1))
    SLL=-sum(dnorm(x=log1p(Capturas), mean=log1p(espe), sd=s, log = T))
    return(SLL)
  }
```
# Loop for simulating "Nreps" depletion experiments with different spatial aggregation and tow overlap
Output is saved to "SalidaSimu.txt"
```
utm = "+proj=utm +zone=21 +south +datum=WGS84"
extP=extent(c(600000, 601000, 5000000, 5001000))

sal=data.frame(q00=NA, sover=NA, rep=NA, Par=NA, Leslie=NA, DeLury=NA, Carle=NA, Soquete=NA,
               No.Leslie=NA, No.DeLury=NA, No.Carle= NA, D.Patch=NA)
sal=sal[-1,]
write.table(sal, "SalidaSimu.txt", quote=F, col.names = T, row.names = F, sep="\t", append=F)

Nreps=100

for(q00 in c(0.05, 0.3, 0.60, 0.95)){
  for(sover in c(30, 2, 0.5)){
    for(rep in 1:Nreps){

Efi=runif(1,0.2,0.8)
Lances=genLances(ext=extP, sover=2, posError=0)
vieiras=getPaisaje(n=4000000, q00=q00)
CapBarri=getCapturas(Lances, vieiras, Efi)
Capturas=CapBarri[[1]]
Tabla=CapBarri[[2]]

Esfuerzo = gLength(Lances[[1]], byid = TRUE)/2
Esfuerzo=Esfuerzo*60/(4*1.852)/1000
Leslie=depletion(catch=Capturas, effort=Esfuerzo, method="Leslie")

AA=area(aggregate(Lances[[1]]))
a = area(Lances[[1]])
q=Leslie$est['q',1]
el=q*mean(Esfuerzo)*AA/mean(a)
  
DeLury=depletion(catch=Capturas+1, effort=Esfuerzo, method="DeLury")
qb=DeLury$est['q',1]
ed=qb*mean(Esfuerzo)*AA/mean(a)
  
barridos=rep(1:4, each=40)
Catch=as.vector(by(Capturas, barridos, sum))
Carle=removal(Catch, method="CarleStrub", alpha=1, beta=1)
pCarle=Carle$est['p']

expo=0:(length(Capturas)-1)
lo=c(0.1,0.0001)
up=c(10, 0.9999)
  
start=runif(2, lo, up)
resu=optim(start, Ajuste, method="L-BFGS-B", lower=lo, upper=up, control=list(maxit=30000))
  
write.table(data.frame(q00=q00, sover=sover, rep=rep, Par=Efi, Leslie=median(el), DeLury=median(ed),
            Carle=pCarle, Patch=resu$par[2], No.Leslie=Leslie$est[1,1], No.DeLury=DeLury$est[1,1], No.Carle= Carle$est[1], D.Patch=resu$par[1]),
            "SalidaSimu.txt", quote=F, col.names = F, row.names = F, sep="\t", append=T)
gc(reset=TRUE)
}
}
}
```
A full script to explore the behavior of different methods for varying scallop clustering, tow overlap and positional error is available in the file [Simulations.R](Simulations.R)
