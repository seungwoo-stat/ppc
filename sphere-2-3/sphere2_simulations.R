source("sphere_funs.R")

## Simulation 1 ----------------------------------------------------------------
## ARC dataset
set.seed(2)
n <- 100
to <- 2*pi*3/4
theta1 <- runif(n, 0, to)
phi1 <- rep(pi/4,n) #+ rnorm(n, 0, 0.1) # sigma = 0.07 also considered
data1 <- sphere2(sph2car(theta1,phi1,radius=1,deg=FALSE))
for(i in 1:n){
  pp <- PARALLEL(sphere2(0,0,1), data1[i,], c(rnorm(2,0,0.1),0))
  data1[i,] <- EXP(data1[i,],pp)
}

## PGA
fm1 <- FRECHET(data1)
geo.pga1 <- PGA(data1)
plot(data1)
text3d(fm1[1], fm1[2], fm1[3], texts = "f", col = "red")
geodesic_full.sphere2(geo.pga1, col = "red")

## PPC
PPC1 <- PPC(data1, lambda=0.1, max.iter = 20, vis = F)

## Simulation 2 ----------------------------------------------------------------
## WAVE dataset
# data generation
n <- 100
set.seed(2)
theta2 <- runif(n,0,to)
phi2 <- sin(3*theta2)/3
data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
for(i in 1:n){
  pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,0.1),0))
  data2[i,] <- EXP(data2[i,],pp)
}

## PGA
fm2 <- FRECHET(data2)
geo.pga2 <- PGA(data2)
plot(data2)
text3d(fm2[1], fm2[2], fm2[3], texts = "f", col = "red")
geodesic_full.sphere2(rbind(c(1,0,0),c(0,1,0)), col = "red",lwd=3)

## PPC
PPC2 <- PPC(data2, lambda=0.003, max.iter = 20, vis = T, geo.pga = rbind(c(1,0,0),c(0,1,0)))

## PLOTS -----------------------------------------------------------------------
## (1 PPC)
plot(data1)
lines3d(sph2car(seq(0,to,length.out = 100),phi1,radius=1,deg=FALSE),lwd=3)
lines3d(PPC1[[1]], col = "red", lwd=3)
geodesic_full.sphere2(geo.pga1, col = "darkgray", lwd=5)
plot1 <- rglwidget(width=1000,height=1000)
htmlwidgets::saveWidget(plot1, paste("data1.html"))

## (2 PPC)
plot(data2)
t <- seq(0,to,length.out = 200)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(PPC2[[1]], col = "red", lwd=3)
geodesic_full.sphere2(rbind(c(1,0,0),c(0,1,0)), col = "darkgray", lwd=5)
plot2 <- rglwidget(width=1000,height=1000)
htmlwidgets::saveWidget(plot2, paste("data2.html"))

## (1-2 PPC)
mfrow3d(1,2, sharedMouse = TRUE)
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data1[,1],data1[,2],data1[,3]), size=6)
lines3d(sph2car(seq(0,to,length.out = 100),phi1,radius=1,deg=FALSE),lwd=3)
lines3d(PPC1[[1]], col = "red", lwd=3)
geodesic_full.sphere2(geo.pga1, col = "darkgray", lwd=5)
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
t <- seq(0,to,length.out = 200)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(PPC2[[1]], col = "red", lwd=3)
geodesic_full.sphere2(rbind(c(1,0,0),c(0,1,0)), col = "darkgray", lwd=5)
plot12 <- rglwidget(width=2000,height=1000)
htmlwidgets::saveWidget(plot12, paste("data_merge.html"))

## (1-2 PPC)
hol <- sphere2(rbind(c(0,0,1),c(1,0,0),c(0,1,0),c(0,0,1)))
(hol.unwrap <- UNROLL(hol))
plot(hol)
geodesics.sphere2(hol, lwd=3)
lines3d(hol.unwrap + matrix(rep(hol[1,],4),byrow=T,ncol=3),lwd=3,col="red")

## HPC (2) ---------------------------------------------------------------------
library(spherepc)
data2_hau <- car2sph(data2[,1],data2[,2],data2[,3])[,1:2]
PPC2_hau_3 <- SPC.Hauberg(data2_hau, q=0.03, deletePoints = TRUE)
PPC2_hau_4 <- SPC.Hauberg(data2_hau, q=0.04, deletePoints = TRUE)
PPC2_hau_5 <- SPC.Hauberg(data2_hau, q=0.05, deletePoints = TRUE)
PPC2_hau_6 <- SPC.Hauberg(data2_hau, q=0.06, deletePoints = TRUE)
PPC2_hau_car_3 <- sph2car(PPC2_hau_3$prin.curves[,1],PPC2_hau_3$prin.curves[,2])[1:78,]
PPC2_hau_car_4 <- sph2car(PPC2_hau_4$prin.curves[,1],PPC2_hau_4$prin.curves[,2])[1:79,]
PPC2_hau_car_5 <- sph2car(PPC2_hau_5$prin.curves[,1],PPC2_hau_5$prin.curves[,2])[1:80,]
PPC2_hau_car_6 <- sph2car(PPC2_hau_6$prin.curves[,1],PPC2_hau_6$prin.curves[,2])[1:81,]

plot(data2, size = 6)
lines3d(PPC2[[1]], col = "red", lwd=3)
q_3 <- lines3d(PPC2_hau_car_3, lwd=3, col = "blue", alpha = 0.2)
q_4 <- lines3d(PPC2_hau_car_4, lwd=3, col = "blue", alpha = 0.4)
q_5 <- lines3d(PPC2_hau_car_5, lwd=3, col = "blue", alpha = 0.6)
q_6 <- lines3d(PPC2_hau_car_6, lwd=3, col = "blue", alpha = 0.8)

plot_hau_toggle <- rglwidget(width=2000,height=1000) %>%
  toggleWidget(ids = q_3) %>%
  toggleWidget(ids = q_4) %>%
  toggleWidget(ids = q_5) %>%
  toggleWidget(ids = q_6) %>%
  asRow(last = 4)
htmlwidgets::saveWidget(plot_hau_toggle, paste("hauberg_toggle.html"))

## Scales ----------------------------------------------------------------------

x=c(0,0.9)
y=c(0,0)
library(scales)
plot(x,y,type="l",col=alpha("blue",.6),xlim=c(0,3),lwd=3)
lines(x,y-0.1,col=alpha("blue",.8),lwd=3)
lines(x,y+0.1,col=alpha("blue",.4),lwd=3)
lines(x,y+0.2,col=alpha("blue",.2),lwd=3)
text(1,y+0.2,"q=0.03",pos=4)
text(1,y+0.1,"q=0.04",pos=4)
text(1,y-0.0,"q=0.05",pos=4)
text(1,y-0.1,"q=0.06",pos=4)

## smoothness plot -------------------------------------------------------------
library(scales)
# PPC1_1 <- PPC(data1, lambda=0.1, max.iter = 20, vis = F)
PPC1_2 <- PPC(data1, lambda=0.05, max.iter = 20, vis = F)
# PPC1_3 <- PPC(data1, lambda=0.01, max.iter = 20, vis = F)
# PPC1_4 <- PPC(data1, lambda=0.005, max.iter = 20, vis = F)
# PPC1_5 <- PPC(data1, lambda=0.001, max.iter = 20, vis = F)
PPC1_6 <- PPC(data1, lambda=0.0005, max.iter = 20, vis = F)

mfrow3d(1,2, sharedMouse = TRUE)
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data1[,1],data1[,2],data1[,3]), size=6)
lines3d(PPC1_2[[1]], col = "red", alpha=1, lwd=5)
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data1[,1],data1[,2],data1[,3]), size=6)
lines3d(PPC1_6[[1]], col = "red", alpha=1, lwd=5)
plot_varsmooth <- rglwidget(width=2000,height=1000)
htmlwidgets::saveWidget(plot_varsmooth, paste("various_lambda.html"))

## PLOT ALL --------------------------------------------------------------------
library(RColorBrewer)

n <- 100
set.seed(2)
theta2 <- runif(n,0,to)
phi2 <- sin(3*theta2)/3
data2 <- sphere2(sph2car(theta2,phi2,radius=1,deg=FALSE))
for(i in 1:n){
  pp <- PARALLEL(sphere2(0,0,1),data2[i,],c(rnorm(2,0,0.1),0))
  data2[i,] <- EXP(data2[i,],pp)
}

PPC2_1 <- PPC(data2,lambda = 0.03, vis = F, geo.pga = rbind(c(1,0,0),c(0,1,0)))
PPC2_2 <- PPC(data2,lambda = 0.003, vis = F, geo.pga = rbind(c(1,0,0),c(0,1,0)))
PPC2_3 <- PPC(data2,lambda = 0.0003, vis = F, geo.pga = rbind(c(1,0,0),c(0,1,0)))

SPC2_1 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.02, deletePoints = TRUE)$prin.curves
SPC2_2 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.03, deletePoints = TRUE)$prin.curves
SPC2_3 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.04, deletePoints = TRUE)$prin.curves
SPC2_4 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.05, deletePoints = TRUE)$prin.curves
SPC2_5 <- SPC(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.06, deletePoints = TRUE)$prin.curves

SPC2_1 <- unique(SPC2_1)
SPC2_1 <- sph2car(SPC2_1[,1],SPC2_1[,2])
SPC2_2 <- unique(SPC2_2)
SPC2_2 <- sph2car(SPC2_2[,1],SPC2_2[,2])
SPC2_3 <- unique(SPC2_3)
SPC2_3 <- sph2car(SPC2_3[,1],SPC2_3[,2])
SPC2_4 <- unique(SPC2_4)
SPC2_4 <- sph2car(SPC2_4[,1],SPC2_4[,2])
SPC2_5 <- unique(SPC2_5)
SPC2_5 <- sph2car(SPC2_5[,1],SPC2_5[,2])

HPC2_1 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.02, deletePoints = TRUE)$prin.curves
HPC2_2 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.03, deletePoints = TRUE)$prin.curves
HPC2_3 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.04, deletePoints = TRUE)$prin.curves
HPC2_4 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.05, deletePoints = TRUE)$prin.curves
HPC2_5 <- SPC.Hauberg(car2sph(data2[,1],data2[,2],data2[,3])[,1:2], q=0.06, deletePoints = TRUE)$prin.curves

HPC2_1 <- unique(HPC2_1)
HPC2_1 <- sph2car(HPC2_1[,1],HPC2_1[,2])
HPC2_2 <- unique(HPC2_2)
HPC2_2 <- sph2car(HPC2_2[,1],HPC2_2[,2])
HPC2_3 <- unique(HPC2_3)
HPC2_3 <- sph2car(HPC2_3[,1],HPC2_3[,2])
HPC2_4 <- unique(HPC2_4)
HPC2_4 <- sph2car(HPC2_4[,1],HPC2_4[,2])
HPC2_5 <- unique(HPC2_5)
HPC2_5 <- sph2car(HPC2_5[,1],HPC2_5[,2])

PNS2 <- PrincipalCircle(car2sph(data2[,1],data2[,2],data2[,3])[,1:2])
PNS2 <- GenerateCircle(PNS2[1:2], radius = PNS2[3])
PNS2 <- sph2car(PNS2[,1],PNS2[,2])

PGA2 <- PGA.sphere2(data2)

#1
mfrow3d(2,3, sharedMouse = TRUE)
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
t <- seq(0,to,length.out = 200)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(PPC2_1[[1]], col = brewer.pal(7, "RdYlGn")[1], lwd=3)
lines3d(PPC2_2[[1]], col = brewer.pal(7, "RdYlGn")[2], lwd=3)
lines3d(PPC2_3[[1]], col = brewer.pal(7, "RdYlGn")[3], lwd=3)

#2
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(SPC2_1[1:76,], col = brewer.pal(11, "BrBG")[7], lwd=3)
lines3d(SPC2_2[1:79,], col = brewer.pal(11, "BrBG")[8], lwd=3)
lines3d(SPC2_3[1:80,], col = brewer.pal(11, "BrBG")[9], lwd=3)
lines3d(SPC2_4[1:81,], col = brewer.pal(11, "BrBG")[10], lwd=3)
lines3d(SPC2_5[1:82,], col = brewer.pal(11, "BrBG")[11], lwd=3)

#3
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(HPC2_1, col = brewer.pal(11, "BrBG")[7], lwd=3)
lines3d(HPC2_2, col = brewer.pal(11, "BrBG")[8], lwd=3)
lines3d(HPC2_3, col = brewer.pal(11, "BrBG")[9], lwd=3)
lines3d(HPC2_4, col = brewer.pal(11, "BrBG")[10], lwd=3)
lines3d(HPC2_5, col = brewer.pal(11, "BrBG")[11], lwd=3)

#4
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
lines3d(PNS2, col = "red", lwd=3)

#5
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
rgl.sphpoints(car2sph(data2[,1],data2[,2],data2[,3]), size=6)
lines3d(sph2car(t,sin(3*t)/3,radius=1,deg=FALSE),lwd=3)
geodesic_full.sphere2(PGA2, col = "red", lwd=3)

#6
next3d()
rgl.sphgrid(col.long='white', col.lat='white', radaxis = FALSE, add= TRUE)
lines3d(sph2car(seq(0,to,length.out = 100),phi1,radius=1,deg=FALSE),lwd=3)
rgl.sphpoints(car2sph(data1[,1],data1[,2],data1[,3]), size=6)
lines3d(PPC1_2[[1]], col = "red", lwd=3)

#save
plot_compare <- rglwidget(width=2000,height=1000)
htmlwidgets::saveWidget(plot_compare, paste("compare.html"))

plot(c(0,0))
legend("topright", paste0("q=0.0",2:6),col=brewer.pal(11, "BrBG")[7:11],lty=1,lwd=3)
legend("topleft", c(expression(paste(lambda,"=",0.0003)),
                    expression(paste(lambda,"=",0.003)),
                    expression(paste(lambda,"=",0.03))),col=brewer.pal(7, "RdYlGn")[3:1],lty=1,lwd=3)
