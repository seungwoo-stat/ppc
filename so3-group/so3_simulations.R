source("so3_funs.R")

## Simulation 1 ----------------------------------------------------------------
set.seed(0)
n <- 100
sigma <- 0.05

true_t <- seq(0,pi,length.out = n)
t <- runif(n,0,pi)
data1 <- vector("list", length = n)
true1 <- vector("list", length = n)
for(i in 1:n){
  rn <- rnorm(3,0,sd=sigma)
  error <- matrix(c(0,rn[1],rn[2], 0,0,rn[3], 0,0,0),ncol=3,byrow=T)
  error <- error - t(error) # skew-symmetric
  if(frob.norm(error) >= sqrt(2)*pi) stop()
  
  data1[[i]] <- matrix(c(0,cos(t[i]),1, 0,0,-sin(t[i]), 0,0,0),ncol=3)
  data1[[i]] <- data1[[i]] - t(data1[[i]])
  data1[[i]] <- exp.so3(data1[[i]])
  data1[[i]] <- EXP.so3(data1[[i]], data1[[i]] %*% error)
  
  true1[[i]] <- matrix(c(0,cos(true_t[i]),1, 0,0,-sin(true_t[i]), 0,0,0),ncol=3)
  true1[[i]] <- true1[[i]] - t(true1[[i]])
  true1[[i]] <- exp.so3(true1[[i]])
}
data1 <- so3(data1)

## FRECHET MEAN and PGA
# PGA1 <- PGA(data1)
PGA1 <- list(FRECHET(data1), matrix(c(-0.3605197, 0.5272191, -1.732824e-16,
                                      -0.7810120, 0.2803752, -1.643185e-16,
                                      0.5099469, 0.8021408, -1.340173e-16), byrow=T, ncol=3))
project_pga1 <- project.so3(data1, PGA1)
# plot(data1, red = project_pga1, gray = list(PGA1[[1]],PGA1[[1]]))

PPC1 <- PPC(data1, lambda = 0.1, vis = F, geo.pga = PGA1)
plot(data1, red = PPC1[[1]], black = true1)
plot(so3(true1), red = PPC1[[1]], threeD = 10)
plot1 <- rglwidget(width=2500,height=600)
htmlwidgets::saveWidget(plot1, paste("set1.html"))

## Fit by PGA, for comparison
f <- project.so3(data1, PGA1)
d <- matrix(0,ncol=n,nrow=n)
for(i in 1:n){
  for(j in 1:n){
    d[i,j] <- DIST.so3(f[[i]],f[[j]])
  }
}
d <- d+1e-10
diag(d) <- 0
fit <- isoMDS(d, k=1)$points
f <- f[order(fit)]
plot(so3(true1), red = f, threeD = 10)

## Simulation 2 ----------------------------------------------------------------
set.seed(1)
n <- 100
sigma <- 0.05

true_t <- seq(0,pi,length.out = n)
t <- runif(n,0,pi)
data2 <- vector("list", length = n)
true2 <- vector("list", length = n)
for(i in 1:n){
  rn <- rnorm(3,0,sd=sigma)
  error <- matrix(c(0,rn[1],rn[2], 0,0,rn[3], 0,0,0),ncol=3,byrow=T)
  error <- error - t(error) # skew-symmetric
  if(frob.norm(error) >= sqrt(2)*pi) stop()
  
  data2[[i]] <- matrix(c(0,cos(t[i]),cos(3*t[i])/2, 0,0,sin(t[i]), 0,0,0),ncol=3)
  data2[[i]] <- data2[[i]] - t(data2[[i]])
  data2[[i]] <- exp.so3(data2[[i]])
  data2[[i]] <- EXP.so3(data2[[i]], data2[[i]] %*% error)
  
  true2[[i]] <- matrix(c(0,cos(true_t[i]),cos(3*true_t[i])/2, 0,0,sin(true_t[i]), 0,0,0),ncol=3)
  true2[[i]] <- true2[[i]] - t(true2[[i]])
  true2[[i]] <- exp.so3(true2[[i]])
}
data2 <- so3(data2)
plot(data2, red = true2)

## FRECHET MEAN and PGA
# PGA2 <- PGA(data2)
PGA2 <- list(FRECHET(data2), matrix(c(-0.04172129,  0.99911008, -2.498208e-16,
                                      -0.75669066, -0.02754767, -8.592065e-17,
                                      -0.65244048, -0.03194023, -7.207782e-17), byrow=T, ncol=3))
project_pga2 <- project.so3(data2, PGA2)
plot(data2, red = project_pga2, gray = list(PGA2[[1]],PGA2[[1]]))

PPC2 <- PPC(data2, lambda = 0.05, vis = T, geo.pga = PGA2)
plot(data2, red = PPC2[[1]], black = true2)

a <- vector(length = n-1)
for(i in seq_along(a)){
  a[i] <- DIST.so3(PPC2[[1]][[i]],PPC2[[1]][[i+1]])
}
a[a==0] <- 1e-5
a <- c(0,cumsum(a))/sum(a)
sample2 <- sapply(1:10,function(i) sum(a <= seq(0,1,length=10)[i])) # sample uniformly

plot(so3(true2), red = PPC2[[1]][sample2], threeD = 10)
plot2 <- rglwidget(width=2500,height=600)
htmlwidgets::saveWidget(plot2, paste("set2.html"))

## Fit by PGA, for comparison
f <- project.so3(data2, PGA2)
d <- matrix(0,ncol=n,nrow=n)
for(i in 1:n){
  for(j in 1:n){
    d[i,j] <- DIST.so3(f[[i]],f[[j]])
  }
}
d <- d+1e-10
diag(d) <- 0
fit <- isoMDS(d, k=1)$points
f <- f[order(fit)]
plot(so3(true2), red = f, threeD = 10)

## PLOTS -----------------------------------------------------------------------
d1 <- matrix(0,ncol=100,nrow=100)
for(i in 1:100){
  for(j in 1:100){
    d1[i,j] <- DIST.so3(project_pga1[[i]],project_pga1[[j]])
  }
}
d1 <- d1 + 1e-10
diag(d1) <- 0
fit1 <- isoMDS(d1,k=1)

d2 <- matrix(0,ncol=100,nrow=100)
for(i in 1:100){
  for(j in 1:100){
    d2[i,j] <- DIST.so3(project_pga2[[i]],project_pga2[[j]])
  }
}
d2 <- d2 + 1e-10
diag(d2) <- 0
fit2 <- isoMDS(d2,k=1)

par(mfrow=c(1,2))
plot(data1, red = PPC1[[1]], black = true1, gray = project_pga1[order(fit1$points)])
plot(data2, red = PPC2[[1]], black = true2, gray = project_pga2[order(fit2$points)])

