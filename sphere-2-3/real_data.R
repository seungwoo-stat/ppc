## real data
source("sphere_funs.R")
source("simplex3_funs.R")

library(readxl)
library(stringr)
library(dplyr)

## Data generation -------------------------------------------------------------
data_real <- tibble()

filenames <- list.files("./data")
for(i in seq_along(filenames)){
  x <- read_excel(paste0("./data/",filenames[i]), skip = 6, col_names = FALSE)[,c(1,4:6)]
  colnames(x) <- c("city","DP","LP","PP")
  m <- nrow(x)
  x[,2] <- sapply(1:m,function(j) str_extract(x[j,2],'(?<=\\().*(?=\\))'))
  x[,3] <- sapply(1:m,function(j) str_extract(x[j,3],'(?<=\\().*(?=\\))'))
  x[,4] <- sapply(1:m,function(j) str_extract(x[j,4],'(?<=\\().*(?=\\))'))
  x <- cbind(district=rep(substring(filenames[i],1,nchar(filenames[i])-5),m),x)
  data_real <- bind_rows(data_real,x)
}
data_real <- as.data.frame(data_real)
data_real$DP <- as.numeric(data_real$DP)/100
data_real$LP <- as.numeric(data_real$LP)/100
data_real$PP <- as.numeric(data_real$PP)/100
data_real$ETC <- 1-rowSums(data_real[,3:5])
write.csv(data_real,file="19pe.csv",row.names = FALSE,fileEncoding = "UTF-8")

## -----------------------------------------------------------------------------
data_real <- read.csv("./19pe.csv")
data_real_sphere <- sphere3(sqrt(data_real[,3:6]))
data_real_simplex <- simplex3(data_real_sphere^2)

# Frechet and PGA --------------------------------------------------------------
f.mean <- FRECHET(data_real_sphere)
geo.pga <- PGA(data_real_sphere)
pga.curve <- t(sapply(1:250, function(i) project.sphere3(data_real_sphere[i,], geo.pga)))
d <- matrix(0,ncol=250,nrow=250)
for(i in 1:250){
  for(j in 1:250){
    d[i,j] <- DIST.sphere3(pga.curve[i,],pga.curve[j,])
  }
}
fit <- isoMDS(d, k=1)$points
pga.curve <- pga.curve[order(fit),]

# PPC --------------------------------------------------------------------------
PPC_real <- PPC(data_real_sphere, lambda = 0.05) # 8 iters

# PLOT -------------------------------------------------------------------------
labs <- c("Democratic Party", "Liberty Party", "People's Party", "ETC")
lab.cols <- c("#004EA2","#C9151E","#006241","black")
plot(data_real_simplex, labs=labs, lab.cols=lab.cols, red = PPC_real[[1]]^2, gray = pga.curve^2)
plot_real <- rglwidget(width = 2000, height = 1500)
htmlwidgets::saveWidget(plot_real,"real_big.html")

## 2 plot in one panel
mfrow3d(1,2)
plot(data_real_simplex, labs=labs, lab.cols=lab.cols, red = PPC_real[[1]]^2, gray = pga.curve^2)
plot(data_real_simplex, labs=labs, lab.cols=lab.cols, red = PPC_real[[1]]^2, gray = pga.curve^2)
plot_real2 <- rglwidget(width = 2000, height = 1500)
htmlwidgets::saveWidget(plot_real2,"real_mfrow.html")

# Dim reduction ----------------------------------------------------------------
PPC.f <- PPC_real[[1]]
PPC.sigma <- PPC_real[[2]]
w <- PPC_real[[3]]
n <- 250

a <- vector(length = n-1)
for(i in seq_along(a)){
  a[i] <- DIST.sphere3(PPC.f[i,],PPC.f[i+1,])
}
a <- c(0,cumsum(a))/sum(a)

PPC.reduction <- sapply(1:n, function(i) sum(w[i,] * a))

districts <- unique(data_real[,1])
mean_districts <- vector(length=length(districts))
for(i in seq_along(districts)){
  mean_districts[i] <- mean(PPC.reduction[data_real[,1]==districts[i]])
}

library(ggplot2)
library(ggrepel)
df <- data.frame(y=rep(0,17),x=mean_districts)
rownames(df) <- districts

# plot(rep(0,17), mean_districts, pch=20, ylab = "Dimension reduction", xaxt='n', xlab="")
# text(rep(0,17), mean_districts, label= districts, pos=c(4,4,4,2, 2,4,4,4, 2,4,4,4, 2,4,4,4, 4))

## save as 3*9 plot
ggplot(df, aes(x, y, label=rownames(df))) + 
  geom_point(size=2) +   
  geom_text_repel(nudge_y = 0.2,  direction = 'x', angle = 45, vjust = 0, segment.size = 0.2) + 
  ylim(-0.1, 0.4) +
  xlim(0.00, 1.0) +
  theme_classic() +  
  labs(x="Dimension Reduction", y="") +
  theme(axis.line.y  = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.title.y = element_blank())
