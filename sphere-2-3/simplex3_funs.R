#' define class `simplex3`
new_simplex3 <- function(x = double(), y = double(), z = double(), w = double()){
  stopifnot(is.double(c(x,y,z,w)) && all(abs(x + y + z + w - 1) < 1e-5))
  structure(cbind(x, y, z, w), class = "simplex3")
}

#' Constructor of class \code{simplex3} 
#' 
#' @param x A matrix with 4 columns.
#' @return n x 4 matrix.
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(100)
#' y <- (1-x)/6
#' z <- (1-x)/6*2
#' w <- (1-x)/6*3
#' s1 <- simplex3(cbind(x,y,z,w))
simplex3 <- function(x){
  return(new_simplex3(as.double(x[,1]), 
                      as.double(x[,2]),
                      as.double(x[,3]),
                      as.double(x[,4])))
}

## plot methods for simplex3 ---------------------------------------------------

simplex_axis <- function(sub=TRUE, labs=NULL, lab.cols=NULL){
  p1 <- c(1/sqrt(3),0,0)
  p2 <- c(-sqrt(3)/6,1/2,0)
  p3 <- c(-sqrt(3)/6,-1/2,0)
  p4 <- c(0,0,sqrt(6)/3)
  plot3d(rbind(p1,p2,p3,p4,p1), type = "l", axes = FALSE,xlab="",ylab="",zlab="")
  lines3d(rbind(p4,p2))
  lines3d(rbind(p3,p1))
  move <- c(0,0,0.02)
  if(is.null(lab.cols)) cols <- rep("black", 4)
  if(!is.null(labs)) text3d(rbind(p1-move,p2-move,p3-move,p4+move),texts=labs,col=lab.cols)
  ## sublines
  if(sub){
    for(k in 0:9){
      for(i in 1:(9-k)){
        lines3d(rbind((i*p1+(10-k-i)*p2+k*p4)/10, ((10-k-i)*p2+i*p3+k*p4)/10),col="grey")
        lines3d(rbind((i*p2+(10-k-i)*p3+k*p4)/10, ((10-k-i)*p3+i*p1+k*p4)/10),col="grey")
        lines3d(rbind((i*p3+(10-k-i)*p1+k*p4)/10, ((10-k-i)*p1+i*p2+k*p4)/10),col="grey")
      }
    }
    for(i in 1:9){
      lines3d(rbind((i*p1+(10-i)*p4)/10, ((10-i)*p4+i*p3)/10),col="grey")
      lines3d(rbind((i*p2+(10-i)*p4)/10, ((10-i)*p4+i*p1)/10),col="grey")
      lines3d(rbind((i*p3+(10-i)*p4)/10, ((10-i)*p4+i*p2)/10),col="grey")
    }
  }
}

plot.simplex3 <- function(data, main=NULL, labs=NULL, lab.cols=NULL, black=NULL, red=NULL, gray=NULL){
  n <- nrow(data)
  simplex_axis(sub=F, labs, lab.cols)
  
  p1 <- c(1/sqrt(3),0,0)
  p2 <- c(-sqrt(3)/6,1/2,0)
  p3 <- c(-sqrt(3)/6,-1/2,0)
  p4 <- c(0,0,sqrt(6)/3)
  
  ## data
  for(i in 1:n){
    coord <- data[i,1]*p1+data[i,2]*p2+data[i,3]*p3+data[i,4]*p4
    points3d(coord[1],coord[2],coord[3],size=4)
  }
  
  if(!is.null(black)){
    coord <- matrix(ncol=3,nrow=nrow(black))
    for(i in 1:n){
      coord[i,] <- black[i,1]*p1+black[i,2]*p2+black[i,3]*p3+black[i,4]*p4
    }
    lines3d(coord[,1],coord[,2],coord[,3],lwd=5,col="black")
  }
  if(!is.null(red)){
    coord <- matrix(ncol=3,nrow=nrow(red))
    for(i in 1:n){
      coord[i,] <- red[i,1]*p1+red[i,2]*p2+red[i,3]*p3+red[i,4]*p4
    }
    lines3d(coord[,1],coord[,2],coord[,3],lwd=5,col="red")
  }
  if(!is.null(gray)){
    coord <- matrix(ncol=3,nrow=nrow(gray))
    for(i in 1:n){
      coord[i,] <- gray[i,1]*p1+gray[i,2]*p2+gray[i,3]*p3+gray[i,4]*p4
    }
    lines3d(coord[,1],coord[,2],coord[,3],lwd=5,col="gray")
  }
  if(!is.null(main)){
    text3d(p4+5*move, texts=main)
  }
  
  # wid <- rglwidget(width = 1000, height = 1000)
  # return(wid)
}
