##################################################################################
# useful functions
require(sp)

model_primicerio=function(dati_filare,saveplot=FALSE,titlefig,namefig,treshshold=40,span_value=0.25) {
  if (length(which(dati_filare$MISSING==1)) == 0) { stop("no missing data!")}
  if (saveplot)   { png(filename = namefig)}
  model=loess(Area ~ PIANTA,span=span_value,data=as.data.frame(dati_filare))
  plot(dati_filare$Area,xlab="Indice Filare",ylab="Feature Area")
  points(predict(model),col="red")
  points(dati_filare$MISSING*mean(dati_filare$Area),col="blue")
  abline(v=which(dati_filare$MISSING==1),col="green")
  if (saveplot) {dev.off()}
  
  candidate=which(as.vector(-model$residuals)>treshshold)
  missing=which(dati_filare$MISSING==1)
  missing_vector=dati_filare$MISSING
  model_residuals=as.vector(-model$residuals)
  vector=rep(0,length(dati_filare$MISSING))
  vector[candidate]=1
  
  res=list(candidate=candidate,
           missing=missing,
           precision=length(intersect(missing,candidate))/length(missing),
           vector=vector,
           missing_vector=missing_vector,
           model_residuals=model_residuals
  )
  
  return(res)
}


localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  return(y)
}

# https://anomaly.io/anomaly-detection-twitter-r/

findPeaks <- 
  function (x, thresh=0.05, span=0.25)
  {
    n <- length(x)
    y <- x
    mu.y.loc <- y
    y.loess <- loess(x~I(1:n), span=span)
    y <- y.loess[[2]]
    sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
    DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
    pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
   return(pks)
  }


getEdges <- function(x) {
  stopifnot(class(x) == "SpatialPolygons")
  lapply(x@polygons, function(y) {
    y@Polygons[[1]]@coords
  })
}


# Splining a polygon.
#
#   The rows of 'xy' give coordinates of the boundary vertices, in order.
#   'vertices' is the number of spline vertices to create.
#              (Not all are used: some are clipped from the ends.)
#   'k' is the number of points to wrap around the ends to obtain
#       a smooth periodic spline.
#
#   Returns an array of points. 
# 

spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}
##################################################################################
