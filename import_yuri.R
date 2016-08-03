#NB: All the series must run from 1749 to 2012 (incl. NA)
load("slp.Rdata")
load("t.Rdata")
load("r.Rdata")
load("ghcn_temp.Rdata")
load("ghcn_precip.Rdata")

x.data <- cbind(slp$data,t$data,r$data,ghcn_temp$data,ghcn_precip$data)
x.lon <- c(slp$lon,t$lon,r$lon,ghcn_temp$lon,ghcn_precip$lon)
x.lat <- c(slp$lat,t$lat,r$lat,ghcn_temp$lat,ghcn_precip$lat)
x.name <- c(slp$name,t$name,r$name,ghcn_temp$name,ghcn_precip$name)
x.height <- c(slp$height,t$height,r$height,ghcn_temp$height,ghcn_precip$height)
x.time <- slp$time


### REMOVE SERIES TOO CLOSE TO EACH OTHERS KEEPING THE ONE WITH MORE DATA
x.min <- 1.0 #min long. distance (degrees)
y.min <- 1.0 #min lat. distance (degrees)
n_new <- n_stations
n_r <- 0
i <- 1
while (i<=n_new) {
  j <- 1
  while (j<=n_new) {
    a <- abs(x.lon-x.lon[i])
    b <- abs(x.lat-x.lat[i])
    if ((a[j]<x.min)&(b[j]<y.min)&(i!=j)&(x.name[i]==x.name[j])) {
      if (sum(!is.na(x.data[,i]))<=sum(!is.na(x.data[,j]))) {
        print(paste(x.lon[i],x.lat[i],sum(!is.na(x.data[,i])),"values: removed"))
        print(paste(x.lon[j],x.lat[j],sum(!is.na(x.data[,j])),"values: kept"))
        n_r <- n_r+1
        j <- n_new
        n_new <- n_new-1
        x.data <- x.data[,-i]
        x.lon <- x.lon[-i]
        x.lat <- x.lat[-i]
        x.name <- x.name[-i]
        x.height <- x.height[-i]
        i <- i-1
      }
      else {
        print(paste(x.lon[i],x.lat[i],sum(!is.na(x.data[,i])),"values: kept"))
        print(paste(x.lon[j],x.lat[j],sum(!is.na(x.data[,j])),"values: removed"))
        n_r <- n_r+1
        n_new <- n_new-1
        x.data <- x.data[,-j]
        x.lon <- x.lon[-j]
        x.lat <- x.lat[-j]
        x.name <- x.name[-j]
        x.height <- x.height[-j]
      }
    }
    else {
      j <- j+1
    }
  }
  i <- i+1
}
print(paste(n_r,"series removed"))


input <- list(data=x.data, lon=x.lon, lat=x.lat, name=x.name, height=x.height, time=x.time)
