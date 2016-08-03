year_min <- 1749
year_max <- 2012
n_stations <- 132

x.data <- array(dim=c((year_max-year_min+1)*12,n_stations))
x.lon <- array(dim=c(n_stations))
x.lat <- array(dim=c(n_stations))
x.name <- array(data="slp",dim=c(n_stations))
x.height <- array(dim=c(n_stations))
x.time <- seq(year_min + 1/24, by=1/12, length=nrow(x.data))
i_staz <- 0


### READ KUETTEL
coords <- read.table("Stations_PP_1722_coordinates.csv",skip=1)
data <- read.table("Stations_PP_1722_monthly_final.csv",skip=1)
for (i_staz in 1:dim(coords)[1]) {
  x.lon[i_staz] <- coords[i_staz,3]
  x.lat[i_staz] <- coords[i_staz,2]
  i_data <- 12
  for (i in 1:dim(data)[1]) {
    x.data[i_data,i_staz] <- data[i,i_staz+2]
    i_data <- i_data+1
  }
}


### READ HISTALP
#files <- read.table("list1")
#for (filename in t(files)) {
#  i_staz <- i_staz+1
#  header <- readLines(filename,n=3)
#  x.lon[i_staz] <- as.numeric(substr(header[3],36,40))
#  x.lat[i_staz] <- as.numeric(substr(header[3],42,46))
#  x.height[i_staz] <- as.integer(substr(header[3],49,52))
#  data <- read.table(filename,na.strings="99999",skip=4)
#  i_data <- (data[1,1]-year_min)*12+1
#  for (i in 1:dim(data)[1]) {
#    for (j in 2:13) {
#      x.data[i_data,i_staz] <- data[i,j]/10
#      i_data <- i_data+1
#    }
#  }
#}


### READ CRU STATIONS
files <- read.table("list2")
for (filename in t(files)) {
  i_staz <- i_staz+1
  header <- readLines(filename,n=2)
  x.lon[i_staz] <- as.numeric(substr(header[2],1,6))
  x.lat[i_staz] <- as.numeric(substr(header[2],8,12))
  data <- read.table(filename,na.strings="-999",skip=2)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]/10
      i_data <- i_data+1
    }
  }
}


### READ GHCN
files <- read.table("list3")
for (filename in t(files)) {
  i_staz <- i_staz+1
  header <- readLines(filename,n=5)
  x.lon[i_staz] <- as.numeric(substr(header[3],25,31))
  x.lat[i_staz] <- as.numeric(substr(header[3],16,21))
  x.height[i_staz] <- as.integer(substr(header[3],35,38))
  data <- read.table(filename,na.strings="-999.9",skip=5)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]
      i_data <- i_data+1
    }
  }
}


### READ OTHERS
files <- read.table("list4")
for (filename in t(files)) {
  i_staz <- i_staz+1
  header <- readLines(filename,n=1)
  x.lon[i_staz] <- as.numeric(substr(header,1,7))
  x.lat[i_staz] <- as.numeric(substr(header,9,14))
#  x.height[i_staz] <- as.integer(substr(header,16,19))
  data <- read.table(filename,na.strings="-999.9",skip=1)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]
      i_data <- i_data+1
    }
  }
}


slp <- list(data=x.data, lon=x.lon, lat=x.lat, name=x.name, height=x.height, time=x.time)
save(slp,file="../slp.Rdata")
write.table(cbind(x.lon,x.lat),file="list_coord",sep="\t",col.names=FALSE,row.names=FALSE)
#system("./plot_stations.gmt")
