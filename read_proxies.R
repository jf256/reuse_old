extract_lonlat <- function(x){
    if (!is.na(x)){
        xout <- gsub('deg.*', '', x)
        if (length(grep('-', xout)) == 1){
            xtmp <- strsplit(xout, '-')
            xout <- mean(as.numeric(unlist(xtmp)))
        }
        xout <- as.numeric(xout)
        if (length(grep('min', x)) == 1){
            xadd <- as.numeric(gsub('min.*', '', gsub('.*deg *', '', x)))
            xout <- xout + xadd/60
        }
    } else {    
        xout <- NA
    }
    xout
}        

compute_lonlat <- function(x){
    dir <- gsub('[-. a-z0-9]*', '', x)
    n.i <- dir[,1] %in% c('N', 'S')
    n.i <- cbind(n.i + 1, (!n.i) + 1)
    s.i <- (dir %in% c('E', 'N'))*2 - 1
    xout <- s.i * apply(x, 1:2, extract_lonlat)
    
    # swap x accordingly
    x2 <- cbind(xout[cbind(1:84, n.i[,1])], xout[cbind(1:84, n.i[,2])])
    # remove lons > 180
    lon.i <- x2[,1] > 180 & !is.na(x2[,1])
    x2[lon.i,1] <- x2[lon.i,1] - 360
    xout <- data.frame(lon=x2[,1], lat=x2[,2])    
    return(xout)
}

extra_lonlat <- function(x){
    if (all(!is.na(x))){
        xout <- gsub('deg.*', '', x)
        if (length(grep('min', x)) > 0){
            xadd <- as.numeric(gsub('min.*', '', gsub('.*deg *', '', x)))
            xout <- as.numeric(xout) + xadd/60
        }
        if (length(grep('-', xout)) > 0){
            xtmp <- strsplit(xout, '-')
            if (is.na(xtmp[[1]][2])) { xtmp[[1]][2]=xtmp[[1]][1] }
            if (is.na(xtmp[[2]][2])) { xtmp[[2]][2]=xtmp[[2]][1] }
            xout <- rbind(xtmp[[1]][c(1,1,2,2,1)], xtmp[[2]][c(1,2,2,1,1)])
        }
        if (is.matrix(xout)) {
            xout <- array(as.numeric(xout), dim(xout))
        } else {
            xout <- as.numeric(xout)
        }
        dir <- gsub(' ', '', gsub('[-. a-z0-9]*', '', x))
        n.i <- dir %in% c('N', 'S') + 1
        s.i <- (dir %in% c('E', 'N'))*2 - 1
        xout <- s.i * xout
        if (is.matrix(xout)){
            xout <- xout[n.i,]
        } else {
            xout <- xout[n.i]
        }
    } else {    
        xout <- c(NA,NA)
    }
    xout
}        
        
tmp <- read.table('~/unibe/data/proxies/temperature_proxies.csv', sep=',', header=TRUE, na.string=c(' -', '-'))

lon <- as.character(tmp[[grep('Longitude', names(tmp))]])
lat <- as.character(tmp[[grep('Latitude', names(tmp))]])

lola <- list()
for (i in seq(along=lon)) lola[[i]] <- c(lon[i], lat[i])

proxy <- lapply(lola, extra_lonlat)
#proxy[[45]][1,3:4] <- 40
proxies <- compute_lonlat(cbind(lon, lat))
# remove duplicates
proxies <- proxies[!duplicated(proxies) & !is.na(proxies[,1]),]


#map(interior=F)
#points(proxies, pch=19, col=2)

#lapply(proxy, function(x) {
#    if (is.matrix(x)){
#        polygon(x[1,], x[2,], lwd=2)
#    } else {
#        points(x[1], x[2], col=4, pch=19)
#    }})

