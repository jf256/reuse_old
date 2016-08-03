### read_GHCN.R --- 
## 
## Filename: read_GHCN.R
## Description: 
## Author: Jonas Bhend
## Maintainer: 
## Created: Wed Aug  4 16:52:19 2010 (+0200)
## Version: 
## Last-Updated: Mon Dec  6 17:32:13 2010 (+1100)
##           By: Jonas
##     Update #: 19
## URL: 
## Keywords: 
## Compatibility: 
## 
######################################################################
## 
### Commentary: 
## 
##  reads GHCN temperature data and provides station information
## 
######################################################################
## 
### Change Log:
## 
## 
######################################################################
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, write to
## the Free Software Foundation, Inc., 51 Franklin Street, Fifth
## Floor, Boston, MA 02110-1301, USA.
## 
######################################################################
## 
### Code:

setwd('/Users/joerg/data/unibe/data/instr/ghcn/')




years <- 1750:1849





stationfile <- 'ghcn_mm_t2m_adj.inv'

#ww <- c(11, -1, 19, -1, 9, -1, 6, -1, 7, -1, 5, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 4, -1, 10, -1, 1)
ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
stations <- read.fwf(stationfile, ww, header=F, fill=T)
#colnames(stations) <- c('number', 'name', 'country', 'lat', 'lon', 'height', 'height2', 'pop',  'npop', 'topo', 'vegetation', 'stloc', 'coastdist', 'airport', 'itowndist', 'gridveg', 'flag')
colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
stations[stations == -9999 ] <- NA

#        ID                 1-11        Integer
#        LATITUDE          13-20        Real
#        LONGITUDE         22-30        Real
#        STNELEV           32-37        Real
#        NAME              39-68        Character
#        GRELEV            70-73        Integer
#        POPCLS            74-74        Character
#        POPSIZ            76-79        Integer
#        TOPO              80-81        Character
#        STVEG             82-83        Character
#        STLOC             84-85        Character
#        OCNDIS            86-87        Integer
#        AIRSTN            88-88        Character
#        TOWNDIS           89-90        Integer
#        GRVEG             91-106       Character
#        POPCSS            107-107      Character

#ww2 <- c(11, 1, 4, rep(5,12))
ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) #5, 1, 1, 1, 82, 5, 1, 1, 1)
tmp <- as.matrix(read.fwf('ghcn_mm_t2m_adj.dat', ww2, header=F, fill=T, na.string='-9999'))
colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')

# ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README

#           ID                 1-11        Integer
#           YEAR              12-15        Integer
#           ELEMENT           16-19        Character
#           VALUE1            20-24        Integer
#           DMFLAG1           25-25        Character
#           QCFLAG1           26-26        Character
#           DSFLAG1           27-27        Character
#             .                 .             .
#             .                 .             .
#             .                 .             .
#           VALUE12          108-112       Integer
#           DMFLAG12         113-113       Character
#           QCFLAG12         114-114       Character
#           DSFLAG12         115-115       Character

# years <- 1750:1850 # now defined at the top
stat.i <- sort(unique(tmp[,1]))

tmp.arr <- array(NA,c(length(years), 12, length(stat.i)),dimnames=c('yr','mon','stat'))
#c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))


# remove gaps in time series by filling with NA
# search for longest contiguous piece with na.contiguous
for (i in 1:length(seq(stat.i))){
   s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
   y.i <- which(years %in% tmp[s.i,2])
   r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
   tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
 }

# original
#    s.i <- which(tmp[,1] == stat.i[i] & tmp[,2] == 0)
#    y.i <- which(years %in% tmp[s.i,3])
#    r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,3] == years[x])))
#    tmp.arr[y.i,,i] <- tmp[s.i[r.i], 4:15]
#  }

#stations <- stations[which(stations[,'number'] %in% stat.i),]
stations.new <- stations[which(stations[,'ID'] %in% stat.i),]

mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])

ghcn.data <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
ghcn <- list(data=t(as.numeric(ghcn.data))/10, lon=stations.new[mask,'lon'], lat=stations.new[mask,'lat'], name=gsub(' *$', '', as.character(stations[mask,'name'])), height=stations.new[mask,'height'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))

save(ghcn, stations, file='ghcn.Rdata')

q(save='no')


######################################################################
### read_GHCN.R ends here
