### read_GHCN_v3.R --- 
## 
## Filename: read_GHCN_v3.R
## Description: 
## Author: Jonas Bhend
## Maintainer: 
## Created: Wed Aug  4 16:52:19 2010 (+0200)
## Version: 
## Last-Updated: Tue May 10 2011
##           By: Jörg
##     Update #: 20
## URL: 
## Keywords: 
## Compatibility: 
## 
######################################################################
## 
### Commentary: 
## 
##  reads GHCN v3 temperature data and provides station information
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

setwd('/Users/joerg/data/unibe/projects/EnSRF/r/')




years <- 1850:1899





stationfile <- '/Users/joerg/data/unibe/data/instr/ghcn/ghcn_mm_t2m_adj.inv'

ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
stations <- read.fwf(stationfile, ww, header=F, fill=T, na.string='-999')
colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
stations[stations == -9999 ] <- NA

ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) #5, 1, 1, 1, 82, 5, 1, 1, 1)
tmp <- as.matrix(read.fwf('/Users/joerg/data/unibe/data/instr/ghcn/ghcn_mm_t2m_adj.dat', ww2, header=F, fill=T, na.string='-9999'))
colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')

# ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README

stat.i <- sort(unique(tmp[,1]))
tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))

for (i in 1:length(seq(stat.i))){
   s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
   y.i <- which(years %in% tmp[s.i,2])
   r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
   tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
 }

stations.new <- stations[which(stations[,'ID'] %in% stat.i),]

# check if at least 80% of data is available and NOT missing (NA)
mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])

ghcn.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
ghcn.data=array(as.numeric(ghcn.data.tmp),c(dim(ghcn.data.tmp)[1],dim(ghcn.data.tmp)[2]))

ghcn <- list(data=ghcn.data/100, lon=stations.new[mask,'LONGITUDE'], lat=stations.new[mask,'LATITUDE'], name=gsub(' *$', '', as.character(stations[mask,'NAME'])), height=stations.new[mask,'STNELEV'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))

save(ghcn, stations, file='ghcn.Rdata')

#q(save='no')


######################################################################
### read_GHCN_v3.R ends here
