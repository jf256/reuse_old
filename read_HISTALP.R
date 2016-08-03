### read_HISTALP.R --- 
## 
## Filename: read_HISTALP.R
## Description: 
## Author: Joerg Franke
## Maintainer: 
## Created: 25 May 2011
## Version: 
## Last-Updated: 
##           By: Joerg
##     Update #: 20
## URL: 
## Keywords: 
## Compatibility: 
## 
######################################################################
## 
### Commentary: 
## 
##  reads HistALP data
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




years <- 1800:1849

stationfile <- '/Users/joerg/data/unibe/data/instr/histalp/histalp_temp.csv'





tmp=read.table(stationfile,header=T)
colnames(tmp)=c("id","k_nam","name","l_code","lon","lat","hoehe","datum","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec","spr","sum","aut","win","apr-sep","oct-mar","ann")

stat.i <- sort(unique(tmp[,1]))
tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))
lon.arr <- array(NA,length(stat.i))
lat.arr <- array(NA,length(stat.i))
elev.arr <- array(NA,length(stat.i))
name.arr <- array(NA,length(stat.i))

for (i in 1:length(seq(stat.i))){
   s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
   y.i <- which(years %in% tmp[s.i,'datum'])
   r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,'datum'] == years[x])))
   tmp.arr[y.i,,i] <- as.matrix(tmp[s.i[r.i],seq(9,20)])
   lon.arr[i] <- tmp[s.i[1],'lon']
   lat.arr[i] <- tmp[s.i[1],'lat']
   elev.arr[i] <- tmp[s.i[1],'hoehe']
   name.arr[i] <- as.character(tmp[s.i[1],'name'])
 }

# check if at least 80% of data is available and NOT missing (NA)
mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])

histalp.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
histalp.data=array(as.numeric(histalp.data.tmp),c(dim(histalp.data.tmp)[1],dim(histalp.data.tmp)[2]))

histalp <- list(data=histalp.data/10, lon=lon.arr[mask], lat=lat.arr[mask], name=name.arr[mask], height=elev.arr[mask], time=seq(min(years) + 1/24, by=1/12, length=nrow(histalp.data)))

save(histalp, file='histalp.Rdata')

#q(save='no')


######################################################################
### read_HISTALP.R ends here
