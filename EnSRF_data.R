#FIX:
# strange no of input data 1641-1643  
# [37,] 1639     0   30   164   35  129
# [38,] 1640     0   30   164   35  129
# [39,] 1641  1668   30    50   19   31
# [40,] 1642  1668   30    50   19   31
# [41,] 1643  1668   30    50   19   31
# [42,] 1644     0   30   164   35  129
# [43,] 1645     0   30   164   35  129


# current code in R only without Fortran!!!
# ATTENTION: to use the faster fortran EnSRF routine execute on the command line:
# R CMD SHLIB EnSRF.f95 uncomment/comment related code in EnSRF_functions.R
# R --arch=x86_64 CMD SHLIB EnSRF.f95 # to compile EnSRF.f95 on 64bit MAC

rm(list=ls())
setwd('~/unibe/projects/EnSRF/r/src')
source('EnSRF_functions.R')
#source('r_functions_joerg.r')

# switches and options
machine="macbook" #"climcal3" # "climpa12" # 
if (machine=="macbook") {
  datadir="/Volumes/DATA/climdata/"
} else {
  datadir="/scratch/joerg/climdata/"
}

syr=1600 #1602               # startyear
  syr2=syr+1
eyr=2005 #1620               # endyear

#!!!ATTENTION: sixmonstatevector year starts in October of previous year (cyr-1)
sixmonstatevector=T    # 6 months of data in state vector for real proxy multiple 
                       # regression approach. ATTENTION: monthly has to be TRUE
if (sixmonstatevector) {
  s <- 2 # only 2 seasons to calculate but still monthly results in long state vector
} else {
  s <- 12
}
monthly_out=F          # if sixmonstatevector=T output is backtransformed to seasonal 
                       # average or monthly data if monthly_out=T 
write_coor=T           # write ascii files with assimilated stations and data per ts
season=c(3,9)          # 3,9 = apr-sep and oct-mar, num=end month of season
                       # season=c(2,5,8,11)
# distances estimated from echam decorrelation excercise
calc_decorr_dist=F     # calculate decorrelation distance for each variable to set L
# stefans recommendation: increase distance, give proxies more weight and instr. less weight
l_dist_temp2=1000*1.5 # factor *1.5 after stefans recommendation
l_dist_slp=1800*1.5
l_dist_precip=300*1.5
l_dist_gph500=1800*1.5
l_dist_gph100=2500*1.5
l_dist_u850=1200*1.5
l_dist_u200=1200*1.5
l_dist_v850=1000*1.5
l_dist_v200=1000*1.5
l_dist_omega500=300*1.5
l_dist_t500=1000*1.5
l_dist_ind=999999

nmem=30
first_prox_per_grid=F  # first proxy per echam grid box ATTENTION: only this 
                       # or second next option (avg_prox_per_grid) can be TRUE
  firstproxres=10      # grid resolution for instr. stations (5 = echamgrid/5)
avg_prox_per_grid=T    # average more than one proxy per echam grid box 
                       # and calc proxy vs echam correlation
reduced_proxies=F      # use every ??th (see code below) proxy record
every2grid=T           # only use every third grid cell of ECHAM, CRU validation, ...
land_only=F            # calc on land only
fasttest=F             # use even less data
tps_only=F             # only use temp, precip and slp in state vector, remove other vars
 no_stream=T           # all echam vars but stream function as there is problem with 5/9 levels
                       # which are in lat dimension
ind_ECHAM=T            # add ECHAM indices to state vector
load_71yr_anom=T       # load 71yr echam anomalies calculated with cdo
anom_reload=F          # anom first calculated in R
anom_save=F            # save anom calculated in R to reload next time
if (load_71yr_anom==T) {
  anom_reload=F
  anom_save=F}
check_assimdata=T      # screen assimilation data before using it

if (no_stream & tps_only) {
  tps_only = F
  print('ACHTUNG: tps_only was set to FALSE')
}

# load or generate data from scratch
generate_ECHAM_1901_70=F # ECHAM data for bias calc with petra's real_prox trw data 
generate_ECHAM=F       # if TRUE -> orig. echam data is read
                       # if FALSE -> echam_syr-eyr.Rdata
generate_ECHAM_103=F   # ECHAM ens. mem. 103 with corrected land use forcing
generate_ECHAM_anom=F  # read echam anom, clim and sd from cdo
# ATTENTION if echam=T the proxy/instr. data have to be generated, too
generate_NCEP=F        # generate NCEP/NCAR reanalysis for independent verification
# next line not included yet: 
generate_20CR=F        # generate 20CR reanalysis for independent verification
#if ((syr > 1900) & (eyr < 2005)) {
# generate_CRUALLVAR=T  # if TRUE -> read CRUTEM3 temp, TS3 precip and 
#} else {               # HADSLP2 gridded instrumentals
generate_CRUALLVAR=F   # if FALSE -> cru_allvar.Rdata 
generate_HadCRU4=F     # HadCRU ens. SD for instr. uncertainty and error-spread ratio
#}
if ((syr < 1901) & (eyr > 1749)) {
 generate_LUTPAULKUT=T # if TRUE -> read Luterbacher, Pauling, 
} else {               # Kuettel's temp. precip. and SLP 
 generate_LUTPAULKUT=F # gridded seasonal recons (1750-1999)
} 
generate_ind_ECHAM=T   # gen echam indices
if (syr > 1900) {
 generate_ind_recon=T   # gen recon indices
} else {
 generate_ind_recon=F
} 
# use scripts in data_yuri to generate .Rdata files 
generate_t_yuri=F      # if TRUE -> yuri's temp. data collection is read
generate_slp_yuri=F    # if TRUE -> yuri's slp data collection is read
generate_GHCN=F        # if TRUE -> orig. GHCN stat. data is read; 
generate_GHCN_precip=F # if FALSE -> ghcn.Rdata 
#generate_HISTALP=F     # histalp temp. is already in yuri's collection 
if (syr < 1960) {
 generate_PROXIES=T     # if FALSE -> real_proxies.Rdata
} else {
 generate_PROXIES=F
} 
  trw_only=F           # Petra's TRW only
  mxd_only=F           # Use only MXD tree ring proxies, NOT Petra's TRW
  schweingr_only=F     # Use Schweingruber MXD grid only
generate_DOCU=F        # if TRUE -> yuri's docu. data collection is read 

# choose data (only 1 of the following options must be TRUE!!!)
if (eyr > 1659) {
 instrumental=T        # all instrumental stations
} else {
 instrumental=F
}
  yuri_temp=T          # yuri's data compilation, SLP always loaded
  yuri_slp=T
  ghcn_temp=T
  ghcn_prec=F
#inst_at_proxy=F        # only few instrumental stations close to proxy sites
#inst_at_proxy_noise=F  # only few instrumental stations close to proxy sites + noise
if (syr < 1960) {
 real_proxies=T         # Proxy data experiment (regression NOT H operator)
} else {
 real_proxies=F
}
if (syr > 1853) {
 docu=F                 # read documentary based data
} else {
 docu=T
}
print(paste("instr:",instrumental, "proxies:",real_proxies, "documentary:",docu))

# other options
scaleprox=T            # scale standardized docu and prox data the echam variance at location
anomaly_assim=T        # work with anomalies to avoid reg. konst in state vector
 shortanom=F # all just 1950-70 because no other data on 2nd_grid to test
 nseas <- 12 # year with 12 months
#subbias=F              # TO BE USED WITH INSTRUMENTAL DATA !!!
                       # substract bias from instrumental data to get close to model data
                       # to avoid problems with correlation between vars of different units
subbias_vali=T         # also debias cru/ncep validation data 
check_dist=F           # test for ideal cut-off distance of spatial correlations
#H_non_lin=F           # new H operator that also allows non-linear functions

# choose validation data set
# ONLY one can be TRUE
# next line not included yet: 
if (eyr < 1750) {
 vali=F                 # switch off prepplot if no vali data selected
} else {
 vali=T
}
twcr_vali=F            # 20CR reanalysis data for validation
ncep_vali=F            # NCEP/NCAR reanalysis data for validation
if ((syr > 1900) & (eyr<2006)) {
  cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
  ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU and NCEP
} else {
  cru_vali=F 
  ind_recon=F
}
if ((syr < 1901) & (eyr > 1749)) {
 recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
} else {
 recon_vali=F
}
# if (eyr<1750) {
#   generate_LUTPAULKUT=F
#   recon_vali=F
#   vali=F
#   print("SET LUTPAULKUT AND VALI TO FALSE BECAUSE EYR < 1750!!!")
# }

# -----------------------------------------------------------------------------------------





##########################################################################################
# start loading data  
##########################################################################################  
if (generate_ECHAM_1901_70){
  print("generate_ECHAM_1901_70")
  # add vars to state vector to calculate clim indices
  # always for every grid cell as just used for proxy calibration
  echam1901_70 <- read_echam_ensmean('EnSRF', timlim=c(1901,1970),small=F)
  save(echam1901_70, file="../data/echam_1911-70.Rdata")
} 

if (generate_ECHAM_103){
  print("generate_ECHAM ens. mem. 103")
    read_echam4(filehead='EnSRF', path=paste0(datadir,'echam/echam103'), timlim=c(syr,eyr), 
                small=every2grid, landonly=land_only)
}  


if (generate_ECHAM){
  print("generate_ECHAM")
  read_echam4('EnSRF', timlim=c(syr,eyr), small=every2grid, landonly=land_only)
} 


if (generate_ECHAM_anom){
  # read echam 70yr anom, clim and sd
  print("generate_ECHAM")
  read_echam4('ano', path=echanompath, timlim=c(1601,2005), small=every2grid, landonly=land_only, anom=T)
  read_echam4('EnSRF', path=echclimpath, timlim=c(1635,1970), small=every2grid, landonly=land_only, clim=T)
  read_echam4('EnSRF', path=echsdpath, timlim=c(1601,2005), small=every2grid, landonly=land_only, std=T)
#
#
}
  
  
if (generate_ind_ECHAM) {
  print("generate_ind_ECHAM")
  if ((shortanom==F) & (anomaly_assim==T)) {
    echind <- read_echam1('EnSRF',timlim=c((syr-35),(eyr+35)), path=echindpath,small=F) # small has to be FALSE
    save(echind, file=paste("../data/indices_echam_",syr,"-",eyr,".Rdata",sep=""))
  } else {
    echind <- read_echam1('EnSRF',timlim=c(syr,eyr), path=echindpath,small=F) # small has to be FALSE
    save(echind, file=paste("../data/indices_echam_",syr,"-",eyr,".Rdata",sep=""))
  }
}

if (generate_ind_recon){
  if (syr<1901){syr_ind=1901} else {syr_ind=syr}
  if (eyr>2004){eyr_ind=2004} else {eyr_ind=eyr}
  ind=read.table(file=paste(datadir,'/indices/stefan/stefan_monthly_indices.txt'
                            ,sep=''),header=T)
  # HC HCL SJ SJLAND Z300 Z100 DIMI NAO PNA PWC
  ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
                           start=ind[1,colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  ind_rec_z100 = window(ts(ind[,colnames(ind) == 'Z100'],start=ind[1, 
                           colnames(ind) == 'yr'],freq=12),syr_ind,freq=12,c(eyr_ind,12))
  ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
                           colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
                           colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
                           colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,
                           colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc, 
                  ind_rec_sj))
  save(indall, file=paste("../data/indices_recon_",syr,"-",eyr,".Rdata",sep=""))
}

if (generate_NCEP){
  print("generate_NCEP")
  if ((!eyr<1948) & (!syr>2009)){
    if (syr<1948){syr_ncep=1948} else {syr_ncep=syr}
    if (eyr>2009){eyr_ncep=2009} else {eyr_ncep=eyr}
    #  see script in EnSRF/script/merge_ncep.sh for regridding and co of orig. ncep data set
    ncepall <- read_echam1('ncep_allvar_1948-2009',timlim=c(syr_ncep,eyr_ncep),
                            path=nceppath,small=every2grid)
#     ncepall0 <- read_echam1('ncep_allvar_1948-2009',timlim=c(syr_ncep,eyr_ncep),
#                            path=nceppath,small=every2grid)
#     if (ind_recon){
#       ind=read.table(file=paste(datadir,'/indices/stefan/stefan_monthly_indices.txt'
#                                 ,sep=''),header=T)
#       # HC HCL SJ SJLAND Z300 Z100 DIMI NAO PNA PWC
#       ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
#                       start=ind[1,colnames(ind) == 'yr'],freq=12),syr_ncep,c(eyr_ncep,12))
#       ind_rec_z100 = window(ts(ind[,colnames(ind) == 'Z100'],start=ind[1, 
#                       colnames(ind) == 'yr'],freq=12),syr_ncep,freq=12,c(eyr_ncep,12))
#       ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_ncep,c(eyr_ncep,12))
#       ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_ncep,c(eyr_ncep,12))
#       ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
#                      colnames(ind) == 'yr'],freq=12),syr_ncep,c(eyr_ncep,12))
#       ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,
#                      colnames(ind) == 'yr'],freq=12),syr_ncep,c(eyr_ncep,12))
#       indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc, 
#                        ind_rec_sj))
#       ncepall <- ncepall0
#       ncepall$data <- rbind(ncepall0$data[1:dim(ncepall0$data)[1],,],indall) #,ncepall0$data[2083:2087,,])
#       ncepall$names <- c(ncepall0$names[1:dim(ncepall0$data)[1]],rownames(indall)) 
#       #,ncepall0$names[2083:2087])
#       ncepall$lon <- c(ncepall$lon,rep(NA,nrow(indall))) 
#       ncepall$lat <- c(ncepall$lat,rep(NA,nrow(indall)))
#     } else {
#      ncepall <- ncepall0
#    }
    if (subbias_vali) {
      load(paste("../data/echam_",syr,"-",eyr,".Rdata",sep=""))
      # cut temp precip slp
      e_tps <- echam
      e_tps$data <- echam$ensmean[((echam$name=="temp2")|(echam$name=="precip")|
                                     (echam$name=="slp")),]
      e_tps$ensmean <- e_tps$data 
      c_tps <- ncepall
      c_tps$data <- ncepall$data[((ncepall$name=="temp2")|(ncepall$name=="precip")|
                                   (ncepall$name=="slp")),]
      c_tps$ensmean <- c_tps$data 
      # split time dim in months, calc monthly mean over time,  both for echam and cruall
      e_c_bias <- bias_fun(c_tps,e_tps,seas=12)
      # split time dim in months and substract echam bias
      tmp <- array(0,c(dim(ncepall$data)[1]-dim(e_c_bias$data)[1],dim(e_c_bias$data)[2:3]))
      e_c_bias$data <- abind(e_c_bias$data,tmp,along=1) #[2082:2087,1,1]
      ncep_tmp <- array(ncepall$data, c(nrow(ncepall$data), 12, 
                                      ncol(ncepall$data)/12))
      cdim <- dim(ncepall$data)
      ncepall$data <-  array((ncep_tmp - as.vector(e_c_bias$data)),cdim)
      # convert 12 21 back to one dimension 252      
    }
    if (subbias_vali) {
      if (every2grid){
        save(ncepall, file=paste("../data/ncep_allvar_debias_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
      } else {
        save(ncepall, file=paste("../data/ncep_allvar_debias_",syr,"-",eyr,".Rdata",sep=""))
      }
    } else {
      if (every2grid) {
        save(ncepall, file=paste("../data/ncep_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
      } else {
        save(ncepall, file=paste("../data/ncep_allvar_",syr,"-",eyr,".Rdata",sep=""))
      }
    }
  } 
}

if (generate_CRUALLVAR){
  print("generate_CRUALLVAR")
  if ((!eyr<1901) & (!syr>2004)){
#    if (syr<1901){syr_cru=1901} else {syr_cru=syr}
    syr_cru=1901
#    if (eyr>2004){eyr_cru=2004} else {eyr_cru=eyr}
    eyr_cru=2004
#  see script in EnSRF/script/merge_cru.sh for regridding and co of orig. cru data set
    cruall <- read_echam1('cru_allvar_abs_1901-2004.nc',timlim=c(syr_cru,eyr_cru),
                       path=crupath,small=every2grid,landonly=land_only)
#     cruall0 <- read_echam1('cru_allvar_abs_1901-2004.nc',timlim=c(syr_cru,eyr_cru),
#                            path=crupath,small=every2grid,landonly=land_only)
#     if (ind_recon){
#       ind=read.table(file=paste(datadir,'/indices/stefan/stefan_monthly_indices.txt'
#                                   ,sep=''),header=T)
#     # HC HCL SJ SJLAND Z300 Z100 DIMI NAO PNA PWC
#       ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
#                       start=ind[1,colnames(ind) == 'yr'],freq=12),syr_cru,c(eyr_cru,12))
#       ind_rec_z100 = window(ts(ind[,colnames(ind) == 'Z100'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_cru,freq=12,c(eyr_cru,12))
#       ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_cru,c(eyr_cru,12))
#       ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_cru,c(eyr_cru,12))
#       ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
#                       colnames(ind) == 'yr'],freq=12),syr_cru,c(eyr_cru,12))
#       ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,colnames(ind) 
#                    == 'yr'],freq=12),syr_cru,c(eyr_cru,12))
#       indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc, 
#                        ind_rec_sj))
#       cruall <- cruall0
#       cruall$data <- rbind(cruall0$data[1:dim(cruall0$data)[1],,],indall) #,cruall0$data[2083:2087,,])
#       cruall$names <- c(cruall0$names[1:dim(cruall0$data)[1]],rownames(indall)) 
#       cruall$lon <- c(cruall$lon,rep(NA,nrow(indall))) 
#       cruall$lat <- c(cruall$lat,rep(NA,nrow(indall))) 
#       #,cruall0$names[2083:2087])
#     } else {
#       cruall <- cruall0
#     }
    if (subbias_vali) {
     if (load_71yr_anom) {
       cruall2 <- read_echam1('cru_allvar_abs_1901-2004.nc',timlim=c(1917,1987),
                             path=crupath,small=every2grid,landonly=land_only)
       if (every2grid) {
         # load 70 year climatology from ECHAM to calc bias
         load("../data/echam_clim/echam_clim_1952-1953_2ndgrid.Rdata")
       } else {
         load("../data/echam_clim/echam_clim_1952-1953.Rdata")
       }
       echam_clim$data[echam_clim$names=='temp2',,] <- echam_clim$data[echam_clim$names=='temp2',,] - 273.15
       echam_clim$ensmean[echam_clim$names=='temp2',] <- echam_clim$ensmean[echam_clim$names=='temp2',] - 273.15
       echam_clim$data[echam_clim$names=='precip',,] <- echam_clim$data[echam_clim$names=='precip',,] * 3600 * 24 * 30
       echam_clim$ensmean[echam_clim$names=='precip',] <- echam_clim$ensmean[echam_clim$names=='precip',] * 3600 * 24 * 30
       echam_clim$data[echam_clim$names=='slp',,] <- echam_clim$data[echam_clim$names=='slp',,] / 100
       echam_clim$ensmean[echam_clim$names=='slp',] <- echam_clim$ensmean[echam_clim$names=='slp',] / 100
       e_tps <- echam_clim
       e_tps$data <- echam_clim$ensmean[((echam_clim$name=="temp2")|(echam_clim$name=="precip")|
                        (echam_clim$name=="slp")),1:12]
       e_tps$ensmean <- e_tps$data 
       c_tps <- cruall2
       c_tps$data <- cruall2$data[((cruall$name=="temp2")|(cruall2$name=="precip")|
                        (cruall2$name=="slp")),,1]
       c_tps$ensmean <- c_tps$data
       # split time dim in months, calc monthly mean over time,  both for echam and cruall
       e_c_bias <- bias_fun(c_tps,e_tps,seas=12)
       # split time dim in months and substract echam bias
       tmp <- array(0,c(dim(cruall$data)[1]-dim(e_c_bias$data)[1],dim(e_c_bias$data)[2:3]))
       e_c_bias$data <- abind(e_c_bias$data,tmp,along=1) #[2082:2087,1,1]
       cru_tmp <- array(cruall$data, c(nrow(cruall$data), 12, 
                                       ncol(cruall$data)/12))
       cdim <- dim(cruall$data)
         #debias <- array((cru_tmp - as.vector(e_c_bias$data)),cdim)
       cruall$data <-  array((cru_tmp - as.vector(e_c_bias$data)),cdim)
     } else {
#      load(paste("../data/echam_",syr,"-",eyr,".Rdata",sep=""))
      for (cyr in syr2:eyr) {
        print(cyr)
        if (every2grid) {
          load(paste("../data/echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep="")) 
        } else {
          load(paste("../data/echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
        }
        # ACHTUNG just works if variables are in same order, not with sixmonstatevector
        tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), 
                    which(echam$names=='slp'))
        if (cyr == syr2) {
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
        } else {
          echam$data <- echam$data[tpspos,1:12,]
          echam$ensmean <- echam$ensmean[tpspos,1:12]
        }
        echam$names <- echam$names[tpspos]
        echam$lon <- echam$lon[tpspos]
        echam$lat <- echam$lat[tpspos]
        if (cyr == syr2) {
          echam.allts=echam
        } else {
          echam.allts$data=abind(echam.allts$data,echam$data,along=2)
          echam.allts$ensmean=cbind(echam.allts$ensmean,echam$ensmean)
          echam.allts$time=c(echam.allts$time,echam$time)
        }
      }      
      echam <- echam.allts
      rm(echam.allts)
      # cut temp precip slp
      e_tps <- echam
      e_tps$data <- echam$ensmean[((echam$name=="temp2")|(echam$name=="precip")|
                                   (echam$name=="slp")),]
      e_tps$ensmean <- e_tps$data 
      c_tps <- cruall
      c_tps$data <- cruall$data[((cruall$name=="temp2")|(cruall$name=="precip")|
                                 (cruall$name=="slp")),,1]
      c_tps$ensmean <- c_tps$data 
      # split time dim in months, calc monthly mean over time,  both for echam and cruall
      e_c_bias <- bias_fun(c_tps,e_tps,seas=12)
      # split time dim in months and substract echam bias
      tmp <- array(0,c(dim(cruall$data)[1]-dim(e_c_bias$data)[1],dim(e_c_bias$data)[2:3]))
      e_c_bias$data <- abind(e_c_bias$data,tmp,along=1) #[2082:2087,1,1]
      cru_tmp <- array(cruall$data, c(nrow(cruall$data), 12, 
                       ncol(cruall$data)/12))
      cdim <- dim(cruall$data)
      #debias <- array((cru_tmp - as.vector(e_c_bias$data)),cdim)
      cruall$data <-  array((cru_tmp - as.vector(e_c_bias$data)),cdim)
      # convert 12 21 back to one dimension 252      
     }
    }
    if (subbias_vali) {
      if (every2grid) {
        save(cruall, file=paste("../data/cru_allvar_debias_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata",sep=""))
      } else {
        save(cruall, file=paste("../data/cru_allvar_debias_",syr_cru,"-",eyr_cru,".Rdata",sep=""))
      }
    } else {
      if (every2grid) {
        save(cruall, file=paste("../data/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata",sep=""))
      } else {
        save(cruall, file=paste("../data/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata",sep=""))  
      }
    }
  }
}


if (generate_HadCRU4){ 
  print("generate_HadCRU4")
  #cru4_may_sep <- read_echam1('HadCRUT4_ens_sd_may-sep_timemean.nc',timlim=c(2001,2002),
  #                              path=cru4path,small=every2grid,landonly=F)
  #cru4_oct_apr <- read_echam1('HadCRUT4_ens_sd_oct-apr_timemean.nc',timlim=c(2001,2002),
  #                            path=cru4path,small=every2grid,landonly=F)
  cru4_may_sep <- read_echam1('HadCRUT4_ens_sd_may-sep_yrmean.nc',timlim=c(1901,2002),
                              path=cru4path,small=every2grid,landonly=F)
  cru4_oct_apr <- read_echam1('HadCRUT4_ens_sd_oct-apr_yrmean.nc',timlim=c(1901,2002),
                              path=cru4path,small=every2grid,landonly=F)
  cru4_may_sep$ensmean <- cru4_may_sep$data
  cru4_oct_apr$ensmean <- cru4_oct_apr$data
  if (every2grid) {
    save(cru4_may_sep,cru4_oct_apr, file="../data/cru4_ens_sd_2ndgrid.Rdata")
  } else {
    save(cru4_may_sep,cru4_oct_apr, file="../data/cru4_ens_sd.Rdata")  
  }
}

if (generate_LUTPAULKUT){ 
  print("generate_LUTPAULKUT")
# ACHTUNG seasonal resolution
  if (syr<1750){syr_recon=1750} else {syr_recon=syr}
  if (eyr>1999){eyr_recon=1999} else {eyr_recon=eyr}
  #  see script in EnSRF/script/merge_recon.sh for regridding and co of orig. recon data set
  reconall <- read_echam1('recon_allvar_1750-1999',xlim=c(-180,180), ylim=c(-90,90), 
                          timlim=c(syr_recon,eyr_recon),path=reconpath,small=every2grid,landonly=land_only)
  reconall$data[reconall$names=="precip",,]<-reconall$data[reconall$names=="precip",,]/3
  if (every2grid) {
    save(reconall, file=paste("../data/recon_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
  } else {
    save(reconall, file=paste("../data/recon_allvar_",syr,"-",eyr,".Rdata",sep=""))  
  }
}

if (generate_GHCN){
  print("generate_GHCN")
  ghcn <- read_ghcn_refyr(1600,2005,1600,1869)
#  ghcn <- read_ghcn_refyr(syr,eyr,1600,1869)
#  ghcn <- read_ghcn_refyr(syr,eyr,syr,eyr)
  ghcn$names <-rep('temp2',length(ghcn$names))
  save(ghcn, file=paste0("../data/ghcn_",syr,"-",eyr,".Rdata"))
}

if (generate_GHCN_precip){
  print("generate_GHCN_precip")
  ghcn_precip <- read_ghcn_refyr_precip(1600,2005,1600,1869)
#  ghcn_precip <- read_ghcn_refyr_precip(syr,eyr,1600,1869)
#  ghcn_precip <- read_ghcn_refyr_precip(syr,eyr,syr,eyr)
  ghcn_precip$data <- ghcn_precip$data / 10 # to make units echam conform
  ghcn_precip$names <-rep('precip',length(ghcn_precip$names))
  save(ghcn_precip, file=paste0("../data/ghcn_precip_",syr,"-",eyr,".Rdata")) 
}

if (generate_t_yuri){
  print("generate_t_yuri")
  source("../data_yuri/t_assimil/read_all.R")
}

if (generate_slp_yuri){
  print("generate_slp_yuri")
  source("../data_yuri/slp_assimil/read_all.R")
}

#if (generate_HISTALP){
#  histalp <- read_histalp(timlim=c(syr,eyr))
#  histalp$names <- rep('temp2',length(histalp$names))
#  save(histalp, file='../data/histalp.Rdata')
#} 

if (generate_DOCU){
  print("generate_DOCU")
  source("../data_yuri/t_docu/read_seas.R")
  source("../data_yuri/t_docu/read_monthly.R")
  source("../data_yuri/t_docu/read_JFMA.R")
  source("../data_yuri/t_docu/read_AMJJA.R")
}

if (generate_PROXIES){
  print("generate_PROXIES")
  # real trw proxy multiple regression approach
  # only with monthly state vector of 6 months
  if (trw_only) {
     realprox <- read_proxy2(syr,eyr)
  } else if (schweingr_only) {
    realprox <- read_proxy_schweingr(syr,eyr)
  } else if (mxd_only) {
    realprox <- read_proxy_mxd(syr,eyr)
  } else {
    schprox <- read_proxy_schweingr(syr,eyr)
    mxdprox <- read_proxy_mxd(syr,eyr)
    trwprox <- read_proxy2(syr,eyr)
    realprox <- list()
    realprox$data <- cbind(mxdprox$data, schprox$data, trwprox$data)
    realprox$lon <- c(mxdprox$lon, schprox$lon, trwprox$lon)
    realprox$lat <- c(mxdprox$lat, schprox$lat, trwprox$lat)
    realprox$time <- mxdprox$time
    realprox$mr <- rbind(mxdprox$mr, schprox$mr, trwprox$mr)
    realprox$var_residu <- c(mxdprox$var_residu, schprox$var_residu, trwprox$var_residu)
  }
  save(realprox, file=paste0("../data/real_proxies_",syr,"-",eyr,".Rdata"))
} 

docu_backup <- docu






                                     

##########################################################################################
# loading echam timeslice (loop over years to reduce size of state vector)
##########################################################################################
if (sixmonstatevector) {syr2=syr+1} else {syr2=syr}
for (cyr in syr2:eyr) {
#for (cyr in syr2:1955) {  
  print(cyr)
  asyr <- cyr-35
  if (asyr < 1601) {asyr = 1601}
  if (shortanom) { asyr = syr }
  aeyr <- cyr+35
  if (aeyr > 2005) {aeyr = 2005}
  if (shortanom) { aeyr = eyr }
  ptm1 <- proc.time()
  if (every2grid) {
    load(paste("../data/echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  } else {
    load(paste("../data/echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
  }
#   if (every2grid) {
#     load(paste("../data/echam_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
#   } else {
#     load(paste("../data/echam_",syr,"-",eyr,".Rdata",sep=""))
#   }
  if (ind_ECHAM) {
    load(file=paste("../data/indices_echam_",syr,"-",eyr,".Rdata",sep=""))
    if ((shortanom==F) & (anomaly_assim==T)) {
      ti=which(floor(echind$time)>=asyr & floor(echind$time)<=aeyr)
      sts=ti[1]
      ets=ti[length(ti)]
      echind$data=echind$data[,sts:ets,]       
      echind$time=echind$time[sts:ets]
    }
  }
  if (subbias_vali) {
    if (every2grid) {
      if (ncep_vali) {load(paste("../data/ncep_allvar_debias__",syr,"-",eyr,"2ndgrid.Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
      } else if (cru_vali) {
#        load(paste("../data/cru_allvar_debias_",syr,"-",eyr,"_2ndgrid.Rdata",sep="")) # cru_vali
        load("../data/cru_allvar_debias_1901-2004_2ndgrid.Rdata") # cru_vali 
        ti=which(floor(cruall$time)>=syr & floor(cruall$time)<=eyr)
        sts=ti[1]
        ets=ti[length(ti)]
        cruall$data=cruall$data[,sts:ets,,drop=F]       
        cruall$time=cruall$time[sts:ets]
      }
    } else {
      if (ncep_vali) {load(paste("../data/ncep_allvar_debias_",syr,"-",eyr,".Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon_allvar_",syr,"-",eyr,".Rdata",sep=""))
      } else if (cru_vali) {load(paste("../data/cru_allvar_debias_",syr,"-",eyr,".Rdata",sep=""))} # cru_vali
    }
  } else {
    if (every2grid) {
      if (ncep_vali) {load(paste("../data/ncep_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))
      } else if (cru_vali) {load(paste("../data/cru_allvar_",syr,"-",eyr,"_2ndgrid.Rdata",sep=""))} # cru_vali
    } else {
      if (ncep_vali) {load(paste("../data/ncep_allvar_",syr,"-",eyr,".Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon_allvar_",syr,"-",eyr,".Rdata",sep=""))
      } else if (cru_vali) {load(paste("../data/cru_allvar_",syr,"-",eyr,".Rdata",sep=""))} # cru_vali
    }
  }
  if (ind_recon) {
    load(file=paste("../data/indices_recon_",syr,"-",eyr,".Rdata",sep=""))
  }
  if (cru_vali) {
    valiall <- cruall
  } else if (ncep_vali) {
    valiall <- ncepall
  } else if (recon_vali) {
    valiall <- reconall
    valiall$data <- valiall$data[,,1]
  } else { vali = F }
  if (instrumental){
    if (yuri_slp) {
      load('../data_yuri/slp.Rdata') # monthly slp collection from yuri
      inst_slp <- slp
      # cut time period as not done yet
      ti=which(floor(inst_slp$time)>=syr & floor(inst_slp$time)<=eyr)
      sts=ti[1]
      ets=ti[length(ti)]
      inst_slp$data=inst_slp$data[sts:ets,]       
      inst_slp$time=inst_slp$time[sts:ets]
    }
    if (yuri_temp) {
      load('../data_yuri/t.Rdata') # monthly temp collection from yuri
      inst_t <- t
      # cut time period as not done yet
      ti=which(floor(inst_t$time)>=syr & floor(inst_t$time)<=eyr)
      sts=ti[1]
      ets=ti[length(ti)]
      inst_t$data=inst_t$data[sts:ets,]       
      inst_t$time=inst_t$time[sts:ets]
    }
    if (ghcn_temp) {
#      load(paste0("../data/ghcn_",syr,"-",eyr,".Rdata"))
      load("../data/ghcn_1600-2005.Rdata")
      # cut time period as not done yet
      ti=which(floor(ghcn$time)>=syr & floor(ghcn$time)<=eyr)
      sts=ti[1]
      ets=ti[length(ti)]
      ghcn$data=ghcn$data[sts:ets,]       
      ghcn$time=ghcn$time[sts:ets]
    }
    if (ghcn_prec) {
#      load(paste0("../data/ghcn_precip_",syr,"-",eyr,".Rdata"))
      load("../data/ghcn_precip_1600-2005.Rdata")
      # cut time period as not done yet
      ti=which(floor(ghcn_precip$time)>=syr & floor(ghcn_precip$time)<=eyr)
      sts=ti[1]
      ets=ti[length(ti)]
      ghcn_precip$data=ghcn_precip$data[sts:ets,]       
      ghcn_precip$time=ghcn_precip$time[sts:ets]
    }
#    load('../data/histalp.Rdata')
  }
  
  if (real_proxies){
    load(paste0("../data/real_proxies_",syr,"-",eyr,".Rdata"))  
  }

  docu <- docu_backup
  if (docu){
#     load('../data_yuri/t_docu_seas.Rdata') 
#     doc_t_seas <- t
#     # cut time period as not done yet
#     ti=which(floor(doc_t_seas$time)>=syr & floor(doc_t_seas$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_seas$data=doc_t_seas$data[sts:ets,]
#     doc_t_seas$time=doc_t_seas$time[sts:ets]
    load('../data_yuri/t_docu_monthly.Rdata') 
    doc_t_mon <- t
    ti=which(floor(doc_t_mon$time)>=syr & floor(doc_t_mon$time)<=eyr)
    sts=ti[1]
    ets=ti[length(ti)]
    doc_t_mon$data=doc_t_mon$data[sts:ets,]
    doc_t_mon$time=doc_t_mon$time[sts:ets]
    if (!any(!is.na(doc_t_mon$data))) { docu=F }
#     load('../data_yuri/t_docu_JFMA.Rdata') 
#     doc_t_JFMA <- t
#     ti=which(floor(doc_t_JFMA$time)>=syr & floor(doc_t_JFMA$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_JFMA$data=doc_t_JFMA$data[sts:ets,]      
#     doc_t_JFMA$time=doc_t_JFMA$time[sts:ets]
#     load('../data_yuri/t_docu_AMJJA.Rdata') 
#     doc_t_AMJJA <- t
#     ti=which(floor(doc_t_AMJJA$time)>=syr & floor(doc_t_AMJJA$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_AMJJA$data=doc_t_AMJJA$data[sts:ets,]       
#     doc_t_AMJJA$time=doc_t_AMJJA$time[sts:ets]
  }

  print('calc time for loading data')
  print(proc.time() - ptm1)




# just leave temp precip slp in state vector  
  if (tps_only) {
    tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), 
      which(echam$names=='slp'), which(echam$names=='bias'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
    if (vali) {
      tpspos2 <- c(which(valiall$names=='temp2'), which(valiall$names=='precip'), 
                   which(valiall$names=='slp'))
      valiall$data <- valiall$data[tpspos2,]
      valiall$names <- valiall$names[tpspos2]
    }
  } else if (no_stream) {
    # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
    tpspos <- c(which(echam$names!='stream'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
  }
  if (fasttest) {
    mulc <- 4 # choose every 4th grid box
    loi <- seq(1:length(echam$lon))
    lai <- seq(1:length(echam$lat))
    di <- seq(1:dim(echam$data)[1])
    ni <- seq(1:length(echam$names))
    loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
    lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    di <- di[seq(ceiling(mulc/2), length(di),mulc)]
    ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
    echam$lon <- echam$lon[loi]
    echam$lat <- echam$lat[lai]
    echam$data <- echam$data[di,,]
    echam$ensmean <- echam$ensmean[di,]
    echam$names <- echam$names[ni]
#     lni1 <- (length(echam$names)-length(echind$names))
#     lni2 <- (length(echam$lon)-length(echind$names))# length no indices
#     pi1 <- seq(lni1+1,length(echam$names)) # pos indices
#     pi2 <- seq(lni2+1,length(echam$lon))
#     loi <- seq(1,lni2)
#     lai <- seq(1,lni2)
#     di <- seq(1,lni1)
#     ni <- seq(1,lni1)
#     loi <- c(loi[seq(ceiling(mulc/2),length(loi),mulc)],pi2)
#     lai <- c(lai[seq(ceiling(mulc/2), length(lai),mulc)],pi2)
#     di <- c(di[seq(ceiling(mulc/2), length(di),mulc)],pi1)
#     ni <- c(ni[seq(ceiling(mulc/2), length(ni),mulc)],pi1)
#     echam$lon <- echam$lon[loi]
#     echam$lat <- echam$lat[lai]
#     echam$data <- echam$data[di,,]
#     echam$ensmean <- echam$ensmean[di,]
#     echam$names <- echam$names[ni]
    if (vali) {
      loi <- seq(1:length(valiall$lon))
      lai <- seq(1:length(valiall$lat))
      di <- seq(1:dim(valiall$data)[1])
      ni <- seq(1:length(valiall$names))
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
      di <- di[seq(ceiling(mulc/2), length(di),mulc)]
      ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
      valiall$lon <- valiall$lon[loi]
      valiall$lat <- valiall$lat[lai]
      valiall$data <- valiall$data[di,,]
#      valiall$ensmean <- valiall$ensmean[di,]
      valiall$names <- valiall$names[ni]
#       lni1 <- (length(valiall$names)-length(echind$names))
#       lni2 <- (length(valiall$lon)-length(echind$names))# length no indices
#       pi1 <- seq(lni1+1,length(valiall$names)) # pos indices
#       pi2 <- seq(lni2+1,length(valiall$lon))
#       loi <- seq(1,lni2)
#       lai <- seq(1,lni2)
#       di <- seq(1,lni1)
#       ni <- seq(1,lni1)
#       loi <- c(loi[seq(ceiling(mulc/2),length(loi),mulc)],pi2)
#       lai <- c(lai[seq(ceiling(mulc/2), length(lai),mulc)],pi2)
#       di <- c(di[seq(ceiling(mulc/2), length(di),mulc)],pi1)
#       ni <- c(ni[seq(ceiling(mulc/2), length(ni),mulc)],pi1)
#       valiall$lon <- valiall$lon[loi]
#       valiall$lat <- valiall$lat[lai]
#       valiall$data <- valiall$data[di,,]
#       valiall$ensmean <- valiall$ensmean[di,]
#       valiall$names <- valiall$names[ni]
# #       lni1 <- (length(valiall$names)-length(echind$names))
# #       lni2 <- (length(valiall$lon)) #-length(echind$names))# length no indices
# #       pi1 <- seq(lni1+1,length(valiall$names)) # pos indices
# # #      pi2 <- seq(lni2,length(valiall$lon))
# #       loi <- seq(1,lni2)
# #       lai <- seq(1,lni2)
# #       di <- seq(1,lni1)
# #       ni <- seq(1,lni1)
# #       loi <- c(loi[seq(ceiling(mulc/2),length(loi),mulc)]) #,pi2)
# #       lai <- c(lai[seq(ceiling(mulc/2), length(lai),mulc)]) #,pi2)
# #       di <- c(di[seq(ceiling(mulc/2), length(di),mulc)],pi1)
# #       ni <- c(ni[seq(ceiling(mulc/2), length(ni),mulc)],pi1)
# # #       loi <- seq(1:length(valiall$lon))
# # #       lai <- seq(1:length(valiall$lat))
# # #       di <- seq(1:dim(valiall$data)[1])
# # #       ni <- seq(1:length(valiall$names))
# # #       loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
# # #       lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
# # #       di <- di[seq(ceiling(mulc/2), length(di),mulc)]
# # #       ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
# #       valiall$lon <- valiall$lon[loi]
# #       valiall$lat <- valiall$lat[lai]
# #       valiall$data <- valiall$data[di,]
# #       valiall$names <- valiall$names[ni]
    }
  }

## calc echam st. dev. for each grid point to scale docu data
## code below now because st. dev is better calculated based on 71 yrs of echam data than on 2 yrs
#   if (cyr == syr2) {
#     est <- array(echam$data,c(nrow(echam$data), 12, 
#                               ncol(echam$data)/12,30))
#     echam.sd.tmp <- array(0, dim(est)[c(1,2,4)])
#     echam.sd <- array(0, dim(est)[1:2])
#     for (i in 1:dim(est)[4]) {
#       for (j in 1:dim(est)[2]) {
#         echam.sd.tmp[,j,i] <- apply(est[,j,,i],1,sd)
#       }
#     }
#     echam.sd <- apply(echam.sd.tmp,1:2,mean)
#     rm(echam.sd.tmp)
#   }


  if (anomaly_assim){
#     nseas <- 12 # year with 12 months
#     asyr <- cyr-35
#     if (asyr < 1601) {asyr = 1601}
#     if (shortanom) { asyr = syr }
#     aeyr <- cyr+35
#     if (aeyr > 2005) {aeyr = 2005}
#     if (shortanom) { aeyr = eyr }
    if (load_71yr_anom) {
      yr1 <- cyr-1
      yr2 <- cyr
      yr3 <- yr1
      yr4 <- yr2
      if (cyr < 1637) {yr3 <- 1636} 
      if (cyr < 1637) {yr4 <- 1637} 
      if (cyr > 1970) {yr3 <- 1969} 
      if (cyr > 1970) {yr4 <- 1970} 
      if (every2grid) {
        load(paste0('../data/echam_anom/echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
        load(paste0('../data/echam_clim/echam_clim_',yr3,'-',yr4,'_2ndgrid.Rdata'))
      } else {
        load(paste0('../data/echam_anom/echam_anom_',yr1,'-',yr2,'.Rdata'))
        load(paste0('../data/echam_clim/echam_clim_',yr3,'-',yr4,'.Rdata'))
      }
      # just leave temp precip slp in state vector  
      if (tps_only) {
        tpspos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'), 
                    which(echam_anom$names=='slp'), which(echam_anom$names=='bias'))
        echam_anom$data <- echam_anom$data[tpspos,,]
        echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
        echam_anom$names <- echam_anom$names[tpspos]
        echam_clim$data <- echam_clim$data[tpspos,,]
        echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
        echam_clim$names <- echam_clim$names[tpspos]
      } else if (no_stream) {
        # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
        tpspos <- c(which(echam_anom$names!='stream'))
        echam_anom$data <- echam_anom$data[tpspos,,]
        echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
        echam_anom$names <- echam_anom$names[tpspos]
        echam_clim$data <- echam_clim$data[tpspos,,]
        echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
        echam_clim$names <- echam_clim$names[tpspos]
      }
      if (fasttest) {
        mulc <- 4 # choose every 4th grid box
        loi <- seq(1:length(echam_anom$lon))
        lai <- seq(1:length(echam_anom$lat))
        di <- seq(1:dim(echam_anom$data)[1])
        ni <- seq(1:length(echam_anom$names))
        loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
        lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
        di <- di[seq(ceiling(mulc/2), length(di),mulc)]
        ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
        echam_anom$lon <- echam_anom$lon[loi]
        echam_anom$lat <- echam_anom$lat[lai]
        echam_anom$data <- echam_anom$data[di,,]
        echam_anom$ensmean <- echam_anom$ensmean[di,]
        echam_anom$names <- echam_anom$names[ni]
        echam_clim$lon <- echam_clim$lon[loi]
        echam_clim$lat <- echam_clim$lat[lai]
        echam_clim$data <- echam_clim$data[di,,]
        echam_clim$ensmean <- echam_clim$ensmean[di,]
        echam_clim$names <- echam_clim$names[ni]
      }
      # error where echam_anom and echam_clim were generated: no unit correction happened, thus here
#      echam_anom$data[echam_anom$names=='temp2',,] <- echam_anom$data[echam_anom$names=='temp2',,] - 273.15
#      echam_anom$ensmean[echam_anom$names=='temp2',] <- echam_anom$ensmean[echam_anom$names=='temp2',] - 273.15
      echam_anom$data[echam_anom$names=='precip',,] <- echam_anom$data[echam_anom$names=='precip',,] * 3600 * 24 * 30
      echam_anom$ensmean[echam_anom$names=='precip',] <- echam_anom$ensmean[echam_anom$names=='precip',] * 3600 * 24 * 30
      echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,] / 100
      echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',] / 100
      echam.anom <- echam_anom
      echam_clim$data[echam_clim$names=='temp2',,] <- echam_clim$data[echam_clim$names=='temp2',,] - 273.15
      echam_clim$ensmean[echam_clim$names=='temp2',] <- echam_clim$ensmean[echam_clim$names=='temp2',] - 273.15
      echam_clim$data[echam_clim$names=='precip',,] <- echam_clim$data[echam_clim$names=='precip',,] * 3600 * 24 * 30
      echam_clim$ensmean[echam_clim$names=='precip',] <- echam_clim$ensmean[echam_clim$names=='precip',] * 3600 * 24 * 30
      echam_clim$data[echam_clim$names=='slp',,] <- echam_clim$data[echam_clim$names=='slp',,] / 100
      echam_clim$ensmean[echam_clim$names=='slp',] <- echam_clim$ensmean[echam_clim$names=='slp',] / 100
      echam.clim <- echam_clim
      rm(echam_anom,echam_clim)
#    } else {
    } else if (anom_reload) {
#        if (machine == "climcal3") {
      load(file=paste0("../data/anom/EnSRF_anom_",cyr,".Rdata")) 
#          if (cyr == syr2) load(file=paste0("../data/EnSRF_sd.Rdata"))
#        } else {  
#          load(file=paste0("../data/anom/EnSRF_anom_",cyr,".Rdata"))
#          if (cyr == syr2) load(file=paste0("../data/EnSRF_sd.Rdata"))
#          ffload(est,file=paste0("../data/EnSRF_est_",cyr,".Rdata"))
#        }
    } else {
    # calc running 71yr anomalies and climatology
    # 2. dimension contains e.g. 71 years of 12 monthly data 
      echam.backup <- echam
      if (machine == "climcal3") {
        etd <- array(0,dim=c(dim(echam$data)[1],((aeyr-asyr+1)*nseas),
                          dim(echam$data)[3]))
      } else {
        etd <- ff(0,dim=c(dim(echam$data)[1],((aeyr-asyr+1)*nseas),
                               dim(echam$data)[3]))
      }
      ete <- array(NA,dim=c(dim(echam$ensmean)[1],((aeyr-asyr+1)*nseas)))
      for (i in seq(0,(aeyr-asyr-1),by=1)) {
        yr1 <- asyr+i
        yr2 <- asyr+i+1
        if (every2grid) {
          load(paste0('../data/echam/echam_',yr1,'-',yr2,'_2ndgrid.Rdata'))
        } else {
          load(paste0('../data/echam/echam_',yr1,'-',yr2,'.Rdata'))
        }
        if (tps_only) {
          tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), 
                      which(echam$names=='slp'), which(echam$names=='bias'))
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
          echam$names <- echam$names[tpspos]
        } else if (no_stream) {
# ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
          tpspos <- c(which(echam$names!='stream'))
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
          echam$names <- echam$names[tpspos]
        } 
        if (fasttest) {
#             lni1 <- (length(echam$names)-length(echind$names))
#             lni2 <- (length(echam$lon)-length(echind$names))# length no indices
#             pi1 <- seq(lni1+1,length(echam$names)) # pos indices
#             pi2 <- seq(lni2+1,length(echam$lon))
#             loi <- seq(1,lni2)
#             lai <- seq(1,lni2)
#             di <- seq(1,lni1)
#             ni <- seq(1,lni1)
#             loi <- c(loi[seq(ceiling(mulc/2),length(loi),mulc)],pi2)
#             lai <- c(lai[seq(ceiling(mulc/2), length(lai),mulc)],pi2)
#             di <- c(di[seq(ceiling(mulc/2), length(di),mulc)],pi1)
#             ni <- c(ni[seq(ceiling(mulc/2), length(ni),mulc)],pi1)
#             echam$lon <- echam$lon[loi]
#             echam$lat <- echam$lat[lai]
#             echam$data <- echam$data[di,,]
#             echam$ensmean <- echam$ensmean[di,]
#             echam$names <- echam$names[ni]
          loi <- seq(1:length(echam$lon))
          lai <- seq(1:length(echam$lat))
          di <- seq(1:dim(echam$data)[1])
          ni <- seq(1:length(echam$names))
          loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
          lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
          di <- di[seq(ceiling(mulc/2), length(di),mulc)]
          ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
          echam$lon <- echam$lon[loi]
          echam$lat <- echam$lat[lai]
          echam$data <- echam$data[di,,]
          echam$ensmean <- echam$ensmean[di,]
          echam$names <- echam$names[ni]
        }
        if (i == (aeyr-asyr-1)) {
          etd[,((i*nseas+1):((i*nseas)+(2*nseas))),] <- echam$data[,1:(2*nseas),]
          ete[,((i*nseas+1):((i*nseas)+(2*nseas)))] <- echam$ensmean[,1:(2*nseas)]
        } else {
          etd[,((i*nseas+1):(i*nseas+nseas)),] <- echam$data[,1:nseas,]
          ete[,((i*nseas+1):(i*nseas+nseas))] <- echam$ensmean[,1:nseas]
        }
      } 
#     ti=which(floor(echam$time)>=(cyr-35) &
#                floor(echam$time)<=(cyr+35))
#     sts=ti[1]
#     ets=ti[length(ti)]
#     echam.tmp=echam
#     echam.tmp$data=echam$data[,sts:ets,]
#     echam.tmp$ensmean=echam$ensmean[,sts:ets]
#     echam.tmp$time=echam$time[sts:ets]
#    etd <- echam.tmp$data # 70 years of echam data
#    ete <- echam.tmp$ensmean # 70 years of echam ensmean
      nclim <- dim(etd)[2]/nseas # climatology is the 71 year mean of echam data
      nclime <- dim(ete)[2]/nseas # climatology is the 71 year mean of echam ensmean
    # calc climatology
      emn <- array(0, dim(etd[,1:nseas,]))
      emne <- array(0, dim(ete[,1:nseas]))
      for (i in 1:nclim){
        emn <- emn + etd[,(1:nseas + (i - 1)*nseas),]/nclim
        emne <- emne + ete[,(1:nseas + (i - 1)*nseas)]/nclime
      }
      eanom <- array(0, dim(etd))
      eanome <- array(0, dim(ete))
      for (i in 1:nclim){
        eanom[,1:nseas + (i - 1)*nseas,] <- etd[,1:nseas + (i - 1)*nseas,] - emn
        eanome[,1:nseas + (i - 1)*nseas] <- ete[,1:nseas + (i - 1)*nseas] - emne
#    eanom <- etd[,1:nseas + (ceiling(nclim/2) - 1)*nseas,] - emn
#    eanome <- ete[,1:nseas + (ceiling(nclime/2) - 1)*nseas,] - emne
       }
       echam.anom <- echam
       echam.anom$data <- eanom
       echam.anom$ensmean <- eanome
       echam.anom$time <- seq(asyr,aeyr+11/12,by=1/12)
    #    echam.anom$data <-    array(echam$data - as.vector(apply(array(echam.tmp$data, 
    #                            c(nrow(echam.tmp$data), 12, ncol(echam.tmp$data)/12,30)), 1:2,
    #                            mean)), c(nrow(echam.tmp$data), 12, ncol(echam.tmp$data)/12, 30))
    #    echam.anom$data <-    (array(echam.anom$data,c(dim(echam.anom$data)[1],                              
    #                              dim(echam.anom$data)[2]*dim(echam.anom$data)[3],
    #                              dim(echam.anom$data)[4])))
    #     echam.anom$ensmean <- array(echam$ensmean - as.vector(apply(array(echam.tmp$ensmean,
    #                             c(nrow(echam.tmp$ensmean), 12, ncol(echam.tmp$ensmean)/12)), 1:2, 
    #                             mean)), c(nrow(echam.tmp$ensmean), 12, ncol(echam.tmp$ensmean)/12))
    #     echam.anom$ensmean <- (array(echam.anom$ensmean,c(dim(echam.anom$ensmean)[1],
    #                             dim(echam.anom$ensmean)[2]*dim(echam.anom$ensmean)[3])))  
      echam.clim <- echam
      echam.clim$data <- emn
      echam.clim$ensmean <- emne
      echam.clim$time <- echam.anom$time
    #    echam.clim$data <-    apply(array(echam.tmp$data, c(nrow(echam.tmp$data), 12, 
    #                            ncol(echam.tmp$data)/12,30)), c(1,2,4), mean)
    #     echam.clim$ensmean <- apply(array(echam.tmp$ensmean, c(nrow(echam.tmp$ensmean),
    #                             12, ncol(echam.tmp$ensmean)/12)), 1:2, mean) 
#        if (machine == "climcal3") {
#  if (cyr == syr2) {
      echam <- echam.backup
      rm(echam.backup) # ;rm(echam.tmp);
      rm(etd);rm(emn);rm(emne);rm(ete);rm(eanom);rm(eanome) # ;rm(est)
      if(anom_save) {
       save(echam.anom,echam.clim,file=paste0("../data/anom/EnSRF_anom_",cyr,".Rdata"))
      }
    }
  }
#    save(echam.sd,file=paste0("../data/EnSRF_sd.Rdata"))
#        } else {
#          save(echam.anom,echam.clim,est,file=paste0("../data/EnSRF_anom_",cyr,".Rdata"))
#  }  
#} else {
#  save(echam.anom,echam.clim,echam.sd,file=paste0("../data/EnSRF_anom_",cyr,".Rdata"))
#  if (cyr == syr2) {
#    save(echam.sd,file=paste0("../data/EnSRF_sd.Rdata"))
#    ffsave(est,file=paste0("../data/EnSRF_est_",cyr,".Rdata"))
#  }

#        est <- array(etd,c(nrow(etd), 12, ncol(etd)/12,30))

  print('calc time for echam anomalies')
  print(proc.time() - ptm1)

# calc echam st. dev. for each grid point and month over ens memb. to scale docu data
#        if (cyr == syr2) {
  echam.sd <- apply(echam$data[,13:24,],1:2,sd)
# 
#   echam.sd.tmp <- array(0, dim(est)[c(1,2,4)])
#   echam.sd <- array(0, dim(est)[1:2])
#   for (i in 1:dim(est)[4]) {
#     for (j in 1:dim(est)[2]) {
#               echam.sd.tmp[,j,i] <- apply(est[,j,,i],1,sd)
#             }
#           }
#           echam.sd <- apply(echam.sd.tmp,1:2,mean)
#           rm(echam.sd.tmp)
#         } else {
#           echam.sd <- array(NA,dim=c(nrow(echam$data),12))
#         }
  print('calc time for standard deviations')
  print(proc.time() - ptm1)


#########################################################################################
# screen proxy/instr. assimilation data 
#########################################################################################
  if ((instrumental) & (check_assimdata)) {
    #    screendata(inst_t,echam)
    if (ghcn_prec) {
      varlist <- c("inst_t","inst_slp","ghcn","ghcn_precip")
    } else {
      varlist <- c("inst_t","inst_slp","ghcn")
    }
    for (varname in varlist) {
      var <- get(varname)
      gpos <- getgridboxnum(var,echam)
      if (!is.na(gpos)) {
        for (i in 1:length(gpos)) {
          m <- gpos[i]
          d <- compute_dist(var$lon[i],var$lat[i],echam$lon[m],echam$lat[m])
          if ((!is.na(d)) & (d > 600)) {
            m=NA
            print(paste('inst data', varname, i, '>600km from echam grid box; set to NA'))
            if (varname=="inst_t") {inst_t$data[,i] <- NA}
            if (varname=="inst_slp") {inst_slp$data[,i] <- NA}
            if (varname=="ghcn") {ghcn$data[,i] <- NA}
            if (ghcn_prec){
              if (varname=="ghcn_precip") {ghcn_precip$data[,i] <- NA}
            }
          }
          if (!is.na(m)) {
            tiv=which(floor(var$time)==cyr)
            stsv=tiv[1]
            vtmp <- array(var$data,c(12, nrow(var$data)/12, dim(var$data)[2]))
            for (j in 1:12) {
              # if bias corrected proxy/inst is outside echam ens range +- 3SD, 
              # data point will not be assimilated at this time step
#               biasm <- mean(est[m,j,,]) - mean(vtmp[j,,i])
              biasm <- echam.clim$ensmean[m,j] - mean(vtmp[j,,i])            
#              print(c(varname,i,j,tiv,biasm))
              if (!is.na(biasm)) {
                if (((vtmp[j,((stsv-1)/12+1),i]+ biasm) < echam$ensmean[m,(j+12)]-5*echam.sd[m,j]) | 
                   ((vtmp[j,((stsv-1)/12+1),i]+ biasm) > echam$ensmean[m,(j+12)]+5*echam.sd[m,j])) {
                  print(paste('inst data', varname, i, j, 'out of range'))
#                  stop()
                  if (cyr == syr2) {
                    write(paste('inst data', varname, i, j, 'out of range'),
                            file='screening_instr.log',append=F)
                  
                  } else {
                    write(paste('inst data', varname, i, j, 'out of range'),
                            file='screening_instr.log',append=T)
                  }
                  write(paste('data lon/lat', c(var$lon[i],var$lat[i])),
                        file='screening_instr.log',append=T)
                  write(paste('echam lon/lat', c(echam$lon[m],echam$lat[m])),
                        file='screening_instr.log',append=T)
                  write(paste('bias corr. data', vtmp[j,((stsv-1)/12+1),i]+ biasm),
                        file='screening_instr.log',append=T)
                  write(paste('echam mean', echam$ensmean[m,(j+12)]),
                        file='screening_instr.log',append=T)
                  write(paste('echam sd', echam.sd[m,j]),
                        file='screening_instr.log',append=T)
                  if (varname=="inst_t") {inst_t$data[,i] <- NA}
                  if (varname=="inst_slp") {inst_slp$data[,i] <- NA}
                  if (varname=="ghcn") {ghcn$data[,i] <- NA}
                  if (ghcn_prec){
                    if (varname=="ghcn_precip") {ghcn_precip$data[,i] <- NA}
                  }
                }
              }  
            }    
          }
        }
      } 
    }
  }

  if (calc_decorr_dist) {
    d <- compute_dist_2d(echam$lon,echam$lat,echam$lon,echam$lat)
    for (i in unique(echam$names)) {
#      tmp <- echam.anom$data[echam$names==i,,]
#      corens <- cor(t(tmp[,,1]))
      tmp <- echam.anom$ensmean[echam$names==i,]
      corens <- cor(t(tmp[,]))
      
#      pdf(paste0('../figures/decorr_',i,'.pdf'), width=9, height=4.5, paper='special')
      png(paste0('../figures/decorr_',i,'.png'), width = 1024, height = 768)
#      plot(as.vector(d[1:100000]),as.vector(corens[1:100000]),col='#4c8bff10', xlim=c(0,5000),ylim=c(0,1))
        plot(as.vector(d),as.vector(corens),col='#4c8bff01',
             xlim=c(0,5000),ylim=c(0,1))
        lines(exp(-1/2 * (seq(1:5000)/get(paste0('l_dist_',i)))**2),col='red')
      dev.off()
    }
  }


  # no check for docu data as already good quality
  if ((real_proxies) & (check_assimdata)) {
  # correlation screening already where multiple regression coefficients are calculated
  # thus, screen if value at current time step is more than 5 std. dev. from mean
  # in this case treated as outlier and set to NA
    for (i in 1:length(realprox$lon)) {
      tiv=which(floor(realprox$time)==cyr)
      rpmean <- mean(realprox$data[,i],na.rm=T)
      rpsd   <- sd(realprox$data[,i],na.rm=T)
      if ((!is.na(rpmean)) & (!is.na(rpsd))) {
        if ((!is.na(realprox$data[tiv,i])) & ((realprox$data[tiv,i] < rpmean-5*rpsd) 
                                           | (realprox$data[tiv,i] > rpmean+5*rpsd))) {
          realprox$data[,i] <- NA
          print(paste('proxy data', i, 'out of range'))
          if (cyr == syr2) {
            write(paste('proxy data', i, 'out of range'),file='screening_proxies.log',append=F)
          } else {
            write(paste('proxy data', i, 'out of range'),file='screening_proxies.log',append=T)
          }
        }  
      }
    }
  }

  print('calc time for screening proxies')
  print(proc.time() - ptm1)

  if (ind_ECHAM) { 
  # add indices to echam$data 
    if (anomaly_assim) {
#      if (load_71yr_anom) { 
      ti=c(which(floor(echind$time)==cyr-1),which(floor(echind$time)==cyr))
      sts=ti[1]
      ets=ti[length(ti)]
      echind$data=echind$data[,sts:ets,]
      echind$ensmean=echind$ensmean[,sts:ets]
      echind$time=echind$time[sts:ets]
 #     }
      echam.anom$data <- abind(echam.anom$data,echind$data,along=1)
      echam.anom$ensmean <- rbind(echam.anom$ensmean,apply(echind$data,1:2,mean))
      echam.anom$names <- c(echam.anom$names,echind$names)
      echam.anom$lon <- c(echam.anom$lon,rep(NA,length(echind$names)))
      echam.anom$lat <- c(echam.anom$lat,rep(NA,length(echind$names)))
      if (load_71yr_anom) { 
        echam.clim$data <- abind(echam.clim$data,array(NA,dim=c(6,24,30)),along=1)
        echam.clim$ensmean <- rbind(echam.clim$ensmean,array(NA,dim=c(6,24)))
      } else {
        echam.clim$data <- abind(echam.clim$data,array(NA,dim=c(6,12,30)),along=1)
        echam.clim$ensmean <- rbind(echam.clim$ensmean,array(NA,dim=c(6,12)))
      }  
      echam.clim$names <- c(echam.clim$names,echind$names)
      echam.clim$lon <- c(echam.clim$lon,rep(NA,length(echind$names)))
      echam.clim$lat <- c(echam.clim$lat,rep(NA,length(echind$names))) 
    }
#     ti=which(floor(echind$time)>=cyr-1 & floor(echind$time)<=cyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     echind$data <- echind$data[,sts:ets,]
#     echind$ensmean <- echind$ensmean[,sts:ets]
#     echind$time=echind$time[sts:ets]
    echam$data <- abind(echam$data,echind$data,along=1)
    echam$ensmean <- rbind(echam$ensmean,apply(echind$data,1:2,mean))
    echam$names <- c(echam$names,echind$names)
    echam$lon <- c(echam$lon,rep(NA,length(echind$names)))
    echam$lat <- c(echam$lat,rep(NA,length(echind$names)))
  }
  if (ind_recon) {
    valiall$data <- rbind(valiall$data[1:dim(valiall$data)[1],,1],indall) #,valiall$data[2083:2087,,])
    valiall$names <- c(valiall$names,rownames(indall)) 
    valiall$lon <- c(valiall$lon,rep(NA,nrow(indall))) 
    valiall$lat <- c(valiall$lat,rep(NA,nrow(indall))) 
  }

# just leave data for one year (max 12 months) and correct time resolution in memory
  if (sixmonstatevector) {
    # change array to have 6 months in state vector for winter and summer
    # first winter starts in oct of syr
    # 6 mon stat vectors for oct-mar and apr and sep
    if (!anomaly_assim) {
      tmp1 <- array(echam$data,c(dim(echam$data)[1]*dim(echam$data)[2],dim(echam$data)[3]))
      tmp2 <- tmp1[((9*dim(echam$data)[1]+1):(dim(tmp1)[1]-(3*dim(echam$data)[1]))),] # cut oct syr to sep eyr
      echam$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(echam$data)[2]/12)-1)*2),
                                 (((dim(echam$data)[2]/12)-1)*2),dim(tmp2)[2])) 
      # reconvert to 2 seasons per year
      tmp11 <- array(echam$ensmean,c(dim(echam$ensmean)[1]*dim(echam$ensmean)[2]))
      tmp12 <- tmp11[((9*dim(echam$ensmean)[1]+1):(dim(tmp11)[1]-(3*dim(echam$ensmean)[1])))] # cut oct syr to sep eyr
      echam$ensmean <- array(tmp12,c(dim(tmp12)[1]/(((dim(echam$ensmean)[2]/12)-1)*2),
                          (((dim(echam$ensmean)[2]/12)-1)*2))) # reconvert to 2 seasons per year
#      echam$time <- seq(syr+1,eyr+0.5,0.5) 
      echam$time <- seq(cyr,cyr+1.5,0.5) 
#      echam$names <- rep(echam$names,6)
      rm(tmp1);rm(tmp2);rm(tmp11);rm(tmp12)
    }
    if (anomaly_assim) {
      tmp41 <- array(echam.anom$data,c(dim(echam.anom$data)[1]*dim(echam.anom$data)[2],
                                       dim(echam.anom$data)[3]))
      tmp42 <- tmp41[((9*dim(echam.anom$data)[1]+1):(dim(tmp41)[1]-(3*dim(echam.anom$data)[1]))),]
      echam.anom$data <- array(tmp42,c(dim(tmp42)[1]/(((dim(echam.anom$data)[2]/12)-1)*2),
                           (((dim(echam.anom$data)[2]/12)-1)*2),dim(tmp42)[2])) 
      tmp51 <- array(echam.anom$ensmean,c(dim(echam.anom$ensmean)[1]*dim(echam.anom$ensmean)[2]))
      tmp52 <- tmp51[((9*dim(echam.anom$ensmean)[1]+1):(dim(tmp51)[1]-
                 (3*dim(echam.anom$ensmean)[1])))] 
      echam.anom$ensmean <- array(tmp52,c(dim(tmp52)[1]/(((dim(echam.anom$ensmean)[2]/12)-1)*2),
                              (((dim(echam.anom$ensmean)[2]/12)-1)*2))) 
      if (load_71yr_anom) { 
        echam.anom$time <- c(cyr,cyr+0.5) 
      } else {  
        echam.anom$time <- seq(asyr+1,aeyr+0.5,0.5) 
      }
#      echam.anom$names <- rep(echam.anom$names,6)
      if (load_71yr_anom) { 
        tmp61 <- array(echam.clim$data,c(dim(echam.clim$data)[1]*dim(echam.clim$data)[2],
                                         dim(echam.clim$data)[3]))
        tmp62 <- tmp61[((9*dim(echam.clim$data)[1]+1):(dim(tmp61)[1]-(3*dim(echam.clim$data)[1]))),]
        echam.clim$data <- array(tmp62,c(dim(tmp62)[1]/(((dim(echam.clim$data)[2]/12)-1)*2),
                                         (((dim(echam.clim$data)[2]/12)-1)*2),dim(tmp62)[2])) 
        tmp71 <- array(echam.clim$ensmean,c(dim(echam.clim$ensmean)[1]*dim(echam.clim$ensmean)[2]))
        tmp72 <- tmp71[((9*dim(echam.clim$ensmean)[1]+1):(dim(tmp71)[1]-
                                                            (3*dim(echam.clim$ensmean)[1])))] 
        echam.clim$ensmean <- array(tmp72,c(dim(tmp72)[1]/(((dim(echam.clim$ensmean)[2]/12)-1)*2),
                                            (((dim(echam.clim$ensmean)[2]/12)-1)*2))) 
        echam.clim$time <- c(cyr,cyr+0.5) 
      } else { 
        tmp61 <- array(echam.clim$data,c(dim(echam.clim$data)[1]*dim(echam.clim$data)[2],
                                       dim(echam.clim$data)[3]))
      # Special because climatology is only 1 year lang -> reorder months
        tmp62 <- tmp61[((9*dim(echam.clim$data)[1]+1):(dim(tmp61)[1])),] # cut oct syr to sep eyr
        tmp63 <- tmp61[(1:(9*dim(echam.clim$data)[1])),]
        tmp64 <- rbind(tmp62,tmp63)
        echam.clim$data <- array(tmp64,c(dim(tmp64)[1]/(((dim(echam.clim$data)[2]/12))*2),
                                       (((dim(echam.clim$data)[2]/12))*2),dim(tmp64)[2])) 
        tmp71 <- array(echam.clim$ensmean,c(dim(echam.clim$ensmean)[1]*dim(echam.clim$ensmean)[2]))
        tmp72 <- tmp71[((9*dim(echam.clim$ensmean)[1]+1):(dim(tmp71)[1]))] # cut oct syr to sep eyr
        tmp73 <- tmp71[(1:(9*dim(echam.clim$ensmean)[1]))]
        tmp74 <- c(tmp72,tmp73)
        echam.clim$ensmean <- array(tmp74,c(length(tmp74)/(((dim(echam.clim$ensmean)[2]/12))*2),
                                          (((dim(echam.clim$ensmean)[2]/12))*2)))  
        echam.clim$time <- seq(syr+1,eyr+0.5,0.5) 
      }
      echam.clim$names <- rep(echam.clim$names,6)
      rm(tmp41);rm(tmp42);rm(tmp51);rm(tmp52);rm(tmp61);rm(tmp62);rm(tmp71);rm(tmp72)
      if (!load_71yr_anom) {rm(tmp63);rm(tmp64);rm(tmp73);rm(tmp74)}
    }
    if (!recon_vali & vali) {
      tmp21 <- array(valiall$data,c(dim(valiall$data)[1]*dim(valiall$data)[2]))
      tmp22 <- tmp21[((9*dim(valiall$data)[1]+1):(dim(tmp21)[1]-(3*dim(valiall$data)[1])))] 
      valiall$data <- array(tmp22,c(dim(tmp22)[1]/(((dim(valiall$data)[2]/12)-1)*2),
                        (((dim(valiall$data)[2]/12)-1)*2))) # reconvert to 2 seasons per year
      valiall$time <- seq(syr+1,eyr+0.5,0.5) 
      valiall$names <- rep(valiall$names,6) 
      rm(tmp21);rm(tmp22)
    } else if (vali) {
      pos <- sort(c(agrep('.042',as.character(valiall$time)), agrep('.542',as.character(valiall$time))))
      valiall$time <- round(valiall$time[pos],digits=1)
      valiall$data <- valiall$data[,pos]      
    }
  }
  
# cut one year time slice  
#  echam.allts=echam
  if (anomaly_assim){
    if (load_71yr_anom) {
      echam=echam.anom
      rm(echam.anom)
    } else {  
      ti=which(floor(echam.anom$time)==cyr)
      sts=ti[1]
      ets=ti[length(ti)]
      echam$data=echam.anom$data[,sts:ets,]
      echam$ensmean=echam.anom$ensmean[,sts:ets]
      echam$time=echam.anom$time[sts:ets]
      rm(echam.anom)
    }
  } else {  
    ti=which(floor(echam$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    echam$data=echam$data[,sts:ets,]
    echam$ensmean=echam$ensmean[,sts:ets]
    echam$time=echam$time[sts:ets]    
  }   
  
  if (vali) {
    valiall.allts=valiall
    if (cru_vali) {
      if (cyr>1900 && cyr<2005){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }
    }
    if (ncep_vali) {
      if (cyr>1947 && cyr<2010){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }   
    }
    if (recon_vali) {
      if (cyr>1750 && cyr<1901){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }  
    }
  }

  if (instrumental) {
    if ((ghcn_temp) & (dim(ghcn$data)[2]>0)) {
      if (anomaly_assim){
        ti=which(floor(ghcn$time)>=(cyr-35) &
                   floor(ghcn$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn.tmp=ghcn
        ghcn.tmp$data=t(ghcn$data[sts:ets,])
        ghcn.tmp$time=ghcn$time[sts:ets]
        ghcn.anom <- ghcn
        ghcn.anom$data <- (t(ghcn$data) - matrix(rep(apply(array(ghcn.tmp$data,
           c(nrow(ghcn.tmp$data), 12, ncol(ghcn.tmp$data)/12)), 1:2, mean,na.rm=T),
           (length(ghcn.tmp$time)/12)), nrow=nrow(t(ghcn$data))))
        ghcn <- ghcn.anom
        ghcn$data <- t(ghcn.anom$data)
      }
    }
    if (ghcn_prec) {
     if (dim(ghcn_precip$data)[2]>0) {  
      if (anomaly_assim){
        ti=which(floor(ghcn_precip$time)>=(cyr-35) &
                   floor(ghcn_precip$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn_precip.tmp=ghcn_precip
        ghcn_precip.tmp$data=t(ghcn_precip$data[sts:ets,])
        ghcn_precip.tmp$time=ghcn_precip$time[sts:ets]
        ghcn_precip.anom <- ghcn_precip
        ghcn_precip.anom$data <- (t(ghcn_precip$data) - matrix(rep(apply(array(ghcn_precip.tmp$data,
          c(nrow(ghcn_precip.tmp$data), 12, ncol(ghcn_precip.tmp$data)/12)), 1:2, mean,na.rm=T),
          (length(ghcn_precip.tmp$time)/12)), nrow=nrow(t(ghcn_precip$data)))) 
        ghcn_precip <- ghcn_precip.anom
        ghcn_precip$data <- t(ghcn_precip.anom$data)
      }
     }  
    }
    if ((yuri_slp) & (dim(inst_slp$data)[2]>0)) {  
      if (anomaly_assim){
        ti=which(floor(inst_slp$time)>=(cyr-35) &
                   floor(inst_slp$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_slp.tmp=inst_slp
        inst_slp.tmp$data=t(inst_slp$data[sts:ets,])
        inst_slp.tmp$time=inst_slp$time[sts:ets]
        inst_slp.anom <- inst_slp
        inst_slp.anom$data <- (t(inst_slp$data) - matrix(rep(apply(array(inst_slp.tmp$data,
          c(nrow(inst_slp.tmp$data), 12, ncol(inst_slp.tmp$data)/12)), 1:2, mean,na.rm=T),
          (length(inst_slp.tmp$time)/12)), nrow=nrow(t(inst_slp$data)))) 
        inst_slp <- inst_slp.anom
        inst_slp$data <- t(inst_slp.anom$data)
      }
    }
    if ((yuri_temp) & (dim(inst_t$data)[2]>0)) {
      if (anomaly_assim){
        ti=which(floor(inst_t$time)>=(cyr-35) &
                   floor(inst_t$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_t.tmp=inst_t
        inst_t.tmp$data=t(inst_t$data[sts:ets,])
        inst_t.tmp$time=inst_t$time[sts:ets]
        inst_t.anom <- inst_t
        inst_t.anom$data <- (t(inst_t$data) - matrix(rep(apply(array(inst_t.tmp$data,
          c(nrow(inst_t.tmp$data), 12, ncol(inst_t.tmp$data)/12)), 1:2, mean,na.rm=T),
          (length(inst_t.tmp$time)/12)), nrow=nrow(t(inst_t$data))))
        inst_t <- inst_t.anom
        inst_t$data <- t(inst_t.anom$data)
      }
    }
    if (ghcn_prec) {  
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0) & 
          (dim(ghcn_precip$data)[2]==0)) {instrumental=F}
    } else {
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0)) {
        instrumental=F}
    }  
  }


  if (docu) {
    if (anomaly_assim){
      ti=which(floor(doc_t_mon$time)>=(cyr-35) &
                 floor(doc_t_mon$time)<=(cyr+35))
      sts=ti[1]
      ets=ti[length(ti)]
      doc_t.tmp=doc_t_mon
      doc_t.tmp$data=t(doc_t_mon$data[sts:ets,])
      doc_t.tmp$time=doc_t_mon$time[sts:ets]
      doc_t_mon.anom <- doc_t_mon
      doc_t_mon.anom$data <- (t(doc_t_mon$data) - matrix(rep(apply(array(doc_t.tmp$data,
        c(nrow(doc_t.tmp$data), 12, ncol(doc_t.tmp$data)/12)), 1:2, mean,na.rm=T),
        (length(doc_t.tmp$time)/12)), nrow=nrow(t(doc_t_mon$data))))
      doc_t_mon <- doc_t_mon.anom
      doc_t_mon$data <- t(doc_t_mon.anom$data)
      doc_t_mon$season <- seq(1,12)
      if (scaleprox){
        # calc echam variability at data location
        echamatprox.arr <- array(NA,c(length(doc_t_mon$lon), length(doc_t_mon$time)))
        for(i in 1:(length(doc_t_mon$lon))){
#          print(i)
          plon <- doc_t_mon$lon[i]
          plat <- doc_t_mon$lat[i]
          clon <- echam$lon[!is.na(echam$lon)]
          clat <- echam$lat[!is.na(echam$lat)]
          k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
          l=which(abs(clat-plat)==min(abs(clat-plat)))
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat)==min(abs(clat-plat)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat)==min(abs(clat-plat)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          m=k[which(match(k,l)>0)]
#          print(m)
          if (length(m) > 0) {
            if (doc_t_mon$names[i]=="precip") { 
              m <- m+length(echam$lon) 
            } else if (doc_t_mon$names[i]=="slp") { 
              m <- m+(2*length(echam$lon)) 
            }
            scalefac <- 1/echam.sd[m,]
          } else {
            scalefac <- NA
          }
        }
        doc_t_mon$data <- t(scale(t(doc_t_mon$data),center=F,scale=rep(scalefac,
                             (dim(doc_t_mon$data)[1]/12))))
      }
      
#       ti=which(floor(doc_t_seas$time)>=(cyr-35) &
#                  floor(doc_t_seas$time)<=(cyr+35))
#       sts=ti[1]
#       ets=ti[length(ti)]
#       doc_t.tmp=doc_t_seas
#       doc_t.tmp$data=t(doc_t_seas$data[sts:ets,])
#       doc_t.tmp$time=doc_t_seas$time[sts:ets]
#       doc_t_seas.anom <- doc_t_seas
#       doc_t_seas.anom$data <- (t(doc_t_seas$data) - matrix(rep(apply(array(doc_t.tmp$data,
#         c(nrow(doc_t.tmp$data), 4, ncol(doc_t.tmp$data)/4)), 1:2, mean,na.rm=T),
#         (length(doc_t.tmp$time)/4)), nrow=nrow(t(doc_t_seas$data))))
#       doc_t_seas <- doc_t_seas.anom
#       doc_t_seas$data <- t(doc_t_seas.anom$data)
#       doc_t_seas$season <- c('win','spr','sum','aut')
#       if (scaleprox){
#         # calc echam variability at data location
#         echamatprox.arr <- array(NA,c(length(doc_t_seas$lon), length(doc_t_seas$time)))
#         for(i in 1:(length(doc_t_seas$lon))){
#           plon <- doc_t_seas$lon[i]
#           plat <- doc_t_seas$lat[i]
#           clon <- echam$lon[!is.na(echam$lon)]
#           clat <- echam$lat[!is.na(echam$lat)]
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           m=k[which(match(k,l)>0)]
#           if (length(m) > 0) {
#             if (doc_t_seas$names[i]=="precip") { 
#               m <- m+length(echam$lon) 
#             } else if (doc_t_seas$names[i]=="slp") { 
#               m <- m+(2*length(echam$lon)) 
#             }
#             scalefac <- 1/c(mean(echam.sd[m,c(12,1,2)]),mean(echam.sd[m,c(3,4,5)]),
#                             mean(echam.sd[m,c(6,7,8)]),mean(echam.sd[m,c(9,10,11)]))
#           } else {
#             scalefac <- NA
#           }
#         }
#         doc_t_seas$data <- t(scale(t(doc_t_seas$data),center=F,scale=rep(scalefac,
#                                  (dim(doc_t_seas$data)[1]/4))))
#       }
# 
# # just scale if one value per year!      
#       ti=which(floor(doc_t_JFMA$time)>=(cyr-35) &
#                  floor(doc_t_JFMA$time)<=(cyr+35))
#       sts=ti[1]
#       ets=ti[length(ti)]
#       doc_t.mean=mean(doc_t_JFMA$data[sts:ets])
#       doc_t_JFMA$data <- doc_t_JFMA$data - doc_t.mean
#       doc_t_JFMA$season <- 'JFMA'
#       if (scaleprox){
#         # calc echam variability at data location
#         echamatprox.arr <- array(NA,c(length(doc_t_JFMA$lon), length(doc_t_JFMA$time)))
#         for(i in 1:(length(doc_t_JFMA$lon))){
#           plon <- doc_t_JFMA$lon[i]
#           plat <- doc_t_JFMA$lat[i]
#           clon <- echam$lon[!is.na(echam$lon)]
#           clat <- echam$lat[!is.na(echam$lat)]
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           m=k[which(match(k,l)>0)]
#           if (length(m) > 0) {
#             if (doc_t_JFMA$names[i]=="precip") { 
#               m <- m+length(echam$lon) 
#             } else if (doc_t_JFMA$names[i]=="slp") { 
#               m <- m+(2*length(echam$lon)) 
#             }
#             scalefac <- 1/mean(echam.sd[m,1:4])
#           } else {
#             scalefac <- NA
#           }
#         }
#         if (!is.na(scalefac)) {
#           doc_t_JFMA$data <- as.vector(scale(doc_t_JFMA$data,center=F,scale=scalefac))
#         }
#       }
#     } 
# 
#     ti=which(floor(doc_t_AMJJA$time)>=(cyr-35) &
#            floor(doc_t_AMJJA$time)<=(cyr+35))
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t.mean=mean(doc_t_AMJJA$data[sts:ets])
#     doc_t_AMJJA$data <- doc_t_AMJJA$data - doc_t.mean
#     doc_t_AMJJA$season <- 'AMJJA'
#     if (scaleprox){
#       # calc echam variability at data location
#       echamatprox.arr <- array(NA,c(length(doc_t_AMJJA$lon), length(doc_t_AMJJA$time)))
#       for(i in 1:(length(doc_t_AMJJA$lon))){
#         plon <- doc_t_AMJJA$lon[i]
#         plat <- doc_t_AMJJA$lat[i]
#         clon <- echam$lon[!is.na(echam$lon)]
#         clat <- echam$lat[!is.na(echam$lat)]
#         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#         l=which(abs(clat-plat)==min(abs(clat-plat)))
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         m=k[which(match(k,l)>0)]
#         if (length(m) > 0) {
#           if (doc_t_AMJJA$names[i]=="precip") { 
#             m <- m+length(echam$lon) 
#           } else if (doc_t_AMJJA$names[i]=="slp") { 
#             m <- m+(2*length(echam$lon)) 
#           }
#           scalefac <- 1/mean(echam.sd[m,4:8])
#         } else {
#           scalefac <- NA
#         }
#       }
#       if (!is.na(scalefac)) {
#         doc_t_AMJJA$data <- as.vector(scale(doc_t_AMJJA$data,center=F,scale=scalefac))
#       }
    }
  }





  if (real_proxies) {
    # no scaling because regresion takes care of it
    #    realprox$data <- scale(realprox$data,center=T,scale=T)
    realprox.allts <- realprox
    tmp1=t(realprox$data)
    tmp2=array(NA,c(dim(tmp1)[1],2,dim(tmp1)[2]))
    tmp2[,2,]=tmp1
    realprox.allts$data=array(tmp2,c(dim(tmp2)[1],dim(tmp2)[2]*dim(tmp2)[3]))
    realprox.allts$data=realprox.allts$data[,3:dim(realprox.allts$data)[2]]
    realprox.allts$time=seq(syr+1,eyr+0.5,0.5) 
    ti=which(floor(realprox.allts$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]      
    realprox$data=realprox.allts$data[,sts:ets]
    realprox$time=realprox.allts$time[sts:ets]
    realprox$names=rep("prox",dim(realprox.allts$data)[1])
    if (anomaly_assim){
      ti=which(floor(realprox.allts$time)>=(cyr-35) &
                 floor(realprox.allts$time)<=(cyr+35))
      sts=ti[1]
      ets=ti[length(ti)]
      realprox.tmp=realprox.allts
      realprox.tmp$data=realprox.allts$data[,sts:ets]
      realprox.tmp$time=realprox.allts$time[sts:ets]
      realprox.anom <- realprox
      realprox.anom$data <- (realprox$data - (apply(array(realprox.tmp$data,
        c(nrow(realprox.tmp$data), 2, ncol(realprox.tmp$data)/2)), 1:2, mean,na.rm=T)))
      realprox <- realprox.anom
    }
#     # write_coor done in EnSRF_prepplot script
#     if (write_coor) {
#       # proxy data coordinates at each time step
#       lonlat_year=list()
#       for (t in 1:length(realprox$time)) {
#         tf=as.logical(!is.na(realprox$data[,t]))
#         outstr <- cbind(realprox$lon[tf],realprox$lat[tf],realprox$names[tf])
#         outprox <- outstr[outstr[,3]=="prox",]
#         outprox2 <- unique(t(apply(outprox[,1:2],1,paste0)))
#         outstat <- outstr[outstr[,3]!="prox",]
#         outstat2 <- unique(t(apply(outstat[,1:2],1,paste0)))
#         write.table(outstat2, file=paste0('../data/coor/stat_coor_',realprox$time[t],'.csv'),
#                     row.names=F,col.names=F)
#         write.table(outprox2, file=paste0('../data/coor/prox_coor_',realprox$time[t],'.csv'),
#                     row.names=F,col.names=F)
#       }
#     }
  }  


#########################################################################################
# prepare proxy/instr. assimilation data
#########################################################################################
  
# real trw proxy multiple regression approach
  if ((real_proxies) & (!instrumental)) {  # & (!docu)) {
    if (reduced_proxies) {
      every <- 12 
      redpos <-seq(1,length(realprox$lon),every)
      realprox$lon <- realprox$lon[redpos]
      realprox$lat <- realprox$lat[redpos]
      realprox$data <- realprox$data[redpos,]
      realprox$names <- realprox$names[redpos]
      realprox.allts$lon <- realprox.allts$lon[redpos]
      realprox.allts$lat <- realprox.allts$lat[redpos]
      realprox.allts$data <- realprox.allts$data[redpos,]
      realprox.allts$names <- realprox.allts$names[redpos]
    }
#     proxies <- list(lon=realprox$lon, lat=realprox$lat, data=realprox$data, 
#                names=realprox$names, time=realprox$time)
#     proxies.allts <- list(lon=realprox.allts$lon, lat=realprox.allts$lat, 
#                      data=realprox.allts$data, names=realprox.allts$names, 
#                      time=realprox.allts$time)
    proxies<-list(data=realprox$data, lon=realprox$lon, 
                    lat=realprox$lat, names=realprox$names, 
                    height=realprox$elevation, time=realprox$time,
                    mr=realprox$mr, var_residu=realprox$var_residu,
                    numavg=rep(1,length(realprox$lon)))
    proxies.allts<-list(data=realprox.allts$data, 
                    lon=realprox.allts$lon, 
                    lat=realprox.allts$lat, 
                    names=realprox.allts$names, 
                    time=realprox.allts$time)  
  } 
  
  if (docu) {
    if (sixmonstatevector) { #allows for merging differnt temp. resolutions
      # change proxy array that only 2 ts instead of 12 monthly ts but 6 times as 
      # many variables, one for each month
      doc_t_mon$data <- t(doc_t_mon$data)
      tmp1 <- array(doc_t_mon$data,c(dim(doc_t_mon$data)[1] *  dim(doc_t_mon$data)[2]))
      tmp2 <- tmp1[((9*dim(doc_t_mon$data)[1]+1):(length(tmp1)[1] - 
                                               (3*dim(doc_t_mon$data)[1])))] # cut oct syr to sep eyr
      doc_t_mon$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(doc_t_mon$data)[2]
                                                /12)-1)*2),(((dim(doc_t_mon$data)[2]/12)-1)*2))) 
      # reconvert to 2 seasons per year
      # first all doc_t_mon for one month, than all prox for next month
      doc_t_mon$time <- seq(syr+1,eyr+0.5,0.5) 
      doc_t_mon$names <- rep(doc_t_mon$names,6)
      doc_t_mon$lon <- rep(doc_t_mon$lon,6)
      doc_t_mon$lat <- rep(doc_t_mon$lat,6)
    }
    ti=which(floor(doc_t_mon$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    doc_t_mon.allts=doc_t_mon
    doc_t_mon$data=doc_t_mon$data[,sts:ets]
    doc_t_mon$time=doc_t_mon$time[sts:ets]
    doc_t_mon$des=rep('mon',nrow(doc_t_mon$data))
  
#    # seasonal documentary data: attribiute JFM to winter season and later create fitting H operator
#     if (sixmonstatevector) { #allows for merging differnt temp. resolutions
#       doc_t_seas$data <- t(doc_t_seas$data)
#       tmp1 <- doc_t_seas$data[,((doc_t_seas$time-trunc(doc_t_seas$time)==0.125) | 
#                                   (doc_t_seas$time-trunc(doc_t_seas$time)==0.625))]
#       tmp2 <- array(tmp1,c(dim(tmp1)[1] *  dim(tmp1)[2]))
#       tmp3 <- tmp2[(2*dim(tmp1)[1]+1):(length(tmp2)[1])] # cut oct syr to sep eyr
#       doc_t_seas$data <- array(tmp3,c(dim(tmp3)[1]/(((dim(tmp1)[2]/2)-1)*2),(((dim(tmp1)[2]/2)-1)*2)))
#       doc_t_seas$time <- seq(syr+1,eyr+0.5,0.5) 
# #      doc_t_seas$names <- rep(doc_t_seas$names,2)
# #      doc_t_seas$lon <- rep(doc_t_seas$lon,2)
# #      doc_t_seas$lat <- rep(doc_t_seas$lat,2)
#     }
#     ti=which(floor(doc_t_seas$time)==cyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_seas.allts=doc_t_seas
#     doc_t_seas$data=doc_t_seas$data[,sts:ets]
#     doc_t_seas$time=doc_t_seas$time[sts:ets]
#     doc_t_seas$des=rep('seas',nrow(doc_t_seas$data))
#     
#     # attribute JFMA documentary data (1 value per year) to winter season and later create fitting H operator
#     ti=which(floor(doc_t_JFMA$time)==cyr)
#     doc_t_JFMA.allts=doc_t_JFMA
#     if (sixmonstatevector) {
#       doc_t_JFMA$data=t(matrix(c(doc_t_JFMA$data[ti],NA)))
#       doc_t_JFMA.allts$data=t(matrix(as.vector(rbind(doc_t_JFMA.allts$data[2:length(doc_t_JFMA.allts$data)],rep(NA,(length(doc_t_JFMA.allts$data)-1))))))
#     } else {
#       doc_t_JFMA$data=doc_t_JFMA$data[ti]
#       print("ACHTUNG: documentary proxies only included with option sixmonstatevector so far!!!")
#     }
#     doc_t_JFMA$time=doc_t_JFMA$time[ti]-0.5
#     doc_t_JFMA$des=rep('JFMA',nrow(doc_t_JFMA$data))
#     
#     # attribute AMJJA documentary data (1 value per year) to summer season and later create fitting H operator
#     ti=which(floor(doc_t_AMJJA$time)==cyr)
#     doc_t_AMJJA.allts=doc_t_AMJJA
#     if (sixmonstatevector) {
#       doc_t_AMJJA$data=t(matrix(c(NA,doc_t_AMJJA$data[ti])))
#       doc_t_AMJJA.allts$data=t(matrix(as.vector(rbind(doc_t_AMJJA.allts$data[2:length(doc_t_AMJJA.allts$data)],rep(NA,(length(doc_t_AMJJA.allts$data)-1))))))
#     } else {
#       doc_t_JFMA$data=doc_t_AMJJA$data[ti]
#       print("ACHTUNG: documentary proxies only included with option sixmonstatevector so far!!!")
#     }
#     doc_t_AMJJA$time=doc_t_AMJJA$time[ti]
#     doc_t_AMJJA$des=rep('AMJJA',nrow(doc_t_AMJJA$data))

    if ((!instrumental) & (!real_proxies)) {
      proxies <- doc_t_mon
#       proxies <- list(lon=c(doc_t_mon$lon,doc_t_seas$lon,doc_t_JFMA$lon,doc_t_AMJJA$lon), 
#                         lat=c(doc_t_mon$lat,doc_t_seas$lat,doc_t_JFMA$lat,doc_t_AMJJA$lat),
#                         data=rbind(doc_t_mon$data,doc_t_seas$data,doc_t_JFMA$data,doc_t_AMJJA$data), 
#                         names=c(doc_t_mon$names,doc_t_seas$names,doc_t_JFMA$names,doc_t_AMJJA$names),
#                         des=c(doc_t_mon$des,doc_t_seas$des,doc_t_JFMA$des,doc_t_AMJJA$des),
#                         time=doc_t_mon$time)
      proxies.allts <- doc_t_mon.allts
#       proxies.allts <- list(lon=c(doc_t_mon.allts$lon,doc_t_seas.allts$lon,doc_t_JFMA.allts$lon,doc_t_AMJJA.allts$lon), 
#                            lat=c(doc_t_mon.allts$lat,doc_t_seas.allts$lat,doc_t_JFMA.allts$lat,doc_t_AMJJA.allts$lat),
#                            data=rbind(doc_t_mon.allts$data,doc_t_seas.allts$data,doc_t_JFMA.allts$data,doc_t_AMJJA.allts$data), 
#                            names=c(doc_t_mon.allts$names,doc_t_seas.allts$names,doc_t_JFMA.allts$names,doc_t_AMJJA.allts$names),
#                            des=c(doc_t_mon.allts$des,doc_t_seas.allts$des,doc_t_JFMA.allts$des,doc_t_AMJJA.allts$des),
#                            time=doc_t_mon.allts$time)        
    } else {
      docall <- doc_t_mon
#       docall <- list(lon=c(doc_t_mon$lon,doc_t_seas$lon,doc_t_JFMA$lon,doc_t_AMJJA$lon), 
#                       lat=c(doc_t_mon$lat,doc_t_seas$lat,doc_t_JFMA$lat,doc_t_AMJJA$lat),
#                       data=rbind(doc_t_mon$data,doc_t_seas$data,doc_t_JFMA$data,doc_t_AMJJA$data), 
#                       names=c(doc_t_mon$names,doc_t_seas$names,doc_t_JFMA$names,doc_t_AMJJA$names),
#                       time=doc_t_mon$time) 
      docall.allts <- doc_t_mon.allts
    }
  }
  
  if (instrumental) {
    if ((ghcn_temp) & (yuri_temp)) {
      inst_t<-list(data=cbind(ghcn$data,inst_t$data), lon=c(ghcn$lon,inst_t$lon), 
                   lat=c(ghcn$lat,inst_t$lat), names=c(ghcn$names,inst_t$names), 
                   height=c(ghcn$height,inst_t$height), time=ghcn$time)
    }
    if ((ghcn_temp) & (!yuri_temp)) {  
      inst_t <- ghcn
    }
    if (ghcn_prec) {
      inst_p <- ghcn_precip
    }
    
        
    if (first_prox_per_grid) {
      res <- firstproxres
#      newgrid <- valiall
#      newgrid$lon <- valiall$lon[seq(1,length(valiall$lon),res)]
#      newgrid$lat <- valiall$lat[seq(1,length(valiall$lon),res)]
      newgrid <- echam
      newgrid$data <- NULL
      newgrid$ensmean <- NULL
      newgrid$lon <- echam$lon[seq(1,length(echam$lon),res)]
      newgrid$lat <- echam$lat[seq(1,length(echam$lat),res)]
      
      
      tmp_t <- inst_t
      tmp_t$data <- t(inst_t$data)
      Hproxy <- compute_H(tmp_t, newgrid, threshold=(700*res))  
      p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
      inst_t$lon <- inst_t$lon[p.i]
      inst_t$lat <- inst_t$lat[p.i]
      inst_t$data <- inst_t$data[,p.i]
      inst_t$names <- inst_t$names[p.i] 

      if (ghcn_prec) {
        tmp_precip <- inst_p
        tmp_precip$data <- t(inst_p$data)
        Hproxy <- compute_H(tmp_precip, newgrid, threshold=(700*res))  
        p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
        inst_p$lon <- inst_p$lon[p.i]
        inst_p$lat <- inst_p$lat[p.i]
        inst_p$data <- inst_p$data[,p.i]
        inst_p$names <- inst_p$names[p.i] 
      }

      if (yuri_slp) {  
        tmp_slp <- inst_slp
        tmp_slp$data <- t(inst_slp$data)
        Hproxy <- compute_H(tmp_slp, newgrid, threshold=(700*res))  
        p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
        inst_slp$lon <- inst_slp$lon[p.i]
        inst_slp$lat <- inst_slp$lat[p.i]
        inst_slp$data <- inst_slp$data[,p.i]
        inst_slp$names <- inst_slp$names[p.i] 
      }

      if (real_proxies) {
        print("DON'T USE FIRST_PROX_PER_GRID IF YOU WANT TO INCLUDE REAL PROXY DATA AND 
              INSTRUMENTALS AT THE SAME TIME")
      }  
    }        
        
    if (avg_prox_per_grid) {
      # average data in same grid box and count number of avg. series as error estimate
      # separate temp, precip, slp
      # makes no sense for realprox: proxies would need to be calibrated before building 
      # regession model
      if (ghcn_prec) {
        varlist <- c("inst_t","inst_p","inst_slp")
      } else {
        varlist <- c("inst_t","inst_slp") 
      }
      for (varname in varlist) {
        var <- get(varname)
        dlist=NA
        for(i in 1:length(var$lon)){
          plon <- var$lon[i]
          plat <- var$lat[i]
          clon <- echam$lon[!is.na(echam$lon)]
          clat <- echam$lat[!is.na(echam$lat)]
          k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
          # +0.001 to avoid find 2 locations with same distance
          l=which(abs(clat-plat)==min(abs(clat-plat))) 
          m=k[which(match(k,l)>0)]
          if (length(m) > 0) {
            dlist[i]=m
          } else {
            dlist[i]=NA
          }
        }
        var.avg=rep(NA,12) 
        var.numavg=NA
        var.lon=NA
        var.lat=NA
        var.names=NA

        for(i in unique(dlist)[(!is.na(unique(dlist)))]){
          no=which(dlist==i)
          mask=as.logical(dlist==i)
          mask[is.na(mask)]=FALSE
          if (length(var$data[1,mask])>1) {
            var.avg=rbind(var.avg,apply(var$data[,mask],1,mean,na.rm=T))
            var.numavg=rbind(var.numavg,(dim(var$data[,mask])[2]))
            var.lon=rbind(var.lon,mean(var$lon[mask],na.rm=T))
            var.lat=rbind(var.lat,mean(var$lat[mask],na.rm=T))
            var.names=rbind(var.names,var$names[mask][1])
          } else {
            var.avg=rbind(var.avg,var$data[,mask])
            var.numavg=rbind(var.numavg,1)
            var.lon=rbind(var.lon,var$lon[mask])
            var.lat=rbind(var.lat,var$lat[mask])
            var.names=rbind(var.names,var$names[mask])
          }
          var.avg[var.avg=='NaN']=NA
        }
        var$lon <- var.lon[2:length(var.lon)]
        var$lat <- var.lat[2:length(var.lat)]
        var$data <- t(var.avg[2:dim(var.avg)[1],])
        var$numavg <- var.numavg[2:length(var.numavg)]
        var$names <- var.names[2:length(var.names)]
        var$ensmean <- NULL
        if (varname == "inst_t") {
          inst_t = var
        }
        if (varname == "inst_p") {
          inst_p = var
        }    
        if (varname == "inst_slp") {
          inst_slp = var
        }
      }
    } 
    
    if (ghcn_prec) {
      inst<-list(data=t(cbind(inst_t$data,inst_p$data,inst_slp$data)), 
                 lon=c(inst_t$lon,inst_p$lon,inst_slp$lon), 
                 lat=c(inst_t$lat,inst_p$lat,inst_slp$lat), 
                 numavg=c(inst_t$numavg,inst_p$numavg,inst_slp$numavg), 
                 names=c(inst_t$names,inst_p$names,inst_slp$names), 
                 height=c(inst_t$height,inst_p$height,inst_slp$height), time=inst_t$time)
    } else { # only t und slp
      inst<-list(data=t(cbind(inst_t$data,inst_slp$data)), lon=c(inst_t$lon,inst_slp$lon), 
                 lat=c(inst_t$lat,inst_slp$lat), names=c(inst_t$names,inst_slp$names), 
                 numavg=c(inst_t$numavg,inst_slp$numavg), 
                 height=c(inst_t$height,inst_slp$height), time=inst_t$time)
    }

    # mask proxy data if there is instrumental data in same grid box at that time
    if ((real_proxies) & (instrumental)) {
      # remove inst data that is all NA
      pos <- apply(!is.na(inst$data),1,any)
      inst$data = inst$data[pos,]
      if (length(which(pos))==1) {
        inst$data = t(matrix(inst$data))
      }
      inst$lon <- inst$lon[pos]
      inst$lat <- inst$lat[pos]
      inst$names <- inst$names[pos] 
      for(i in 1:length(realprox$lon)){
        plon <- realprox$lon[i]
        plat <- realprox$lat[i]
        clon <- inst$lon
        clat <- inst$lat
        dist <- compute_dist(plon, plat, clon, clat)
        if (min(dist)<700){
#          print(paste('set proxy',i,'to NA'))
          realprox$data[i,] = NA
          realprox$lon[i] <- NA
          realprox$lat[i] <- NA
          realprox$mr[i,] <- NA
          realprox$var_residu[i] <- NA
        }
      }
      pos <- apply(!is.na(realprox$data),1,any)
      realprox$data = realprox$data[pos,]
      realprox$lon <- realprox$lon[pos]
      realprox$lat <- realprox$lat[pos]
      realprox$mr <- realprox$mr[pos,]
      realprox$var_residu <- realprox$var_residu[pos]
      realprox$names <- realprox$names[pos] 
    }
    if ((docu) & (instrumental)) {
      ti <- which(floor(inst$time)==cyr)
      sts <- ti[1]
      ets <- ti[length(ti)]
      if (any(!is.na(inst$data[,sts:ets]))) {
        itmp <- apply(inst$data[,sts:ets],1,mean,na.rm=T)
        itmp[inst$names!="temp2"] <- NaN
        clon <- inst$lon[!is.nan(itmp)]
        clat <- inst$lat[!is.nan(itmp)]
        for(i in 1:length(docall$lon)){
          plon <- docall$lon[i]
          plat <- docall$lat[i]
          if ((!is.na(plon)) & (!is.na(plat))) {    
            dist <- compute_dist(plon, plat, clon, clat)
            if (min(dist)<700){
#            print(paste('set doc',i,'to NA'))
              docall$data[i,] = NA
              docall$lon[i] <- NA
              docall$lat[i] <- NA
            }
          }
#           pos <- apply(!is.na(docall$data),1,any)
#           docall$data = docall$data[pos,]
#           docall$lon <- docall$lon[pos]
#           docall$lat <- docall$lat[pos]
#           docall$height <- docall$height[pos]
#           docall$season <- docall$season[pos]
#           docall$des <- docall$des[pos]
#           docall$names <- docall$names[pos]
        }
      }       
    }
    
    if (reduced_proxies) {
      every <- 12 
      if (real_proxies) {
        redpos <-seq(1,length(realprox$lon),every)
        realprox$lon <- realprox$lon[redpos]
        realprox$lat <- realprox$lat[redpos]
        realprox$data <- realprox$data[redpos,]
        realprox$names <- realprox$names[redpos]
        realprox.allts$lon <- realprox.allts$lon[redpos]
        realprox.allts$lat <- realprox.allts$lat[redpos]
        realprox.allts$data <- realprox.allts$data[redpos,]
        realprox.allts$names <- realprox.allts$names[redpos]
      }
      if (instrumental) {
        redpos <-seq(1,length(inst$lon),every)
        inst$lon <- inst$lon[redpos]
        inst$lat <- inst$lat[redpos]
        inst$data <- inst$data[redpos,]
        inst$names <- inst$names[redpos]
      }
    }
    # add data source
    if (instrumental) {
    inst$sour <- rep('inst',length(inst$lon))    }
    if (docu) {
    docall$sour <- rep('doc',length(docall$lon))
    docall.allts$sour <- rep('doc',length(docall$lon))    }
    if (real_proxies) {
      realprox$sour <- rep('prox',length(realprox$lon))
      realprox.allts$sour <- rep('prox',length(realprox$lon))    }

    if (sixmonstatevector) { 
      # change proxy array that only 2 ts instead of 12 monthly ts but 6 times as 
      # many variables, one for each month
      tmp1 <- array(inst$data,c(dim(inst$data)[1] *  dim(inst$data)[2]))
      tmp2 <- tmp1[((9*dim(inst$data)[1]+1):(length(tmp1)[1] - 
                    (3*dim(inst$data)[1])))] # cut oct syr to sep eyr
      inst$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(inst$data)[2]
                                 /12)-1)*2),(((dim(inst$data)[2]/12)-1)*2))) 
      # reconvert to 2 seasons per year
      # first all inst for one month, than all prox for next month
      inst$time <- seq(syr+1,eyr+0.5,0.5) 
      inst$names <- rep(inst$names,6)
      inst$sour <- rep(inst$sour,6)
      inst$lon <- rep(inst$lon,6)
      inst$lat <- rep(inst$lat,6)
      if (avg_prox_per_grid) {
        inst$numavg <- rep(inst$numavg,6)
      }
    }
    ti=which(floor(inst$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    inst.allts=inst
    inst$data=inst$data[,sts:ets]
    inst$time=inst$time[sts:ets]
    
    
    if (!real_proxies) {        
      proxies <- inst
      proxies.allts <- inst.allts
    }

    if (real_proxies) {
      tmpmr <- matrix(NA,nrow=length(inst$lon),ncol=ncol(realprox$mr))
      tmpres <- rep(NA,length(realprox$var_residu))
      if (avg_prox_per_grid) {
        tmpnum1 <- inst$numavg
        tmpnum2 <- rep(1,length(realprox$lon))
      } else {
        tmpnum1 <- rep(1,length(inst$lon))
        tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
      }
      if (real_proxies & instrumental) {
        proxies<-list(data=rbind(inst$data,realprox$data), lon=c(inst$lon,realprox$lon), 
                    lat=c(inst$lat,realprox$lat), names=c(inst$names,realprox$names), 
                    sour=c(inst$sour,realprox$sour), 
                    height=c(inst$height,realprox$elevation), time=inst$time,
                    mr=rbind(tmpmr,realprox$mr), var_residu=c(tmpres,realprox$var_residu),
                    numavg=c(tmpnum1,tmpnum2))
        proxies.allts<-list(data=rbind(inst.allts$data,realprox.allts$data), 
                          lon=c(inst.allts$lon,realprox.allts$lon), 
                          lat=c(inst.allts$lat,realprox.allts$lat), 
                          names=c(inst.allts$names,realprox.allts$names), 
                          sour=c(inst.allts$sour,realprox.allts$sour), 
                          time=inst.allts$time)  
                          # WHAT TO DO WITH $MR AND $VAR_RESIDU and INST.NUMAVG???      
      }
    }

    if (docu) {
      tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
      tmpres <- rep(NA,length(proxies$var_residu))
      tmpelev <- rep(NA,length(proxies$height))
      tmpnum3 <- rep(1,length(docall$lon))
      proxies<-list(data=rbind(proxies$data,docall$data), lon=c(proxies$lon,docall$lon), 
                    lat=c(proxies$lat,docall$lat), names=c(proxies$names,docall$names), 
                    sour=c(proxies$sour,docall$sour), 
                    height=c(proxies$height,tmpelev), time=proxies$time,
                    mr=rbind(proxies$mr,tmpmr), var_residu=c(realprox$var_residu,tmpres),
                    numavg=c(proxies$numavg,tmpnum3))
      proxies.allts<-list(data=rbind(proxies.allts$data,docall.allts$data), 
                    lon=c(proxies.allts$lon,docall.allts$lon), 
                    lat=c(proxies.allts$lat,docall.allts$lat), 
                    names=c(proxies.allts$names,docall.allts$names), 
                    sour=c(proxies.allts$sour,docall.allts$sour), 
                    time=proxies.allts$time)  
    }
  }

  if (!instrumental & real_proxies & docu) {
    docall$sour <- rep('doc',length(docall$lon))
    docall.allts$sour <- rep('doc',length(docall$lon))    
    realprox$sour <- rep('prox',length(realprox$lon))
    realprox.allts$sour <- rep('prox',length(realprox$lon)) 
    tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
    tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
    tmpres <- rep(NA,length(proxies$var_residu))
    tmpelev <- rep(NA,length(proxies$height))
    tmpnum3 <- rep(1,length(docall$lon))
    proxies<-list(data=rbind(realprox$data,docall$data), lon=c(realprox$lon,docall$lon), 
                  lat=c(realprox$lat,docall$lat), names=c(realprox$names,docall$names), 
                  sour=c(realprox$sour,docall$sour), 
                  height=c(realprox$height,tmpelev), time=realprox$time,
                  mr=rbind(realprox$mr,tmpmr), var_residu=c(realprox$var_residu,tmpres),
                  numavg=c(tmpnum2,tmpnum3))
    proxies.allts<-list(data=rbind(realprox.allts$data,docall.allts$data), 
                        lon=c(realprox.allts$lon,docall.allts$lon), 
                        lat=c(realprox.allts$lat,docall.allts$lat), 
                        names=c(realprox.allts$names,docall.allts$names), 
                        sour=c(realprox.allts$sour,docall.allts$sour), 
                        time=realprox.allts$time)  
  }

  if (real_proxies) {
    if (dim(realprox$data)[1]==0) { real_proxies=F }
  }
  
  calibrate <- proxies
  calibrate.allts <- proxies.allts
  if (sum(c(ncep_vali,cru_vali,recon_vali))>1) {
    print("WARNING: more than 1 validation data set selected!")
  }
  if (vali) {
    validate=valiall
    validate$ensmean=validate$data
  }
  
  print(paste('number of proxies/observations:',dim(calibrate$data)[1]))
  print("calc time preparing proxies")
  print(proc.time() - ptm1)




## convert to new data format
# ACHTUNG only 11 vars withOUT stream function included so far
  if (sixmonstatevector) {
    numvar <- length(unique(echam$names)) 
    if (ind_ECHAM) {
      numvar <- numvar-6
      echam$lon <-  rep(c(rep(echam$lon[!is.na(echam$lon)], numvar),rep(NA,length(echind$names))),6)
      echam$lat <- rep(c(rep(echam$lat[!is.na(echam$lat)], numvar),rep(NA,length(echind$names))),6)
      echam$names <- rep(echam$names, 6)
#    }
#     if (numvar > 3) {
#       echam$lon <- c(rep(echam$lon, (numvar*6)), unlist(echam[paste0('lon', 
#                      unique(echam$names)[-(1:numvar)])]))
#       echam$lat <- c(rep(echam$lat, (numvar*6)), unlist(echam[paste0('lat', 
#                      unique(echam$names)[-(1:numvar)])]))
    } else {
      echam$lon <-  c(rep(echam$lon, (numvar*6)))
      echam$lat <- c(rep(echam$lat, (numvar*6)))
      echam$names <- c(rep(echam$names, 6))
    }
#    echam$height <- c(rep(echam$height, (numvar*6)))
#    echam$height <- c(rep(echam$height, (numvar*6)), rep(NA, length=nrow(echam$data) - 
#                                                            length(echam$height)*(numvar*6)))
    echam <- echam[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
  }
  
  ## set up the cutoff distance for localisation
  lvec <- rep(l_dist_ind, length(unique(echam$names)))
  names(lvec) <- unique(echam$names)
  lvec['slp'] <- l_dist_slp
  lvec['precip'] <- l_dist_precip
  lvec['temp2'] <- l_dist_temp2 #l_dist
  lvec['gph500'] <- l_dist_gph500
  lvec['gph100'] <- l_dist_gph100
  lvec['u850'] <- l_dist_u850
  lvec['u200'] <- l_dist_u200
  lvec['v850'] <- l_dist_v850
  lvec['v200'] <- l_dist_v200
  lvec['omega500'] <- l_dist_omega500
  lvec['t500'] <- l_dist_t500
#  lvec <- c(lvec,l_dist_prox)
#  names(lvec) <- c(unique(echam$names),'prox')

# ------------------------------------------------------------------------------
# Compute H (forward operator)
# ------------------------------------------------------------------------------
  if ((real_proxies) & (!instrumental) & (!docu) & (!sixmonstatevector)) {
     Hcal <- Matrix(compute_H_proxy(realprox, echam, realprox$mr, threshold=700), 
                    sparse=T)
  }
#   } else if (instrumental & sixmonstatevector & !real_proxies & !docu) {
#     Hcal <- Matrix(compute_H_sixmonstatevector(calibrate, echam, threshold=700), 
#                    sparse=T)
#   } else if (docu & sixmonstatevector & !real_proxies & !instrumental) {
#     Hcal <- Matrix(compute_H_sixmonstatevector_docu(calibrate, echam, threshold=700), 
#                    sparse=T)
#   } else if (instrumental & sixmonstatevector & real_proxies & !docu) {  
#     Hcal1 <- compute_H_proxy(realprox, echam, realprox$mr, threshold=700)
#     Hcal2 <- compute_H_sixmonstatevector(inst, echam, threshold=700)
#     Hcal <- Matrix(rbind(Hcal2,Hcal1), sparse=T)
#     rm(Hcal1);rm(Hcal2)
#   } else 
  if (sixmonstatevector) {  
#    if (dim(echam$data)[1] > 10000) {
#      Hcal1 <- array(0,dim=c(dim(realprox$data)[1],dim(echam$data)[1]))
#      Hcal2 <- array(0,dim=c(dim(docall$data)[1],dim(echam$data)[1]))
#      Hcal3 <- ff(0,dim=c(dim(inst$data)[1],dim(echam$data)[1]))
#      Hcal <- ff(0,dim=c((dim(Hcal1)[1])+(dim(Hcal2)[1])+(dim(Hcal3)[1]),dim(echam$data)[1]))
#    } else {
#      Hcal1 <- array(0,dim=c(dim(realprox$data)[1],dim(echam$data)[1]))
#      Hcal2 <- array(0,dim=c(dim(docall$data)[1],dim(echam$data)[1]))
#      Hcal3 <- array(0,dim=c(dim(inst$data)[1],dim(echam$data)[1]))
#      Hcal <- array(0,dim=c((dim(Hcal1)[1])+(dim(Hcal2)[1])+(dim(Hcal3)[1]),dim(echam$data)[1]))
#    }
    # next 3 lines: solve problem of distance calc with NA in compute_dist function
    etmp <- echam
    etmp$lon[is.na(etmp$lon)] <- 0
    etmp$lat[is.na(etmp$lat)] <- -90
    if (instrumental) {
      Hcal1 <- array(NA,dim=c(dim(inst$data)[1],2))
      Hcal1 <- compute_Hi_Hredux_sixmonstatevector(inst, etmp, threshold=700)
    }
    if (real_proxies) {
      Hcal2 <- array(NA,dim=c(dim(realprox$data)[1],14))
      Hcal2 <- compute_Hi_Hredux_proxy(realprox, etmp, realprox$mr, threshold=700)
    }
    if (docu) {
      Hcal3 <- array(NA,dim=c(dim(docall$data)[1],2))
      Hcal3 <- compute_Hi_Hredux_sixmonstatevector(docall, etmp, threshold=700)
      Hcal3[Hcal3==0] <- NA
    }
    rm(etmp)
#    Hcal1[,] <- compute_H_proxy(realprox, echam, realprox$mr, threshold=700)
#    Hcal2[,] <- compute_H_sixmonstatevector(docall, echam, threshold=700)
#    Hcal3[,] <- compute_H_sixmonstatevector(inst, echam, threshold=700)
#    Hcal[1:(dim(Hcal3)[1]),] <- Hcal3[,]
#    Hcal[1:(dim(Hcal2)[1]),] <- Hcal2[,]
#    Hcal[1:(dim(Hcal1)[1]),] <- Hcal1[,]
    # Matrix(rbind(Hcal3,Hcal2,Hcal1), sparse=T)
#    rm(Hcal1);rm(Hcal2);rm(Hcal3)
#  } else {
#    print('ACHTUNG, compute_H has NOT been updated to Jonas new code')
#    Hcal <- Matrix(compute_H_temp_precip_slp(calibrate, echam),sparse=T)
#  }
    if (instrumental & !docu & sixmonstatevector & !real_proxies) {
    #    print('ACHTUNG, compute_H has NOT been updated to Jonas new code')
    #    H.i <- array(NA,c(nrow(Hcal),length(which(Hcal[1,]!=0))))
      H.i <- array(NA,c(nrow(Hcal1),1))
    } else if (instrumental & docu & sixmonstatevector & !real_proxies) {
#    print('ACHTUNG, compute_H has NOT been updated to Jonas new code')
#    H.i <- array(NA,c(nrow(Hcal),length(which(Hcal[1,]!=0))))
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal3)),1))
    } else if (!instrumental & docu & sixmonstatevector & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal2)+nrow(Hcal3)),7))
    } else if (instrumental & !docu & sixmonstatevector & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)),7))
    } else { 
#    H.i <- array(NA,c(nrow(Hcal),max(sapply(1:dim(Hcal)[1],function(x) length(which(Hcal[x,]!=0))))))
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)+nrow(Hcal3)),7))
    }
    Hredux <- H.i
    if (instrumental){
      H.i[1:nrow(Hcal1),1] <- Hcal1[,1] 
      Hredux[1:nrow(Hcal1),1] <- Hcal1[,2] 
    } 
    if (instrumental & real_proxies & docu) {
      H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
      Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(2,4,6,8,10,12,14)]   
      H.i[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
      Hredux[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
    }
    if (!instrumental & real_proxies & docu) {
      H.i[1:nrow(Hcal2),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
      Hredux[1:nrow(Hcal2),] <- Hcal2[,c(2,4,6,8,10,12,14)]   
      H.i[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
      Hredux[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
    }
    if (instrumental & real_proxies & !docu) {
      H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
      Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(2,4,6,8,10,12,14)] 
    }

    H.i[H.i==0] <- NA
    Hredux[Hredux==0] <- NA
  }
#   notna <- apply(!is.na(Hredux),1,any)
#   calibrate$data <- calibrate$data[notna,]
#   calibrate$mr <- calibrate$mr[notna,]
#   calibrate$names <- calibrate$names[notna]
#   calibrate$lon <- calibrate$lon[notna]
#   calibrate$lon <- calibrate$lat[notna]
#   calibrate$numavg <- calibrate$numavg[notna]
#   calibrate.allts$data <- calibrate.allts$data[notna,]
#   calibrate.allts$mr <- calibrate.allts$mr[notna,]
#   calibrate.allts$names <- calibrate.allts$names[notna]
#   calibrate.allts$lon <- calibrate.allts$lon[notna]
#   calibrate.allts$lon <- calibrate.allts$lat[notna]
#   calibrate.allts$numavg <- calibrate.allts$numavg[notna]
#   H.i <- H.i[notna,]
#   Hredux <- Hredux[notna,]
#   for (i in 1:nrow(Hcal)) {
# #    print(length(which(Hcal[i,]!=0)))
#     if ((length(which(Hcal[i,]!=0))) > 0){
#       if ((instrumental | docu) & sixmonstatevector & !real_proxies) {
#         H.i[i,] <- which(Hcal[i,]!=0)
#         Hredux[i,] <- Hcal[i,(Hcal[i,]!=0)]
# #      } else if (docu & sixmonstatevector & !real_proxies) {
# #        l <- 5-(length(which(Hcal[i,]!=0)))
# #        H.i[i,] <- c(which(Hcal[i,]!=0),rep(NA,l))
# #        Hredux[i,] <- c(Hcal[i,(Hcal[i,]!=0)],rep(NA,l))
#       } else {  
#         l <- 7-(length(which(Hcal[i,]!=0)))
#         H.i[i,] <- c(which(Hcal[i,]!=0),rep(NA,l))
#         Hredux[i,] <- c(Hcal[i,(Hcal[i,]!=0)],rep(NA,l))
#       }
# #      if ((length(which(Hcal[i,]!=0))) == 1){
# #        H.i[i,] <- c(which(Hcal[i,]!=0),rep(NA,6))
# #        Hredux[i,] <- c(Hcal[i,(Hcal[i,]!=0)],rep(NA,6))
# #      } else {
# #        H.i[i,] <- which(Hcal[i,]!=0)
# #        Hredux[i,] <- Hcal[i,(Hcal[i,]!=0)]
# #      }
#     }
#   }
#  Hval <- Matrix(compute_H_temp_precip(validate, echam),sparse=T)
#  Hcal <- compute_H(calibrate, echam)
#  Hval <- compute_H(validate, echam)
#
#   if (sixmonstatevector) {
#     tmpmat <- matrix(NA, nrow=nrow(Hcal1), ncol=ncol(Hcal1))
#     Hcalwin <- Matrix(rbind(Hcal2,tmpmat), sparse=T)
#     Hcalsum <- Hcal
#   }

  print("calc time H operator")
  print(proc.time() - ptm1)


# ------------------------------------------------------------------------------
# run analysis
# ------------------------------------------------------------------------------
  if ((instrumental) || (real_proxies) || (docu)) {
# R for perfect data would be
# R <- rep(0, nrow(Hcal))
# R is simply the squared random error assumed for instr. data because we assume 
# the error would be spatially uncorrelated
# R <- apply(array(disturb2, c(nrow(Hcal), 2, ncol(pdata)/2))**2, c(1,2), mean)
# set squared error R to 1 for 1degC measurement error
    if (sixmonstatevector) {
#      R1 = array(1,c(nrow(Hcal1), 1))
#      R2 = array(1,c(nrow(Hcal2), 1))
#      R1 = array(1,c(nrow(Hcalwin), 1))
#      R2 = array(1,c(nrow(Hcalsum), 1))
      if ((real_proxies) & ((instrumental) | (docu))) { #< 1970
#        Rcal <- c(temp2=0.5, precip=50, slp=10)[calibrate$names]
        Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
#        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
#        Rcal[calibrate$names=="prox"] <- realprox$var_residu
        Rcal[calibrate$names=="prox"] <- realprox$var_residu/2
        Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
      } else if ((real_proxies) & (!instrumental) & (!docu)) { # should not happen
#        R1[calibrate$names=='temp',]<-try(abs(calibrate$data[calibrate$names=='temp',1]
#                                             *0.0001))
        Rcal <- realprox$var_residu/2
        Rcal[which(is.na(Rcal))] <- 0
#        R1[,1] <- realprox$var_residu #/100 #to check why there is hardly any correction
#        R1[which(is.na(R1))] <- 0
#        R2[calibrate$names=='temp',]<-try(abs(calibrate$data[calibrate$names=='temp',2]
#                                             *0.0001))
#        R2[,1] <- realprox$var_residu
#        R2[which(is.na(R2))] <- 0 #/100 #to check why there is hardly any correction
      } else if (((instrumental) | (docu)) & (!real_proxies)) { # >1970
        # Jonas for pseudoproxies
#        Rcal <- c(temp2=0.5, precip=50, slp=20)[calibrate$names[1:(length(calibrate$names)/6)]]
        Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
        Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
#        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
#        R1 <- R2 <- Rcal
#        R = array(1,c(nrow(Hcal), s))
#      # set squared error R for precip measurements to 25% of data value
#      # R <- abs(calibrate$data*0.25)
#        R1[calibrate$names=='precip',]<-try(abs(calibrate$data[calibrate$names=='precip',1]
#                                            *0.5))
#        R1[calibrate$names=='slp',]<- 1
#        R1[which(is.na(R1))] <- 0
#        R2[calibrate$names=='precip',]<-try(abs(calibrate$data[calibrate$names=='precip',2]
#                                             *0.5))
#        R2[calibrate$names=='slp',]<- 1
#        R2[which(is.na(R2))] <- 0
      }  
#       if (fortran_ensrf) {calibrate$data[is.na(calibrate$data)]=9e36 #}
#         C1 <- calibrate
#         C1$data <- calibrate$data[,1,drop=F]
#         C1$time <- calibrate$time[1]
#         C2 <- calibrate
#         C2$data <- calibrate$data[,2,drop=F]
#         C2$time <- calibrate$time[2]
#         E1 <- echam
#         E1$data <- echam$data[,1,,drop=F]
#         E1$ensmean <- echam$ensmean[,1,drop=F]
#         E1$time <- echam$time[1]
#         E2 <- echam
#         E2$data <- echam$data[,2,,drop=F]
#         E2$ensmean <- echam$ensmean[,2,drop=F]
#         E2$time <- echam$time[2]
#         A1 <- EnSRF_new(E1, C1,  R=R1, Hcal=Hcal1, weights=d.weights_all)
#         A2 <- EnSRF_new(E2, C2,  R=R2, Hcal=Hcal2, weights=d.weights_all)
#         analysis <- A1
#         analysis$data <- abind(A1$data,A2$data,along=2)
#         analysis$ensmean <- abind(A1$ensmean,A2$ensmean,along=2)
#         analysis$time <- c(A1$time,A2$time)
#       } else {

## run the analysis from within R -- optimised for size of matrices rather than speed
        analysis <- echam
## take anomalies
        analysis$data <- echam$data - as.vector(echam$ensmean)
        nmonths <- 6
        ndim <- nrow(analysis$data)
        ntim <- ncol(analysis$data)
        nens <- dim(analysis$data)[3]
        nprox <- nrow(calibrate$data) #/nmmonths
        ndimold <- length(echam$lon) / nmonths
        itime <- rep(1:nmonths, each=ndimold)
        i=1
## compute loop over obs first to minimize repetition of computation (e.g. H and weights)
        for (j in 1:nprox){
#        for (j in which(!is.na(calibrate$data[,i]))) {
          if (j %% 100 == 0) {
            print(paste('Assimilating observation', j))
          }
          ## assume proxy location in state vector is known
          ## otherwise specify / compute h.i and reduced-size H here
          hisna <- is.na(H.i[j,])
          h.i <- H.i[j,!hisna,drop=F]
          H <-  t(as.matrix(Hredux[j,!hisna]))
        
#          print(paste("H:",H))
          if (!is.na(h.i[1])) {
            pos <- which(is.na(echam$lon))
            echam$lon[pos] <- 0
            echam$lat[pos] <- -90
            dist <- compute_dist_2d(echam$lon, echam$lat, echam$lon[h.i], echam$lat[h.i]) 
## weights are a matrix of ndim x nh (non-null elements in H, here different months)
#          wgt <- corr_function(dist, pmin(lvec[calibrate$names[j]], lvec[echam$names])) *
#                   outer(itime, 1:nmonths, function(x,y) exp(-0.5 * sqrt((x - y)**2)))
# NO temporal correlation for use of monthly instr. data and 1 col H operator          
            #wgt <- matrix(corr_function(dist,pmin(lvec[calibrate$names[j]],lvec[echam$names])))
            # wgt <- corr_function(dist,pmin(lvec[calibrate$names[j]],lvec[echam$names]))
            wgt <- corr_function(dist,outer(lvec[echam$names], lvec[echam$names[h.i]], pmin))
            wgt[which(echam$names %in% c("DIMI","z100","z300","PWC","HC","SJ")),] <- 1
## temporal correlation quickly drops to zero (~ 0.6 for 1-month lag, 
## ~0.4 for two months, etc.)
          } else {
            dist <- NA
            wgt <- NA
          }
          if ((!is.na(dist[1])) & (!is.na(wgt[1]))){
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
#              print(i)
                x2tmp <- analysis$data[,i,] # entire state vector at time step i, all members
                x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
                PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                K <- PH / HPHR
                Ktilde <- K / (1 + sqrt(Rcal[j]/HPHR))
                analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                          H %*% analysis$ensmean[h.i,i])
                analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% 
                                          analysis$data[h.i,i,]
              }
            }
          }  
#        print(analysis$data[c(1,100,1000,10000),1,1])
#        print(analysis$ensmean[c(1,100,1000,10000),1])
#        print(analysis$data[which(analysis$names=="HC"),1,1]
        }    
        ## add ensemble mean analysis back in
        analysis$data <- analysis$data + as.vector(analysis$ensmean)
        if (anomaly_assim){
          analysis.anom <- analysis
          analysis.abs <- analysis
          analysis.abs$data <- analysis$data + as.vector(echam.clim$data)
          analysis.abs$ensmean <- analysis$ensmean + as.vector(echam.clim$ensmean)
          echam.anom <- echam
          echam.abs <- echam
          echam.abs$data <- echam$data + as.vector(echam.clim$data)
          echam.abs$ensmean <- echam$ensmean + as.vector(echam.clim$ensmean)
        }
        if (vali){
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate,file=paste0('../data/analysis/analysis_',cyr,'.Rdata'))
        } else {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,calibrate,file=paste0('../data/analysis/analysis_',cyr,'.Rdata'))
        }
#        print(analysis$data[1,i,1])
#        if (is.na(analysis$data[1,i,1])) stop("analysis became NA")
#        if (analysis$data[1,i,1] > 10) stop("analysis became NA")
#      }
    } else {      
      R = array(1,c(nrow(Hcal), s))
      R[calibrate$names=='precip',]<-try(abs(calibrate$data[calibrate$names=='precip',]
                                            *0.5))  
      R[calibrate$names=='slp',]<-try(abs(calibrate$data[calibrate$names=='slp',]
                                            *0.1))  
      R[which(is.na(R))] <- 0
      if (fortran_ensrf) {calibrate$data[is.na(calibrate$data)]=9e36}  
      analysis <- EnSRF_new(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
#      analysis <- EnSRF4(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
      if (vali) {
        if (every2grid){
          save(analysis,echam,validate,calibrate,calibrate.allts,
               file=paste0('../data/analysis/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis,echam,validate,calibrate,calibrate.allts,
             file=paste0('../data/analysis/analysis_',cyr,'.Rdata'))
        }
      } else {
        if (every2grid){
          save(analysis,echam,calibrate,calibrate.allts,
               file=paste0('../data/analysis/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {  
          save(analysis,echam,calibrate,calibrate.allts,
             file=paste0('../data/analysis/analysis_',cyr,'.Rdata'))
        }
      }
    }
  }
print("calc time for a year")
print(proc.time() - ptm1)
} # end of yearly calc loop




# ---------------------------------------------------------------------------------------------
# Save output to file for display with EnSRF.Rnw
# ------------------------------------------------------------------------------
#save.image(paste("../data/EnSRF_",syr,"-",eyr,".Rdata",sep=""))

warnings()
#quit(save='no')
