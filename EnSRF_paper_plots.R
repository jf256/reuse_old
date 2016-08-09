# THINGS TO DO!

# check if previously calc seasonal avg can be used to speed it up, otherwise do seasonal means and save them to file!

# READ EXTRA VALIDATION DATA AND MAKE PLOTS

# NOTES:
# NO OCEAN INDICES IN LAND_ONLY VERSION

rm(list=ls())

machine="macbook" #"climcal3" # "climpa12" # 
if (machine=="macbook") {
  datadir="/Volumes/DATA/climdata/"
} else {
  datadir="/scratch/joerg/climdata/"
}

save_prepplot=F
  # load already transformed data to seasonal resolution, embedded in save_prepplot
  load_precalc=T 
mergetime=F
load_prepplot=T
  # monthly or seasonal analysis
  # '../data/prepplot_month/' or '../data/prepplot_season'
  prepplotpath=paste0(datadir,"ENSRF_analysis/prepplot_v3_seasonal/")
timeseriesplots=T
drought_1790=F
plots1816=F
multiann=F # check multiannual drought
warmcold_decades=F # warm, cold periods global
warmcold_us=F # warm, cold periods US
warmcold_eu=F # warm, cold periods EU
warmcold_nh=F # warm, cold periods NH
dryhumid_us=F # dry, humid periods US
dryhumid_eu=F # dry, humid periods EU
landusebugbias=F # check bias due to land use bug, see also Marco's analysis

# Olivia Nov. 2015: 
# add wind vectors to Z500!
# show assimilation gain over model
# Z500 and temp. anom. at same place = radiation is cause
# Z500 and temp. anom. shifted = advection is cause
# precip anom. in % instead of mm, otherwise just tropics visible

# if (timeseriesplots) {
#   syr=1603     
#   eyr=2004
# } else 
if (drought_1790) {
  syr <- 1790
  eyr <- 1799  
} else if (plots1816) {
  syr <- 1816
  eyr <- 1817  
} else {
  syr=1603     
  eyr=2004
}

setwd('~/unibe/projects/EnSRF/src')
source('EnSRF_functions.R')
#source('r_functions_joerg.r')

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

nmem=30
first_prox_per_grid=F  # first proxy per echam grid box ATTENTION: only this 
# or second next option (avg_prox_per_grid) can be TRUE
firstproxres=10       # grid resolution for instr. stations (5 = echamgrid/5)
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

# choose data (only 1 of the following options must be TRUE!!!)
if (eyr > 1659) {
  instrumental=T        # all instrumental stations
} else {
  instrumental=F
}
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
#print(paste("instr:",instrumental, "proxies:",real_proxies, "documentary:",docu))

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
twcr_vali=F            # 20CR reanalysis data for validation
ncep_vali=F            # NCEP/NCAR reanalysis data for validation





if (timeseriesplots) { # start timeseriesplots

if (save_prepplot) {
  ptm1 <- proc.time()
  for (cyr in syr:eyr) {
    if ((cyr < 1751) | (cyr > 1979)) {
      vali=F                 # switch off prepplot if no vali data selected
    } else {
      vali=T
    }
    if ((cyr < 1901) & (cyr > 1750)) {
      recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
    } else {
      recon_vali=F
    }
    
    print(cyr)
    if (load_precalc) {
      # NEW load already calculated seasonal data
      load(file=paste0('/Volumes/data/climdata/ENSRF_analysis/prepplot_v3_seasonal/analysis_',cyr,'.Rdata'))
#      load(file=paste0('../data/prepplot_season/analysis_',cyr,'.Rdata'))
      echam.abs <- echam
      analysis.abs <- analysis
      bemlist <- read.table(file='../data/bem/bem.txt',header=F)
      bem <- bemlist[which(bemlist[,1]==cyr),2]
      echam.abs$bem <- echam.abs$data[,,bem]
      echam.anom$bem <- echam.anom$data[,,bem]
      analysis.abs$bem <- analysis.abs$data[,,bem]
      analysis.anom$bem <- analysis.anom$data[,,bem]
      ech_ind$bem <- ech_ind$data[,,bem]
      ana_ind$bem <- ana_ind$data[,,bem]
    } else {   
      load(file=paste0('../data/analysis/analysis_',cyr,'.Rdata'))
      bemlist <- read.table(file='../data/bem/bem.txt',header=F)
      bem <- bemlist[which(bemlist[,1]==cyr),2]
      echam.abs$bem <- echam.abs$data[,,bem]
      echam.anom$bem <- echam.anom$data[,,bem]
      analysis.abs$bem <- analysis.abs$data[,,bem]
      analysis.anom$bem <- analysis.anom$data[,,bem]
      #  if (anomaly_assim){
      #    echam <- echam.abs
      #    analysis <- analysis.abs
      #    echam.abs <- NULL
      #    analysis.abs <- NULL   
      #   }
      
      if (ind_ECHAM) {
        # extract the indices from output 
        # 1. Analysis indices
        ech_ind <- echam.anom
        ech_ind$data <- echam.anom$data[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ"),,]
        ech_ind$ensmean <- echam.anom$ensmean[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ"),]
        ech_ind$bem <- echam.anom$bem[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ"),]
        ech_ind$names <- echam.anom$names[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ")]
        ech_ind$lon <- echam.anom$lon[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ")]
        ech_ind$lat <- echam.anom$lat[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                             echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                             echam.anom$names=="SJ")]
        ana_ind <- analysis.anom
        ana_ind$data <- analysis.anom$data[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ"),,]
        ana_ind$ensmean <- analysis.anom$ensmean[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ"),]
        ana_ind$bem <- analysis.anom$bem[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ"),]
        ana_ind$names <- analysis.anom$names[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ")]
        ana_ind$lon <- analysis.anom$lon[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ")]
        ana_ind$lat <- analysis.anom$lat[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                             analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                             analysis.anom$names=="SJ")]
      }
      
      
      # remove precalc indices from echam
      #  echam <- echam.anom
      echam.anom$data <- echam.anom$data[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ"),,]
      echam.anom$ensmean <- echam.anom$ensmean[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ"),]
      echam.anom$bem <- echam.anom$bem[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ"),]
      echam.anom$names <- echam.anom$names[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ")]
      echam.anom$lon <- echam.anom$lon[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ")]
      echam.anom$lat <- echam.anom$lat[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                             echam.anom$names!="z300" & echam.anom$names!="PWC" & echam.anom$names!="HC" & 
                             echam.anom$names!="SJ")]
      #  analysis <- analysis.anom
      analysis.anom$data <- analysis.anom$data[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ"),,]
      analysis.anom$ensmean <- analysis.anom$ensmean[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ"),]
      analysis.anom$bem <- analysis.anom$bem[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ"),]
      analysis.anom$names <- analysis.anom$names[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ")]
      analysis.anom$lon <- analysis.anom$lon[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ")]
      analysis.anom$lat <- analysis.anom$lat[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                             analysis.anom$names!="SJ")]
      echam.abs$data <- echam.abs$data[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ"),,]
      echam.abs$ensmean <- echam.abs$ensmean[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ"),]
      echam.abs$bem <- echam.abs$bem[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ"),]
      echam.abs$names <- echam.abs$names[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ")]
      echam.abs$lon <- echam.abs$lon[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ")]
      echam.abs$lat <- echam.abs$lat[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
                             echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
                             echam.abs$names!="SJ")]
      analysis.abs$data <- analysis.abs$data[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ"),,]
      analysis.abs$ensmean <- analysis.abs$ensmean[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ"),]
      analysis.abs$bem <- analysis.abs$bem[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ"),]
      analysis.abs$names <- analysis.abs$names[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ")]
      analysis.abs$lon <- analysis.abs$lon[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ")]
      analysis.abs$lat <- analysis.abs$lat[which(analysis.abs$names!="DIMI" & analysis.abs$names!="z100" & 
                             analysis.abs$names!="z300" & analysis.abs$names!="PWC" & analysis.abs$names!="HC" & 
                             analysis.abs$names!="SJ")]
      if (vali) {  
        vali_ind <- validate
        vali_ind$data <- validate$data[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" |
                             validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                             validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj"),]
        vali_ind$names <- validate$names[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" |
                             validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                             validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
        vali_ind$lon <- validate$lon[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" | 
                             validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                             validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
        vali_ind$lat <- validate$lat[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" | 
                             validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                             validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
        validate$data <- validate$data[which(validate$names=="temp2" | validate$names=="precip" |
                             validate$names=="slp"),]
        validate$names <- validate$names[which(validate$names=="temp2" | validate$names=="precip" |
                             validate$names=="slp")]
        validate$lon <- validate$lon[which(validate$names=="temp2" | validate$names=="precip" |
                             validate$names=="slp")]
        validate$lat <- validate$lat[which(validate$names=="temp2" | validate$names=="precip" |
                             validate$names=="slp")]
      }
    }
    
    if (land_only) {
      xlim=c(-180,180)
      ylim=c(-90,90)
      nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
      lon <- nc$dim$lon$vals
      lat <- nc$dim$lat$vals
      lon[lon > 180] <- lon[lon > 180] - 360
      loi <- which(lon >= xlim[1] & lon <= xlim[2])
      lai <- which(lat >= ylim[1] & lat <= ylim[2])
      # mulc for reading each 2rd grid cell to avoid memory problems
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
      lsm <- ncvar_get(nc)[loi, lai]
      nc_close(nc)
      lonlist=round(rep(lon[loi], length(lai))[lsm > 0.5],digits=2)
      latlist=round(rep(lat[lai], each=length(loi))[lsm > 0.5],digits=2)
      landll <- paste(lonlist,latlist)
      
      e.abs.ll <- paste(round(echam.abs$lon,digits=2),round(echam.abs$lat,digits=2))
      landpos <- match(landll,e.abs.ll)
      landpos_init <- landpos
      if (!load_precalc) {  
        ngrid <- dim(echam.abs$data)[1]/(nseas/s)/length(unique(echam.abs$names))
        nrep <-  (nseas/s)*length(unique(echam.abs$names))
      } else {
        # NEW if prepplot data are loaded
        ngrid <- dim(echam.abs$data)[1]/length(unique(echam.abs$names))
        nrep <-  length(unique(echam.abs$names))
      }
      for (i in 1:(nrep-1)) {
        landpos <- c(landpos,(landpos_init+ngrid*i))
      }
      echam.anom$data <- echam.anom$data[landpos,,] 
      echam.anom$ensmean <- echam.anom$ensmean[landpos,] 
      echam.anom$bem <- echam.anom$bem[landpos,] 
      echam.anom$lon <- echam.anom$lon[landpos] 
      echam.anom$lat <- echam.anom$lat[landpos] 
      echam.anom$names <- echam.anom$names[landpos] 
      echam.abs$data <- echam.abs$data[landpos,,] 
      echam.abs$ensmean <- echam.abs$ensmean[landpos,] 
      echam.abs$bem <- echam.abs$bem[landpos,] 
      echam.abs$lon <- echam.abs$lon[landpos] 
      echam.abs$lat <- echam.abs$lat[landpos] 
      echam.abs$names <- echam.abs$names[landpos] 
      analysis.anom$data <- analysis.anom$data[landpos,,] 
      analysis.anom$ensmean <- analysis.anom$ensmean[landpos,] 
      analysis.anom$bem <- analysis.anom$bem[landpos,] 
      analysis.anom$lon <- analysis.anom$lon[landpos] 
      analysis.anom$lat <- analysis.anom$lat[landpos] 
      analysis.anom$names <- analysis.anom$names[landpos] 
      analysis.abs$data <- analysis.abs$data[landpos,,] 
      analysis.abs$ensmean <- analysis.abs$ensmean[landpos,] 
      analysis.abs$bem <- analysis.abs$bem[landpos,] 
      analysis.abs$lon <- analysis.abs$lon[landpos] 
      analysis.abs$lat <- analysis.abs$lat[landpos] 
      analysis.abs$names <- analysis.abs$names[landpos] 
    }
    
    
    if (!load_precalc) {    
      # change to seasonal resolution
      if (sixmonstatevector)
        #  echam.anom.mon <- echam.anom        
        #  analysis.anom.mon <- analysis.anom
        #  echam.abs.mon <- echam.abs        
        #  analysis.abs.mon <- analysis.abs
        #  validate.mon <- validate
        #  if (landonly) {
        #    for (p in 1:length(which(!is.na(echam$lat))))
        #    vpos <- which(!is.na(validate$data[,1]))
        #    Hind2 <- Hind[,vpos]
        #  }
        
      etmp <- array(echam.abs$data,c((dim(echam.abs$data)[1]/6),6,dim(echam.abs$data)[2],
                                     dim(echam.abs$data)[3]))
      echam.abs$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
      for (ensmem in 1:(dim(etmp)[4])) {
        print(paste('ECHAM member',ensmem))
        echam.abs$data[,,ensmem] <- apply(etmp[,,,ensmem],c(1,3),mean)
      }
      #      echam.abs$data <- apply(etmp,c(1,3,4),mean)
      etmp2 <- array(echam.abs$ensmean,c((dim(echam.abs$ensmean)[1]/6),6,dim(echam.abs$ensmean)[2]))
      echam.abs$ensmean <- apply(etmp2,c(1,3),mean)
      etmp3 <- array(echam.abs$bem,c((dim(echam.abs$bem)[1]/6),6,dim(echam.abs$bem)[2]))
      echam.abs$bem <- apply(etmp3,c(1,3),mean)
      #    echam.abs$names <- echam.abs$names[1:dim(echam.abs$data)[1]/length(unique(echam.abs$names))]
      echam.abs$names <- rep(unique(echam.abs$names),each=dim(echam.abs$data)[1]/
                               length(unique(echam.abs$names)))  
      echam.abs$lon <- echam.abs$lon[1:(dim(echam.abs$data)[1])] #/length(unique(echam.abs$names)))]
      echam.abs$lat <- echam.abs$lat[1:(dim(echam.abs$data)[1])] #/length(unique(echam$names)))]
      etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
                                   dim(ech_ind$data)[3]))  
      ech_ind$data <- apply(etmp,c(1,3,4),mean)
      etmp2 <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),6,dim(ech_ind$ensmean)[2]))
      ech_ind$ensmean <- apply(etmp2,c(1,3),mean)
      etmp3 <- array(ech_ind$bem,c((dim(ech_ind$bem)[1]/6),6,dim(ech_ind$bem)[2]))
      ech_ind$bem <- apply(etmp3,c(1,3),mean)
      #    ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]/length(unique(ech_ind$names))]
      ech_ind$names <- rep(unique(ech_ind$names),each=dim(ech_ind$data)[1]/
                             length(unique(ech_ind$names)))  
      ech_ind$lon <- ech_ind$lon[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]
      ech_ind$lat <- ech_ind$lat[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))] 
      etmp <- NULL
      etmp2 <- NULL
      
      atmp <- array(analysis.abs$data,c((dim(analysis.abs$data)[1]/6),6,dim(analysis.abs$data)[2],
                                        dim(analysis.abs$data)[3]))  
      analysis.abs$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
      for (ensmem in 1:(dim(atmp)[4])) {
        print(paste('Analysis member',ensmem))
        analysis.abs$data[,,ensmem] <- apply(atmp[,,,ensmem],c(1,3),mean)
      }
      #      analysis.abs$data <- apply(atmp,c(1,3,4),mean)
      #    analysis.abs$names <- analysis.abs$names[1:dim(analysis.abs$data)[1]]
      analysis.abs$names <- rep(unique(analysis.abs$names),each=dim(analysis.abs$data)[1]/
                                  length(unique(analysis.abs$names)))
      analysis.abs$lon <- analysis.abs$lon[1:(dim(analysis.abs$data)[1])] #/length(unique(analysis.abs$names)))]
      analysis.abs$lat <- analysis.abs$lat[1:(dim(analysis.abs$data)[1])] #/length(unique(analysis.abs$names)))]
      atmp2 <- array(analysis.abs$ensmean,c((dim(analysis.abs$ensmean)[1]/6),6,
                                            dim(analysis.abs$ensmean)[2]))
      analysis.abs$ensmean <- apply(atmp2,c(1,3),mean)
      atmp3 <- array(analysis.abs$bem,c((dim(analysis.abs$bem)[1]/6),6,
                                            dim(analysis.abs$bem)[2]))
      analysis.abs$bem <- apply(atmp3,c(1,3),mean)
      atmp <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),6,dim(ana_ind$data)[2],
                                   dim(ana_ind$data)[3]))  
      ana_ind$data <- apply(atmp,c(1,3,4),mean)
      #    ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
      ana_ind$names <- rep(unique(ana_ind$names),each=dim(ana_ind$data)[1]/
                             length(unique(ana_ind$names)))
      ana_ind$lon <- ana_ind$lon[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
      ana_ind$lat <- ana_ind$lat[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
      atmp2 <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),6,
                                       dim(ana_ind$ensmean)[2]))
      ana_ind$ensmean <- apply(atmp2,c(1,3),mean)
      atmp3 <- array(ana_ind$bem,c((dim(ana_ind$bem)[1]/6),6,
                                       dim(ana_ind$bem)[2]))
      ana_ind$bem <- apply(atmp3,c(1,3),mean)
      atmp <- NULL
      atmp2 <- NULL
      
      # same for anomalies
      etmp <- array(echam.anom$data,c((dim(echam.anom$data)[1]/6),6,dim(echam.anom$data)[2],
                                      dim(echam.anom$data)[3]))
      echam.anom$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
      for (ensmem in 1:(dim(etmp)[4])) {
        print(paste('ECHAM anomaly member',ensmem))
        echam.anom$data[,,ensmem] <- apply(etmp[,,,ensmem],c(1,3),mean)
      }
      #    echam.anom$data <- apply(etmp,c(1,3,4),mean)
      etmp2 <- array(echam.anom$ensmean,c((dim(echam.anom$ensmean)[1]/6),6,dim(echam.anom$ensmean)[2]))
      echam.anom$ensmean <- apply(etmp2,c(1,3),mean)
      etmp3 <- array(echam.anom$bem,c((dim(echam.anom$bem)[1]/6),6,dim(echam.anom$bem)[2]))
      echam.anom$bem <- apply(etmp3,c(1,3),mean)
      #    echam.anom$names <- echam.anom$names[1:dim(echam.anom$data)[1]/length(unique(echam.anom$names))]
      echam.anom$names <- rep(unique(echam.anom$names),each=dim(echam.anom$data)[1]/
                                length(unique(echam.anom$names)))  
      echam.anom$lon <- echam.anom$lon[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
      echam.anom$lat <- echam.anom$lat[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
      etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
                                   dim(ech_ind$data)[3]))  
      etmp <- NULL
      
      atmp <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),6,dim(analysis.anom$data)[2],
                                         dim(analysis.anom$data)[3])) 
      analysis.anom$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
      for (ensmem in 1:(dim(atmp)[4])) {
        print(paste('Analysis anomaly member',ensmem))
        analysis.anom$data[,,ensmem] <- apply(atmp[,,,ensmem],c(1,3),mean)
      }
      #      analysis.anom$data <- apply(atmp,c(1,3,4),mean)
      #    analysis.anom$names <- analysis.anom$names[1:dim(analysis.anom$data)[1]]
      analysis.anom$names <- rep(unique(analysis.anom$names),each=dim(analysis.anom$data)[1]/
                                   length(unique(analysis.anom$names)))
      analysis.anom$lon <- analysis.anom$lon[1:(dim(analysis.anom$data)[1])] #/length(unique(analysis.anom$names)))]
      analysis.anom$lat <- analysis.anom$lat[1:(dim(analysis.anom$data)[1])] #/length(unique(analysis.anom$names)))]
      atmp2 <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),6,
                                             dim(analysis.anom$ensmean)[2]))
      analysis.anom$ensmean <- apply(atmp2,c(1,3),mean)
      atmp3 <- array(analysis.anom$bem,c((dim(analysis.anom$bem)[1]/6),6,
                                             dim(analysis.anom$bem)[2]))
      analysis.anom$bem <- apply(atmp3,c(1,3),mean)
      atmp <- NULL
      
      if (vali) {    
        if (!recon_vali) {
          vtmp <- array(validate$data,c((dim(validate$data)[1]/6),6,dim(validate$data)[2]))  
          validate$data <- apply(vtmp,c(1,3),mean)
          validate$ensmean <- validate$data
          #      validate$names <- validate$names[1:dim(validate$data)[1]]
          validate$names <- rep(unique(validate$names),each=dim(validate$data)[1]/
                                  length(unique(validate$names)))
          validate$lon <- validate$lon[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
          validate$lat <- validate$lat[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
          vtmp <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),6,dim(vali_ind$data)[2]))  
          #      vali_ind=validate
          vali_ind$data <- apply(vtmp,c(1,3),mean)
          vali_ind$ensmean <- vali_ind$data
          #      vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
          vali_ind$names <- rep(unique(vali_ind$names),each=dim(vali_ind$data)[1]/
                                  length(unique(vali_ind$names)))
          vali_ind$lon <- vali_ind$lon[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
          vali_ind$lat <- vali_ind$lat[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
        }
      }
      if (!real_proxies) {
        # ERROR, NOT ADJUSTED FOR seasonal/annual DOCUMENTARY DATA, YET !!!
        ctmp <- array(calibrate$data,c((dim(calibrate$data)[1]/6),6,
                                       dim(calibrate$data)[2]))  
        calibrate$data <- apply(ctmp,c(1,3),mean)
        calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]] 
        calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
        calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
      }
    }
    
  





  
# calculate the indices from output 
    H.giorgi <- compute_giorgi_H_v2(giorgi, echam.abs) #, numvar=3) # 3 vars temp, precip, slp
#    H.giorgi.anom <- compute_giorgi_H_sixmon(giorgi, echam.anom) #, numvar=3) # 3 vars temp, precip, slp
#    indices <- c('NH.temp2', 'NH.precip', 'NH.slp', 'NEU.temp2', 'NEU.precip', 'NEU.slp',
#                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
#    Hind <- matrix(0,nrow=12,ncol=nrow(echam$data))
    indices_tmp <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                     'PV1.calc', 'PV2.calc', 'PWC1.calc', 'PWC2.calc', 
                     'DIMI1.calc', 'DIMI2.calc', 'NAO1.calc', 'NAO2.calc',  
                     'PNA1.calc', 'PNA2.calc', 'PNA3.calc', 'PNA4.calc')
    indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                 'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                 'NEU.temp2','MED.temp2',
                 'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                 'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                 'NEU.precip','MED.precip',
                 'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                 'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                 'NEU.slp','MED.slp',
                 'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                 'HC.calc', 'ITCZ.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc', 
                 'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc',
                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
# "MC" kein z300 in output! 'HCL.calc',
# SO WIRD INDEX NUR MITTELWERT ¨UBER REGION UND ZEIT! 
# Deshalb Indicies die min max etc basiert sind extra rechnen
    Hind <- matrix(0,nrow=length(indices_tmp),ncol=nrow(echam.abs$data))
    Hind[1,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                    which(echam.abs$names=="temp2")]
    Hind[2,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                    which(echam.abs$names=="temp2")]
    Hind[3,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                    which(echam.abs$names=="temp2")]
    Hind[4,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                        which(echam.abs$names=="temp2")]
    Hind[5,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                        which(echam.abs$names=="temp2")]
    Hind[6,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                        which(echam.abs$names=="temp2")]
    Hind[7,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                        which(echam.abs$names=="temp2")]
    Hind[8,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                        which(echam.abs$names=="temp2")]
    Hind[9,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                        which(echam.abs$names=="temp2")]
    Hind[10,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                         which(echam.abs$names=="temp2")]
    Hind[11,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                          which(echam.abs$names=="temp2")]
    Hind[12,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                          which(echam.abs$names=="precip")]
    Hind[13,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                          which(echam.abs$names=="precip")]
    Hind[14,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                          which(echam.abs$names=="precip")]
    Hind[15,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                          which(echam.abs$names=="precip")]
    Hind[16,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                          which(echam.abs$names=="precip")]
    Hind[17,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                          which(echam.abs$names=="precip")]
    Hind[18,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                          which(echam.abs$names=="precip")]
    Hind[19,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                          which(echam.abs$names=="precip")]
    Hind[20,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                          which(echam.abs$names=="precip")]
    Hind[21,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                       which(echam.abs$names=="temp2")]
    Hind[22,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                       which(echam.abs$names=="slp")]
    Hind[23,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                       which(echam.abs$names=="slp")]
    Hind[24,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                       which(echam.abs$names=="slp")]
    Hind[25,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                       which(echam.abs$names=="slp")]
    Hind[26,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                       which(echam.abs$names=="slp")]
    Hind[27,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                       which(echam.abs$names=="slp")]
    Hind[28,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                       which(echam.abs$names=="slp")]
    Hind[29,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                       which(echam.abs$names=="slp")]
    Hind[30,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                       which(echam.abs$names=="slp")]  
    Hind[31,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                    which(echam.abs$names=="temp2")]
    Hind[32,which(echam.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                    which(echam.abs$names=="temp2")]
    Hind[33,which(echam.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                     which(echam.abs$names=="precip")]
    Hind[34,which(echam.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                     which(echam.abs$names=="slp")] 
# calc indices from model fields
# midlatitude circulation
    #    Hind[5,which(echam.abs$names=="gph300")] <- H.giorgi[which(giorgi.short == 'MC'),
    #                                              which(echam.abs$names=="gph300")] 
    # cdo -s fldmean -sellonlatbox,0,360,30,60 -sellevel,3000 -selvar,geopoth $f ${STOR_DIR}/z300.$f

# stratospheric polar vortex
    Hind[35,which(echam.abs$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV1'),
                                               which(echam.abs$names=="gph100")] 
    Hind[36,which(echam.abs$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV2'),
                                                     which(echam.abs$names=="gph100")] 
    # cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
    # cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f

# Pacific walker circulation
    Hind[37,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
                                                    which(echam.abs$names=="omega500")] 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f

# Pacific walker circulation
    Hind[38,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC2'),
                                                    which(echam.abs$names=="omega500")] 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f

# Dynamic Indian Monsoon index
    Hind[39,which(echam.abs$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI1'),
                                                   which(echam.abs$names=="u850")] 
    # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f

# Dynamic Indian Monsoon index
    Hind[40,which(echam.abs$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI2'),
                                                which(echam.abs$names=="u850")] 
    # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f

# NAO index
    Hind[41,which(echam.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO1'),
                                               which(echam.anom$names=="slp")] 
    Hind[42,which(echam.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO2'),
                                               which(echam.anom$names=="slp")] 
    Hind[43,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA1'),
                                                  which(echam.anom$names=="gph500")] 
    Hind[44,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA2'),
                                                  which(echam.anom$names=="gph500")] 
    Hind[45,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA3'),
                                                  which(echam.anom$names=="gph500")] 
    Hind[46,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA4'),
                                                  which(echam.anom$names=="gph500")] 
    eind_tmp <- list(ensmean=Hind%*%echam.abs$ensmean, names=indices_tmp)

# Hadley cell strength
    #    Hind[5,which(echam.abs$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
    #                                              which(echam.abs$names=="stream")] 
    # Hind[5,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
    #                                                   which(echam.abs$names=="omega500")] 
    c=1
    zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
    for (l in unique(echam.abs$lat)) {
      pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(echam.abs$names=='omega500')
      pos <- intersect(pos1,pos2)
      if (length(echam.abs$ensmean[pos,])==2) { zmean[c,] <- echam.abs$ensmean[pos,] }
      if (length(echam.abs$ensmean[pos,])>2) { zmean[c,] <- apply(echam.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% which(unique(echam.abs$lat) < 30))]
    pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
    # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f

# Hadley cell poleward extend
    # Hind[6,which(echam.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
    #                                               which(echam.abs$names=="u200")] 

# ITCZ location
    #Hind[7,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
    #                                                   which(echam.abs$names=="omega500")] 
    c=1
    zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
    for (l in unique(echam.abs$lat)) {
      pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(echam.abs$names=='omega500')
      pos <- intersect(pos1,pos2)
      if (length(echam.abs$ensmean[pos,])==2) { zmean[c,] <- echam.abs$ensmean[pos,] }
      if (length(echam.abs$ensmean[pos,])>2) { zmean[c,] <- apply(echam.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% which(unique(echam.abs$lat) < 30))]
    pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(zmean[gridminpos1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(zmean[gridminpos2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    itcz <- c(latmin1,latmin2)

# subtropical jet from subtrop. jet (u200) max
    #Hind[8,which(echam.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
    #                                               which(echam.abs$names=="u200")] 
    c=1
    # zonal means:
    zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
    for (l in unique(echam.abs$lat)) {
      pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(echam.abs$names=='u200')
      pos <- intersect(pos1,pos2)
      if (length(echam.abs$ensmean[pos,])==2) { zmean[c,] <- echam.abs$ensmean[pos,] }
      if (length(echam.abs$ensmean[pos,])>2) { zmean[c,] <- apply(echam.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% which(unique(echam.abs$lat) < 57))]
    pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    # pos1 for season 1; would need adjustment for monthly data
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    # pos2 for season 2
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
    sj <- sj * -1
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    sj_u200 <- c(latmin1,latmin2)
    # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f

    # subtropical jet from zonal SLP max
    c=1
    zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
    for (l in unique(echam.abs$lat)) {
      pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(echam.abs$names=='slp')
      pos <- intersect(pos1,pos2)
      if (length(echam.abs$ensmean[pos,])==2) { zmean[c,] <- echam.abs$ensmean[pos,] }
      if (length(echam.abs$ensmean[pos,])>2) { zmean[c,] <- apply(echam.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% which(unique(echam.abs$lat) < 57))]
    pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
    sj <- sj * -1
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                    (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                    (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    sj_slp <- c(latmin1,latmin2)
    # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
    
    
# stratospheric polar vortex
    pv_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PV1.calc"),]
    pv_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PV2.calc"),]
    pv <- pv_reg1 - pv_reg2 
    #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
    #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f

# Pacific walker circulation
    #Hind[10,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
    #                                                    which(echam.abs$names=="omega500")] 
    pwc_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PWC1.calc"),]
    pwc_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PWC2.calc"),]
    pwc <- pwc_reg1 - pwc_reg2 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f

# Dynamic Indian Monsoon index
    dimi_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="DIMI1.calc"),]
    dimi_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="DIMI2.calc"),]
    dimi <- dimi_reg1 - dimi_reg2 
    #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f


# NAO based on anomalies
    nao_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="NAO1.calc"),]
    nao_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="NAO2.calc"),]
    nao <- nao_reg1 - nao_reg2 

# PNA based on anomalies
    pna_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA1.calc"),]
    pna_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA2.calc"),]
    pna_reg3 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA3.calc"),]
    pna_reg4 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA4.calc"),]
    pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
    # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 

    eind <- list(ensmean=matrix(0,nrow=length(indices),ncol=2), names=indices)  
    eind$ensmean[which(indices=='ENH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ENH.temp2'),]
    eind$ensmean[which(indices=='NAM.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NAM.temp2'),]
    eind$ensmean[which(indices=='SAM.temp2'),] <- eind_tmp$ensmean[which(eind$names=='SAM.temp2'),]
    eind$ensmean[which(indices=='AFR.temp2'),] <- eind_tmp$ensmean[which(eind$names=='AFR.temp2'),]
    eind$ensmean[which(indices=='ASI.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ASI.temp2'),]
    eind$ensmean[which(indices=='AUS.temp2'),] <- eind_tmp$ensmean[which(eind$names=='AUS.temp2'),]
    eind$ensmean[which(indices=='ARC.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ARC.temp2'),]
    eind$ensmean[which(indices=='ANT.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ANT.temp2'),]
    eind$ensmean[which(indices=='NEU.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NEU.temp2'),]
    eind$ensmean[which(indices=='MED.temp2'),] <- eind_tmp$ensmean[which(eind$names=='MED.temp2'),]
    eind$ensmean[which(indices=='GLO.temp2'),] <- eind_tmp$ensmean[which(eind$names=='GLO.temp2'),]
    eind$ensmean[which(indices=='NAM.precip'),] <- eind_tmp$ensmean[which(eind$names=='NAM.precip'),]
    eind$ensmean[which(indices=='SAM.precip'),] <- eind_tmp$ensmean[which(eind$names=='SAM.precip'),]
    eind$ensmean[which(indices=='AFR.precip'),] <- eind_tmp$ensmean[which(eind$names=='AFR.precip'),]
    eind$ensmean[which(indices=='ASI.precip'),] <- eind_tmp$ensmean[which(eind$names=='ASI.precip'),]
    eind$ensmean[which(indices=='AUS.precip'),] <- eind_tmp$ensmean[which(eind$names=='AUS.precip'),]
    eind$ensmean[which(indices=='ARC.precip'),] <- eind_tmp$ensmean[which(eind$names=='ARC.precip'),]
    eind$ensmean[which(indices=='ANT.precip'),] <- eind_tmp$ensmean[which(eind$names=='ANT.precip'),]
    eind$ensmean[which(indices=='NEU.precip'),] <- eind_tmp$ensmean[which(eind$names=='NEU.precip'),]
    eind$ensmean[which(indices=='MED.precip'),] <- eind_tmp$ensmean[which(eind$names=='MED.precip'),]
    eind$ensmean[which(indices=='SH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='SH.temp2'),]
    eind$ensmean[which(indices=='NAM.slp'),] <- eind_tmp$ensmean[which(eind$names=='NAM.slp'),]
    eind$ensmean[which(indices=='SAM.slp'),] <- eind_tmp$ensmean[which(eind$names=='SAM.slp'),]
    eind$ensmean[which(indices=='AFR.slp'),] <- eind_tmp$ensmean[which(eind$names=='AFR.slp'),]
    eind$ensmean[which(indices=='ASI.slp'),] <- eind_tmp$ensmean[which(eind$names=='ASI.slp'),]
    eind$ensmean[which(indices=='AUS.slp'),] <- eind_tmp$ensmean[which(eind$names=='AUS.slp'),]
    eind$ensmean[which(indices=='ARC.slp'),] <- eind_tmp$ensmean[which(eind$names=='ARC.slp'),]
    eind$ensmean[which(indices=='ANT.slp'),] <- eind_tmp$ensmean[which(eind$names=='ANT.slp'),]
    eind$ensmean[which(indices=='NEU.slp'),] <- eind_tmp$ensmean[which(eind$names=='NEU.slp'),]
    eind$ensmean[which(indices=='MED.slp'),] <- eind_tmp$ensmean[which(eind$names=='MED.slp'),]
    eind$ensmean[which(indices=='NH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NH.temp2'),]
    eind$ensmean[which(indices=='EU.temp2'),] <- eind_tmp$ensmean[which(eind$names=='EU.temp2'),]
    eind$ensmean[which(indices=='EU.precip'),] <- eind_tmp$ensmean[which(eind$names=='EU.precip'),]
    eind$ensmean[which(indices=='EU.slp'),] <- eind_tmp$ensmean[which(eind$names=='EU.slp'),]
    eind$ensmean[which(indices=='HC.calc'),] <- hc #eind_tmp$ensmean[which(eind$names=='HC.calc'),]
    eind$ensmean[which(indices=='ITCZ.calc'),] <- itcz #eind_tmp$ensmean[which(eind$names=='ITCZ.calc'),]
    eind$ensmean[which(indices=='SJ_u200.calc'),] <- sj_u200 #eind_tmp$ensmean[which(eind$names=='SJ.calc'),]
    eind$ensmean[which(indices=='SJ_slp.calc'),] <- sj_slp
    eind$ensmean[which(indices=='PV.calc'),] <- pv #eind_tmp$ensmean[which(eind$names=='PV.calc'),]
    eind$ensmean[which(indices=='PWC.calc'),] <- pwc #eind_tmp$ensmean[which(eind$names=='PWC.calc'),]
    eind$ensmean[which(indices=='DIMI.calc'),] <- dimi #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),]
    eind$ensmean[which(indices=='NAO.calc'),] <- nao #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),] #c(0,0)
    eind$ensmean[which(indices=='PNA.calc'),] <- pna #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),] #c(0,0)
    eind$ensmean[which(indices=='DIMI'),] <- ech_ind$ensmean[which(ech_ind$names=='DIMI'),]
    eind$ensmean[which(indices=='HC'),] <- ech_ind$ensmean[which(ech_ind$names=='HC'),]
    eind$ensmean[which(indices=='SJ'),] <- ech_ind$ensmean[which(ech_ind$names=='SJ'),]
    eind$ensmean[which(indices=='Z100'),] <- ech_ind$ensmean[which(ech_ind$names=='z100'),]
    eind$ensmean[which(indices=='Z300'),] <- ech_ind$ensmean[which(ech_ind$names=='z300'),]
    eind$ensmean[which(indices=='PWC'),] <- ech_ind$ensmean[which(ech_ind$names=='PWC'),]



# BEM
eind_tmp <- list(bem=Hind%*%echam.abs$bem, names=indices_tmp)

# Hadley cell strength
#    Hind[5,which(echam.abs$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
#                                              which(echam.abs$names=="stream")] 
# Hind[5,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
#                                                   which(echam.abs$names=="omega500")] 
c=1
zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
for (l in unique(echam.abs$lat)) {
  pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(echam.abs$names=='omega500')
  pos <- intersect(pos1,pos2)
  if (length(echam.abs$bem[pos,])==2) { zmean[c,] <- echam.abs$bem[pos,] }
  if (length(echam.abs$bem[pos,])>2) { zmean[c,] <- apply(echam.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% which(unique(echam.abs$lat) < 30))]
pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
# cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f

# Hadley cell poleward extend
# Hind[6,which(echam.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
#                                               which(echam.abs$names=="u200")] 

# ITCZ location
#Hind[7,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
#                                                   which(echam.abs$names=="omega500")] 
c=1
zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
for (l in unique(echam.abs$lat)) {
  pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(echam.abs$names=='omega500')
  pos <- intersect(pos1,pos2)
  if (length(echam.abs$bem[pos,])==2) { zmean[c,] <- echam.abs$bem[pos,] }
  if (length(echam.abs$bem[pos,])>2) { zmean[c,] <- apply(echam.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% which(unique(echam.abs$lat) < 30))]
pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(itcz[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(itcz[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
itcz <- c(latmin1,latmin2)

# subtropical jet (u200)
#Hind[8,which(echam.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
#                                               which(echam.abs$names=="u200")] 
c=1
zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
for (l in unique(echam.abs$lat)) {
  pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(echam.abs$names=='u200')
  pos <- intersect(pos1,pos2)
  if (length(echam.abs$bem[pos,])==2) { zmean[c,] <- echam.abs$bem[pos,] }
  if (length(echam.abs$bem[pos,])>2) { zmean[c,] <- apply(echam.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% which(unique(echam.abs$lat) < 57))]
pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
sj <- sj * -1
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
sj_u200 <- c(latmin1,latmin2)
# cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f

# subtropical jet (slp)
c=1
zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
for (l in unique(echam.abs$lat)) {
  pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(echam.abs$names=='slp')
  pos <- intersect(pos1,pos2)
  if (length(echam.abs$bem[pos,])==2) { zmean[c,] <- echam.abs$bem[pos,] }
  if (length(echam.abs$bem[pos,])>2) { zmean[c,] <- apply(echam.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% which(unique(echam.abs$lat) < 57))]
pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
sj <- sj * -1
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
sj_slp <- c(latmin1,latmin2)

# stratospheric polar vortex
pv_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PV1.calc"),]
pv_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PV2.calc"),]
pv <- pv_reg1 - pv_reg2 
#cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
#cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f

# Pacific walker circulation
#Hind[10,which(echam.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
#                                                    which(echam.abs$names=="omega500")] 
pwc_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PWC1.calc"),]
pwc_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PWC2.calc"),]
pwc <- pwc_reg1 - pwc_reg2 
# cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
# cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f

# Dynamic Indian Monsoon index
dimi_reg1 <- eind_tmp$bem[which(eind_tmp$names=="DIMI1.calc"),]
dimi_reg2 <- eind_tmp$bem[which(eind_tmp$names=="DIMI2.calc"),]
dimi <- dimi_reg1 - dimi_reg2 
#cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
#cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f


# NAO based on anomalies
nao_reg1 <- eind_tmp$bem[which(eind_tmp$names=="NAO1.calc"),]
nao_reg2 <- eind_tmp$bem[which(eind_tmp$names=="NAO2.calc"),]
nao <- nao_reg1 - nao_reg2 

# PNA based on anomalies
pna_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PNA1.calc"),]
pna_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PNA2.calc"),]
pna_reg3 <- eind_tmp$bem[which(eind_tmp$names=="PNA3.calc"),]
pna_reg4 <- eind_tmp$bem[which(eind_tmp$names=="PNA4.calc"),]
pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
# PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 

eind$bem <- matrix(0,nrow=length(indices),ncol=2)
#eind <- list(bem=matrix(0,nrow=length(indices),ncol=2), names=indices)  
eind$bem[which(indices=='ENH.temp2'),] <- eind_tmp$bem[which(eind$names=='ENH.temp2'),]
eind$bem[which(indices=='NAM.temp2'),] <- eind_tmp$bem[which(eind$names=='NAM.temp2'),]
eind$bem[which(indices=='SAM.temp2'),] <- eind_tmp$bem[which(eind$names=='SAM.temp2'),]
eind$bem[which(indices=='AFR.temp2'),] <- eind_tmp$bem[which(eind$names=='AFR.temp2'),]
eind$bem[which(indices=='ASI.temp2'),] <- eind_tmp$bem[which(eind$names=='ASI.temp2'),]
eind$bem[which(indices=='AUS.temp2'),] <- eind_tmp$bem[which(eind$names=='AUS.temp2'),]
eind$bem[which(indices=='ARC.temp2'),] <- eind_tmp$bem[which(eind$names=='ARC.temp2'),]
eind$bem[which(indices=='ANT.temp2'),] <- eind_tmp$bem[which(eind$names=='ANT.temp2'),]
eind$bem[which(indices=='NEU.temp2'),] <- eind_tmp$bem[which(eind$names=='NEU.temp2'),]
eind$bem[which(indices=='MED.temp2'),] <- eind_tmp$bem[which(eind$names=='MED.temp2'),]
eind$bem[which(indices=='GLO.temp2'),] <- eind_tmp$bem[which(eind$names=='GLO.temp2'),]
eind$bem[which(indices=='NAM.precip'),] <- eind_tmp$bem[which(eind$names=='NAM.precip'),]
eind$bem[which(indices=='SAM.precip'),] <- eind_tmp$bem[which(eind$names=='SAM.precip'),]
eind$bem[which(indices=='AFR.precip'),] <- eind_tmp$bem[which(eind$names=='AFR.precip'),]
eind$bem[which(indices=='ASI.precip'),] <- eind_tmp$bem[which(eind$names=='ASI.precip'),]
eind$bem[which(indices=='AUS.precip'),] <- eind_tmp$bem[which(eind$names=='AUS.precip'),]
eind$bem[which(indices=='ARC.precip'),] <- eind_tmp$bem[which(eind$names=='ARC.precip'),]
eind$bem[which(indices=='ANT.precip'),] <- eind_tmp$bem[which(eind$names=='ANT.precip'),]
eind$bem[which(indices=='NEU.precip'),] <- eind_tmp$bem[which(eind$names=='NEU.precip'),]
eind$bem[which(indices=='MED.precip'),] <- eind_tmp$bem[which(eind$names=='MED.precip'),]
eind$bem[which(indices=='SH.temp2'),] <- eind_tmp$bem[which(eind$names=='SH.temp2'),]
eind$bem[which(indices=='NAM.slp'),] <- eind_tmp$bem[which(eind$names=='NAM.slp'),]
eind$bem[which(indices=='SAM.slp'),] <- eind_tmp$bem[which(eind$names=='SAM.slp'),]
eind$bem[which(indices=='AFR.slp'),] <- eind_tmp$bem[which(eind$names=='AFR.slp'),]
eind$bem[which(indices=='ASI.slp'),] <- eind_tmp$bem[which(eind$names=='ASI.slp'),]
eind$bem[which(indices=='AUS.slp'),] <- eind_tmp$bem[which(eind$names=='AUS.slp'),]
eind$bem[which(indices=='ARC.slp'),] <- eind_tmp$bem[which(eind$names=='ARC.slp'),]
eind$bem[which(indices=='ANT.slp'),] <- eind_tmp$bem[which(eind$names=='ANT.slp'),]
eind$bem[which(indices=='NEU.slp'),] <- eind_tmp$bem[which(eind$names=='NEU.slp'),]
eind$bem[which(indices=='MED.slp'),] <- eind_tmp$bem[which(eind$names=='MED.slp'),]
eind$bem[which(indices=='NH.temp2'),] <- eind_tmp$bem[which(eind$names=='NH.temp2'),]
eind$bem[which(indices=='EU.temp2'),] <- eind_tmp$bem[which(eind$names=='EU.temp2'),]
eind$bem[which(indices=='EU.precip'),] <- eind_tmp$bem[which(eind$names=='EU.precip'),]
eind$bem[which(indices=='EU.slp'),] <- eind_tmp$bem[which(eind$names=='EU.slp'),]
eind$bem[which(indices=='HC.calc'),] <- hc #eind_tmp$bem[which(eind$names=='HC.calc'),]
eind$bem[which(indices=='ITCZ.calc'),] <- itcz #eind_tmp$bem[which(eind$names=='ITCZ.calc'),]
eind$bem[which(indices=='SJ_u200.calc'),] <- sj_u200 #eind_tmp$bem[which(eind$names=='SJ.calc'),]
eind$bem[which(indices=='SJ_slp.calc'),] <- sj_slp
eind$bem[which(indices=='PV.calc'),] <- pv #eind_tmp$bem[which(eind$names=='PV.calc'),]
eind$bem[which(indices=='PWC.calc'),] <- pwc #eind_tmp$bem[which(eind$names=='PWC.calc'),]
eind$bem[which(indices=='DIMI.calc'),] <- dimi #eind_tmp$bem[which(eind$names=='DIMI.calc'),]
eind$bem[which(indices=='NAO.calc'),] <- nao #eind_tmp$bem[which(eind$names=='DIMI.calc'),] #c(0,0)
eind$bem[which(indices=='PNA.calc'),] <- pna #eind_tmp$bem[which(eind$names=='DIMI.calc'),] #c(0,0)
eind$bem[which(indices=='DIMI'),] <- ech_ind$bem[which(ech_ind$names=='DIMI'),]
eind$bem[which(indices=='HC'),] <- ech_ind$bem[which(ech_ind$names=='HC'),]
eind$bem[which(indices=='SJ'),] <- ech_ind$bem[which(ech_ind$names=='SJ'),]
eind$bem[which(indices=='Z100'),] <- ech_ind$bem[which(ech_ind$names=='z100'),]
eind$bem[which(indices=='Z300'),] <- ech_ind$bem[which(ech_ind$names=='z300'),]
eind$bem[which(indices=='PWC'),] <- ech_ind$bem[which(ech_ind$names=='PWC'),]







# data for ensemble member
#    eind$data <- array(Hind %*% array(echam$data, c(nrow(echam$data), 
#                       length(echam$data)/nrow(echam$data))), c(nrow(Hind), 
#                       dim(echam$data)[2:3]))
    eind_tmp <- list(data=array(Hind %*% array(echam.abs$data, c(nrow(echam.abs$data), 
                       length(echam.abs$data)/nrow(echam.abs$data))), c(nrow(Hind), 
                       dim(echam.abs$data)[2:3])), names=indices_tmp)
# Hadley cell strength
    hc <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
      for (l in unique(echam.abs$lat)) {
        pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(echam.abs$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam.abs$data[pos,,m])==2) { zmean[c,] <- echam.abs$data[pos,,m] }
        if (length(echam.abs$data[pos,,m])>2) { zmean[c,] <- apply(echam.abs$data[pos,,m],2,mean) }
#        zmean[c,] <- apply(echam.abs$data[pos,,m],2,mean)
        c=c+1
      }
      latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% 
                   which(unique(echam.abs$lat) < 30))]
      pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      hc[m,] <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
    }
# Hadley cell poleward extend

# ITCZ location
    itcz <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
      for (l in unique(echam.abs$lat)) {
        pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(echam.abs$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam.abs$data[pos,,m])==2) { zmean[c,] <- echam.abs$data[pos,,m] }
        if (length(echam.abs$data[pos,,m])>2) { zmean[c,] <- apply(echam.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > -30) %in% 
                    which(unique(echam.abs$lat) < 30))]
      pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      itcz_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                   (2*(abs(itcz_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                   (2*(abs(itcz_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz[m,] <- c(latmin1,latmin2)
    }

# subtropical jet (u200)
    sj_u200 <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
      for (l in unique(echam.abs$lat)) {
        pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(echam.abs$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(echam.abs$data[pos,,m])==2) { zmean[c,] <- echam.abs$data[pos,,m] }
        if (length(echam.abs$data[pos,,m])>2) { zmean[c,] <- apply(echam.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% 
                    which(unique(echam.abs$lat) < 57))]
      pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj_tmp <- sj_tmp * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                  (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                  (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200[m,] <- c(latmin1,latmin2)
    }
    
    # subtropical jet (slp)
    sj_slp <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(echam.abs$lat)),ncol=2)
      for (l in unique(echam.abs$lat)) {
        pos1 <- which(trunc(echam.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(echam.abs$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(echam.abs$data[pos,,m])==2) { zmean[c,] <- echam.abs$data[pos,,m] }
        if (length(echam.abs$data[pos,,m])>2) { zmean[c,] <- apply(echam.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(echam.abs$lat)[which(which(unique(echam.abs$lat) > 0) %in% 
                                               which(unique(echam.abs$lat) < 57))]
      pos3 <- match(round(latlist,digits=5),round(unique(echam.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj_tmp <- sj_tmp * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp[m,] <- c(latmin1,latmin2)
    }    

# stratospheric polar vortex
    pv <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pv_reg1 <- eind_tmp$data[which(eind_tmp$names=="PV1.calc"),,m]
      pv_reg2 <- eind_tmp$data[which(eind_tmp$names=="PV2.calc"),,m]
      pv[m,] <- pv_reg1 - pv_reg2 
    }

# Pacific walker circulation
    pwc <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pwc_reg1 <- eind_tmp$data[which(eind_tmp$names=="PWC1.calc"),,m]
      pwc_reg2 <- eind_tmp$data[which(eind_tmp$names=="PWC2.calc"),,m]
      pwc[m,] <- pwc_reg1 - pwc_reg2 
    }

# Dynamic Indian Monsoon index
    dimi <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      dimi_reg1 <- eind_tmp$data[which(eind_tmp$names=="DIMI1.calc"),,m]
      dimi_reg2 <- eind_tmp$data[which(eind_tmp$names=="DIMI2.calc"),,m]
      dimi[m,] <- dimi_reg1 - dimi_reg2 
    }

# NAO based on anomalies
    nao <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      nao_reg1 <- eind_tmp$data[which(eind_tmp$names=="NAO1.calc"),,m]
      nao_reg2 <- eind_tmp$data[which(eind_tmp$names=="NAO2.calc"),,m]
      nao[m,] <- nao_reg1 - nao_reg2 
    }
      
# PNA based on anomalies
    pna <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pna_reg1 <- eind_tmp$data[which(eind_tmp$names=="PNA1.calc"),,m]
      pna_reg2 <- eind_tmp$data[which(eind_tmp$names=="PNA2.calc"),,m]
      pna_reg3 <- eind_tmp$data[which(eind_tmp$names=="PNA3.calc"),,m]
      pna_reg4 <- eind_tmp$data[which(eind_tmp$names=="PNA4.calc"),,m]
      pna[m,] <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
    }

    eind$data <- array(0,c(length(indices),2,nmem))
    for (m in 1:nmem) {
      eind$data[which(indices=='ENH.temp2'),,m] <- eind_tmp$data[which(eind$names=='ENH.temp2'),,m]
      eind$data[which(indices=='NAM.temp2'),,m] <- eind_tmp$data[which(eind$names=='NAM.temp2'),,m]
      eind$data[which(indices=='SAM.temp2'),,m] <- eind_tmp$data[which(eind$names=='SAM.temp2'),,m]
      eind$data[which(indices=='AFR.temp2'),,m] <- eind_tmp$data[which(eind$names=='AFR.temp2'),,m]
      eind$data[which(indices=='ASI.temp2'),,m] <- eind_tmp$data[which(eind$names=='ASI.temp2'),,m]
      eind$data[which(indices=='AUS.temp2'),,m] <- eind_tmp$data[which(eind$names=='AUS.temp2'),,m]
      eind$data[which(indices=='ARC.temp2'),,m] <- eind_tmp$data[which(eind$names=='ARC.temp2'),,m]
      eind$data[which(indices=='ANT.temp2'),,m] <- eind_tmp$data[which(eind$names=='ANT.temp2'),,m]
      eind$data[which(indices=='NEU.temp2'),,m] <- eind_tmp$data[which(eind$names=='NEU.temp2'),,m]
      eind$data[which(indices=='MED.temp2'),,m] <- eind_tmp$data[which(eind$names=='MED.temp2'),,m]
      eind$data[which(indices=='GLO.temp2'),,m] <- eind_tmp$data[which(eind$names=='GLO.temp2'),,m]
      eind$data[which(indices=='NAM.precip'),,m] <- eind_tmp$data[which(eind$names=='NAM.precip'),,m]
      eind$data[which(indices=='SAM.precip'),,m] <- eind_tmp$data[which(eind$names=='SAM.precip'),,m]
      eind$data[which(indices=='AFR.precip'),,m] <- eind_tmp$data[which(eind$names=='AFR.precip'),,m]
      eind$data[which(indices=='ASI.precip'),,m] <- eind_tmp$data[which(eind$names=='ASI.precip'),,m]
      eind$data[which(indices=='AUS.precip'),,m] <- eind_tmp$data[which(eind$names=='AUS.precip'),,m]
      eind$data[which(indices=='ARC.precip'),,m] <- eind_tmp$data[which(eind$names=='ARC.precip'),,m]
      eind$data[which(indices=='ANT.precip'),,m] <- eind_tmp$data[which(eind$names=='ANT.precip'),,m]
      eind$data[which(indices=='NEU.precip'),,m] <- eind_tmp$data[which(eind$names=='NEU.precip'),,m]
      eind$data[which(indices=='MED.precip'),,m] <- eind_tmp$data[which(eind$names=='MED.precip'),,m]
      eind$data[which(indices=='SH.temp2'),,m] <- eind_tmp$data[which(eind$names=='SH.temp2'),,m]
      eind$data[which(indices=='NAM.slp'),,m] <- eind_tmp$data[which(eind$names=='NAM.slp'),,m]
      eind$data[which(indices=='SAM.slp'),,m] <- eind_tmp$data[which(eind$names=='SAM.slp'),,m]
      eind$data[which(indices=='AFR.slp'),,m] <- eind_tmp$data[which(eind$names=='AFR.slp'),,m]
      eind$data[which(indices=='ASI.slp'),,m] <- eind_tmp$data[which(eind$names=='ASI.slp'),,m]
      eind$data[which(indices=='AUS.slp'),,m] <- eind_tmp$data[which(eind$names=='AUS.slp'),,m]
      eind$data[which(indices=='ARC.slp'),,m] <- eind_tmp$data[which(eind$names=='ARC.slp'),,m]
      eind$data[which(indices=='ANT.slp'),,m] <- eind_tmp$data[which(eind$names=='ANT.slp'),,m]
      eind$data[which(indices=='NEU.slp'),,m] <- eind_tmp$data[which(eind$names=='NEU.slp'),,m]
      eind$data[which(indices=='MED.slp'),,m] <- eind_tmp$data[which(eind$names=='MED.slp'),,m]
      eind$data[which(indices=='NH.temp2'),,m] <- eind_tmp$data[which(eind$names=='NH.temp2'),,m]
      eind$data[which(indices=='EU.temp2'),,m] <- eind_tmp$data[which(eind$names=='EU.temp2'),,m]
      eind$data[which(indices=='EU.precip'),,m] <- eind_tmp$data[which(eind$names=='EU.precip'),,m]
      eind$data[which(indices=='EU.slp'),,m] <- eind_tmp$data[which(eind$names=='EU.slp'),,m]
      eind$data[which(indices=='HC.calc'),,m] <- hc[m,] #eind_tmp$data[which(eind$names=='HC.calc'),,m]
      eind$data[which(indices=='ITCZ.calc'),,m] <- itcz[m,] #eind_tmp$data[which(eind$names=='ITCZ.calc'),,m]
      eind$data[which(indices=='SJ_u200.calc'),,m] <- sj_u200[m,] #eind_tmp$data[which(eind$names=='SJ.calc'),,m]
      eind$data[which(indices=='SJ_slp.calc'),,m] <- sj_slp[m,]
      eind$data[which(indices=='PV.calc'),,m] <- pv[m,] #eind_tmp$data[which(eind$names=='PV.calc'),,m]
      eind$data[which(indices=='PWC.calc'),,m] <- pwc[m,] #eind_tmp$data[which(eind$names=='PWC.calc'),,m]
      eind$data[which(indices=='DIMI.calc'),,m] <- dimi[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m]
      eind$data[which(indices=='NAO.calc'),,m] <- nao[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m] #c(0,0)
      eind$data[which(indices=='PNA.calc'),,m] <- pna[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m] #c(0,0)
      eind$data[which(indices=='DIMI'),,m] <- ech_ind$data[which(ech_ind$names=='DIMI'),,m]
      eind$data[which(indices=='HC'),,m] <- ech_ind$data[which(ech_ind$names=='HC'),,m]
      eind$data[which(indices=='SJ'),,m] <- ech_ind$data[which(ech_ind$names=='SJ'),,m]
      eind$data[which(indices=='Z100'),,m] <- ech_ind$data[which(ech_ind$names=='z100'),,m]
      eind$data[which(indices=='Z300'),,m] <- ech_ind$data[which(ech_ind$names=='z300'),,m]
      eind$data[which(indices=='PWC'),,m] <- ech_ind$data[which(ech_ind$names=='PWC'),,m]
    }

    eind$data[which(indices=='DIMI'),,] <- ech_ind$data[which(ech_ind$names=='DIMI'),,]
    eind$data[which(indices=='HC'),,] <- ech_ind$data[which(ech_ind$names=='HC'),,]
    eind$data[which(indices=='SJ'),,] <- ech_ind$data[which(ech_ind$names=='SJ'),,]
    eind$data[which(indices=='Z100'),,] <- ech_ind$data[which(ech_ind$names=='z100'),,]
    eind$data[which(indices=='Z300'),,] <- ech_ind$data[which(ech_ind$names=='z300'),,]
    eind$data[which(indices=='PWC'),,] <- ech_ind$data[which(ech_ind$names=='PWC'),,]




    # SAME FOR ANALYSIS
    H.giorgi <- compute_giorgi_H_v2(giorgi, analysis.abs) #, numvar=3) # 3 vars temp, precip, slp
#    H.giorgi <- compute_giorgi_H_sixmon(giorgi, analysis.abs) #, numvar=3) # 3 vars temp, precip, slp
#    H.giorgi.anom <- compute_giorgi_H_sixmon(giorgi, analysis.anom) #, numvar=3) # 3 vars temp, precip, slp
    #    indices <- c('NH.temp2', 'NH.precip', 'NH.slp', 'NEU.temp2', 'NEU.precip', 'NEU.slp',
    #                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
    #    Hind <- matrix(0,nrow=12,ncol=nrow(analysis$data))
    indices_tmp <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                     'PV1.calc', 'PV2.calc', 'PWC1.calc', 'PWC2.calc', 
                     'DIMI1.calc', 'DIMI2.calc', 'NAO1.calc', 'NAO2.calc',  
                     'PNA1.calc', 'PNA2.calc', 'PNA3.calc', 'PNA4.calc')
    indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                 'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                 'NEU.temp2','MED.temp2',
                 'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                 'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                 'NEU.precip','MED.precip',
                 'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                 'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                 'NEU.slp','MED.slp',
                 'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                 'HC.calc', 'ITCZ.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc', 
                 'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc',
                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
    # "MC" kein z300 in output! 'HCL.calc',
    # SO WIRD INDEX NUR MITTELWERT ¨UBER REGION UND ZEIT! 
    # Deshalb Indicies die min max etc basiert sind extra rechnen
    Hind <- matrix(0,nrow=length(indices_tmp),ncol=nrow(analysis.abs$data))
    Hind[1,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[2,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[3,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[4,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[5,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[6,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[7,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[8,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[9,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                           which(analysis.abs$names=="temp2")]
    Hind[10,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                            which(analysis.abs$names=="temp2")]
    Hind[11,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                             which(analysis.abs$names=="temp2")]
    Hind[12,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                             which(analysis.abs$names=="precip")]
    Hind[13,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                             which(analysis.abs$names=="precip")]
    Hind[14,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                             which(analysis.abs$names=="precip")]
    Hind[15,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                             which(analysis.abs$names=="precip")]
    Hind[16,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                             which(analysis.abs$names=="precip")]
    Hind[17,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                             which(analysis.abs$names=="precip")]
    Hind[18,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                             which(analysis.abs$names=="precip")]
    Hind[19,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                             which(analysis.abs$names=="precip")]
    Hind[20,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                             which(analysis.abs$names=="precip")]
    Hind[21,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                          which(analysis.abs$names=="temp2")]
    Hind[22,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                          which(analysis.abs$names=="slp")]
    Hind[23,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                          which(analysis.abs$names=="slp")]
    Hind[24,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                          which(analysis.abs$names=="slp")]
    Hind[25,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                          which(analysis.abs$names=="slp")]
    Hind[26,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                          which(analysis.abs$names=="slp")]
    Hind[27,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                          which(analysis.abs$names=="slp")]
    Hind[28,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                          which(analysis.abs$names=="slp")]
    Hind[29,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                          which(analysis.abs$names=="slp")]
    Hind[30,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                          which(analysis.abs$names=="slp")]  
    
    Hind[31,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                               which(analysis.abs$names=="temp2")]
    Hind[32,which(analysis.abs$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                               which(analysis.abs$names=="temp2")]
    Hind[33,which(analysis.abs$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                                which(analysis.abs$names=="precip")]
    Hind[34,which(analysis.abs$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                             which(analysis.abs$names=="slp")] 
    # calc indices from model fields
    # midlatitude circulation
    #    Hind[5,which(analysis.abs$names=="gph300")] <- H.giorgi[which(giorgi.short == 'MC'),
    #                                              which(analysis.abs$names=="gph300")] 
    # cdo -s fldmean -sellonlatbox,0,360,30,60 -sellevel,3000 -selvar,geopoth $f ${STOR_DIR}/z300.$f
    
    # stratospheric polar vortex
    Hind[35,which(analysis.abs$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV1'),
                                                                which(analysis.abs$names=="gph100")] 
    Hind[36,which(analysis.abs$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV2'),
                                                                which(analysis.abs$names=="gph100")] 
    # cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
    # cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
    
    # Pacific walker circulation
    Hind[37,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
                                                                  which(analysis.abs$names=="omega500")] 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
    
    # Pacific walker circulation
    Hind[38,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC2'),
                                                                  which(analysis.abs$names=="omega500")] 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
    
    # Dynamic Indian Monsoon index
    Hind[39,which(analysis.abs$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI1'),
                                                              which(analysis.abs$names=="u850")] 
    # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
    
    # Dynamic Indian Monsoon index
    Hind[40,which(analysis.abs$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI2'),
                                                               which(analysis.abs$names=="u850")] 
    # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
    
    # NAO index
    Hind[41,which(analysis.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO1'),
                                                               which(analysis.anom$names=="slp")] 
    Hind[42,which(analysis.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO2'),
                                                               which(analysis.anom$names=="slp")] 
    Hind[43,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA1'),
                                                                  which(analysis.anom$names=="gph500")] 
    Hind[44,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA2'),
                                                                  which(analysis.anom$names=="gph500")] 
    Hind[45,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA3'),
                                                                  which(analysis.anom$names=="gph500")] 
    Hind[46,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA4'),
                                                                  which(analysis.anom$names=="gph500")] 
    aind_tmp <- list(ensmean=Hind%*%analysis.abs$ensmean, names=indices_tmp)
    
    # Hadley cell strength
    #    Hind[5,which(analysis.abs$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
    #                                              which(analysis.abs$names=="stream")] 
    # Hind[5,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
    #                                                   which(analysis.abs$names=="omega500")] 
    c=1
    zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
    for (l in unique(analysis.abs$lat)) {
      pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(analysis.abs$names=='omega500')
      pos <- intersect(pos1,pos2)
      if (length(analysis.abs$ensmean[pos,])==2) { zmean[c,] <- analysis.abs$ensmean[pos,] }
      if (length(analysis.abs$ensmean[pos,])>2) { zmean[c,] <- apply(analysis.abs$ensmean[pos,],2,mean) }
#      zmean[c,] <- apply(analysis.abs$ensmean[pos,],2,mean)
      c=c+1
    }
    latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                                which(unique(analysis.abs$lat) < 30))]
    pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
    # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
    
    # Hadley cell poleward extend
    # Hind[6,which(analysis.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
    #                                               which(analysis.abs$names=="u200")] 
    
    # ITCZ location
    #Hind[7,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
    #                                                   which(analysis.abs$names=="omega500")] 
    c=1
    zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
    for (l in unique(analysis.abs$lat)) {
      pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(analysis.abs$names=='omega500')
      pos <- intersect(pos1,pos2)
      if (length(analysis.abs$ensmean[pos,])==2) { zmean[c,] <- analysis.abs$ensmean[pos,] }
      if (length(analysis.abs$ensmean[pos,])>2) { zmean[c,] <- apply(analysis.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                                which(unique(analysis.abs$lat) < 30))]
    pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                    (2*(abs(itcz[1]-min(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                    (2*(abs(itcz[2]-min(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    itcz <- c(latmin1,latmin2)
    
    # subtropical jet (u200)
    #Hind[8,which(analysis.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
    #                                               which(analysis.abs$names=="u200")] 
    c=1
    zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
    for (l in unique(analysis.abs$lat)) {
      pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(analysis.abs$names=='u200')
      pos <- intersect(pos1,pos2)
      if (length(analysis.abs$ensmean[pos,])==2) { zmean[c,] <- analysis.abs$ensmean[pos,] }
      if (length(analysis.abs$ensmean[pos,])>2) { zmean[c,] <- apply(analysis.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% which(unique(analysis.abs$lat) < 57))]
    pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
    sj <- sj * -1
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                    (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                    (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    sj_u200 <- c(latmin1,latmin2)
    # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
 
    # subtropical jet (slp)
    c=1
    zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
    for (l in unique(analysis.abs$lat)) {
      pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
      pos2 <- which(analysis.abs$names=='slp')
      pos <- intersect(pos1,pos2)
      if (length(analysis.abs$ensmean[pos,])==2) { zmean[c,] <- analysis.abs$ensmean[pos,] }
      if (length(analysis.abs$ensmean[pos,])>2) { zmean[c,] <- apply(analysis.abs$ensmean[pos,],2,mean) }
      c=c+1
    }
    latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% which(unique(analysis.abs$lat) < 57))]
    pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
    zmean <- zmean[pos3,]
    sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
    dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
    gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
    gridmin1 <- latlist[gridminpos1]
    gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
    gridmin2 <- latlist[gridminpos2]
    zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
    sj <- sj * -1
    latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                    (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
    latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                    (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
    if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
    if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
    sj_slp <- c(latmin1,latmin2)    
       
    # stratospheric polar vortex
    pv_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PV1.calc"),]
    pv_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PV2.calc"),]
    pv <- pv_reg1 - pv_reg2 
    #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
    #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
    
    # Pacific walker circulation
    #Hind[10,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
    #                                                    which(analysis.abs$names=="omega500")] 
    pwc_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PWC1.calc"),]
    pwc_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PWC2.calc"),]
    pwc <- pwc_reg1 - pwc_reg2 
    # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
    # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
    
    # Dynamic Indian Monsoon index
    dimi_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="DIMI1.calc"),]
    dimi_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="DIMI2.calc"),]
    dimi <- dimi_reg1 - dimi_reg2 
    #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
    #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
    
    
    # NAO based on anomalies
    nao_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="NAO1.calc"),]
    nao_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="NAO2.calc"),]
    nao <- nao_reg1 - nao_reg2 
    
    # PNA based on anomalies
    pna_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA1.calc"),]
    pna_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA2.calc"),]
    pna_reg3 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA3.calc"),]
    pna_reg4 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA4.calc"),]
    pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
    # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 
    
    aind <- list(ensmean=matrix(0,nrow=length(indices),ncol=2), names=indices)
    aind$ensmean[which(indices=='ENH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ENH.temp2'),]
    aind$ensmean[which(indices=='NAM.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NAM.temp2'),]
    aind$ensmean[which(indices=='SAM.temp2'),] <- aind_tmp$ensmean[which(aind$names=='SAM.temp2'),]
    aind$ensmean[which(indices=='AFR.temp2'),] <- aind_tmp$ensmean[which(aind$names=='AFR.temp2'),]
    aind$ensmean[which(indices=='ASI.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ASI.temp2'),]
    aind$ensmean[which(indices=='AUS.temp2'),] <- aind_tmp$ensmean[which(aind$names=='AUS.temp2'),]
    aind$ensmean[which(indices=='ARC.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ARC.temp2'),]
    aind$ensmean[which(indices=='ANT.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ANT.temp2'),]
    aind$ensmean[which(indices=='NEU.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NEU.temp2'),]
    aind$ensmean[which(indices=='MED.temp2'),] <- aind_tmp$ensmean[which(aind$names=='MED.temp2'),]
    aind$ensmean[which(indices=='GLO.temp2'),] <- aind_tmp$ensmean[which(aind$names=='GLO.temp2'),]
    aind$ensmean[which(indices=='NAM.precip'),] <- aind_tmp$ensmean[which(aind$names=='NAM.precip'),]
    aind$ensmean[which(indices=='SAM.precip'),] <- aind_tmp$ensmean[which(aind$names=='SAM.precip'),]
    aind$ensmean[which(indices=='AFR.precip'),] <- aind_tmp$ensmean[which(aind$names=='AFR.precip'),]
    aind$ensmean[which(indices=='ASI.precip'),] <- aind_tmp$ensmean[which(aind$names=='ASI.precip'),]
    aind$ensmean[which(indices=='AUS.precip'),] <- aind_tmp$ensmean[which(aind$names=='AUS.precip'),]
    aind$ensmean[which(indices=='ARC.precip'),] <- aind_tmp$ensmean[which(aind$names=='ARC.precip'),]
    aind$ensmean[which(indices=='ANT.precip'),] <- aind_tmp$ensmean[which(aind$names=='ANT.precip'),]
    aind$ensmean[which(indices=='NEU.precip'),] <- aind_tmp$ensmean[which(aind$names=='NEU.precip'),]
    aind$ensmean[which(indices=='MED.precip'),] <- aind_tmp$ensmean[which(aind$names=='MED.precip'),]
    aind$ensmean[which(indices=='SH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='SH.temp2'),]
    aind$ensmean[which(indices=='NAM.slp'),] <- aind_tmp$ensmean[which(aind$names=='NAM.slp'),]
    aind$ensmean[which(indices=='SAM.slp'),] <- aind_tmp$ensmean[which(aind$names=='SAM.slp'),]
    aind$ensmean[which(indices=='AFR.slp'),] <- aind_tmp$ensmean[which(aind$names=='AFR.slp'),]
    aind$ensmean[which(indices=='ASI.slp'),] <- aind_tmp$ensmean[which(aind$names=='ASI.slp'),]
    aind$ensmean[which(indices=='AUS.slp'),] <- aind_tmp$ensmean[which(aind$names=='AUS.slp'),]
    aind$ensmean[which(indices=='ARC.slp'),] <- aind_tmp$ensmean[which(aind$names=='ARC.slp'),]
    aind$ensmean[which(indices=='ANT.slp'),] <- aind_tmp$ensmean[which(aind$names=='ANT.slp'),]
    aind$ensmean[which(indices=='NEU.slp'),] <- aind_tmp$ensmean[which(aind$names=='NEU.slp'),]
    aind$ensmean[which(indices=='MED.slp'),] <- aind_tmp$ensmean[which(aind$names=='MED.slp'),]
    aind$ensmean[which(indices=='NH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NH.temp2'),]
    aind$ensmean[which(indices=='EU.temp2'),] <- aind_tmp$ensmean[which(aind$names=='EU.temp2'),]
    aind$ensmean[which(indices=='EU.precip'),] <- aind_tmp$ensmean[which(aind$names=='EU.precip'),]
    aind$ensmean[which(indices=='EU.slp'),] <- aind_tmp$ensmean[which(aind$names=='EU.slp'),]
    aind$ensmean[which(indices=='HC.calc'),] <- hc #aind_tmp$ensmean[which(aind$names=='HC.calc'),]
    aind$ensmean[which(indices=='ITCZ.calc'),] <- itcz #aind_tmp$ensmean[which(aind$names=='ITCZ.calc'),]
    aind$ensmean[which(indices=='SJ_u200.calc'),] <- sj_u200 #aind_tmp$ensmean[which(aind$names=='SJ.calc'),]
    aind$ensmean[which(indices=='SJ_slp.calc'),] <- sj_slp
    aind$ensmean[which(indices=='PV.calc'),] <- pv #aind_tmp$ensmean[which(aind$names=='PV.calc'),]
    aind$ensmean[which(indices=='PWC.calc'),] <- pwc #aind_tmp$ensmean[which(aind$names=='PWC.calc'),]
    aind$ensmean[which(indices=='DIMI.calc'),] <- dimi #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),]
    aind$ensmean[which(indices=='NAO.calc'),] <- nao #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),] #c(0,0)
    aind$ensmean[which(indices=='PNA.calc'),] <- pna #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),] #c(0,0)
    aind$ensmean[which(indices=='DIMI'),] <- ana_ind$ensmean[which(ana_ind$names=='DIMI'),]
    aind$ensmean[which(indices=='HC'),] <- ana_ind$ensmean[which(ana_ind$names=='HC'),]
    aind$ensmean[which(indices=='SJ'),] <- ana_ind$ensmean[which(ana_ind$names=='SJ'),]
    aind$ensmean[which(indices=='Z100'),] <- ana_ind$ensmean[which(ana_ind$names=='z100'),]
    aind$ensmean[which(indices=='Z300'),] <- ana_ind$ensmean[which(ana_ind$names=='z300'),]
    aind$ensmean[which(indices=='PWC'),] <- ana_ind$ensmean[which(ana_ind$names=='PWC'),]
    
   

# BEM
aind_tmp <- list(bem=Hind%*%analysis.abs$bem, names=indices_tmp)

# Hadley cell strength
#    Hind[5,which(analysis.abs$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
#                                              which(analysis.abs$names=="stream")] 
# Hind[5,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
#                                                   which(analysis.abs$names=="omega500")] 
c=1
zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
for (l in unique(analysis.abs$lat)) {
  pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(analysis.abs$names=='omega500')
  pos <- intersect(pos1,pos2)
  if (length(analysis.abs$bem[pos,])==2) { zmean[c,] <- analysis.abs$bem[pos,] }
  if (length(analysis.abs$bem[pos,])>2) { zmean[c,] <- apply(analysis.abs$bem[pos,],2,mean) }
  #      zmean[c,] <- apply(analysis.abs$bem[pos,],2,mean)
  c=c+1
}
latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                            which(unique(analysis.abs$lat) < 30))]
pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
# cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f

# Hadley cell poleward extend
# Hind[6,which(analysis.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
#                                               which(analysis.abs$names=="u200")] 

# ITCZ location
#Hind[7,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
#                                                   which(analysis.abs$names=="omega500")] 
c=1
zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
for (l in unique(analysis.abs$lat)) {
  pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(analysis.abs$names=='omega500')
  pos <- intersect(pos1,pos2)
  if (length(analysis.abs$bem[pos,])==2) { zmean[c,] <- analysis.abs$bem[pos,] }
  if (length(analysis.abs$bem[pos,])>2) { zmean[c,] <- apply(analysis.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                            which(unique(analysis.abs$lat) < 30))]
pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(itcz[1]-min(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(itcz[2]-min(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
itcz <- c(latmin1,latmin2)

# subtropical jet (u200)
#Hind[8,which(analysis.abs$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
#                                               which(analysis.abs$names=="u200")] 
c=1
zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
for (l in unique(analysis.abs$lat)) {
  pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(analysis.abs$names=='u200')
  pos <- intersect(pos1,pos2)
  if (length(analysis.abs$bem[pos,])==2) { zmean[c,] <- analysis.abs$bem[pos,] }
  if (length(analysis.abs$bem[pos,])>2) { zmean[c,] <- apply(analysis.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% which(unique(analysis.abs$lat) < 57))]
pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
zmean <- zmean * -1
sj <- sj * -1
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
sj_u200 <- c(latmin1,latmin2)
# cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f

# subtropical jet (slp)
c=1
zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
for (l in unique(analysis.abs$lat)) {
  pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
  pos2 <- which(analysis.abs$names=='slp')
  pos <- intersect(pos1,pos2)
  if (length(analysis.abs$bem[pos,])==2) { zmean[c,] <- analysis.abs$bem[pos,] }
  if (length(analysis.abs$bem[pos,])>2) { zmean[c,] <- apply(analysis.abs$bem[pos,],2,mean) }
  c=c+1
}
latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% which(unique(analysis.abs$lat) < 57))]
pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
zmean <- zmean[pos3,]
sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj[1],digits=5)) 
gridmin1 <- latlist[gridminpos1]
gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj[2],digits=5)) 
gridmin2 <- latlist[gridminpos2]
zmean <- zmean * -1
sj <- sj * -1
latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
sj_slp <- c(latmin1,latmin2)

# stratospheric polar vortex
pv_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PV1.calc"),]
pv_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PV2.calc"),]
pv <- pv_reg1 - pv_reg2 
#cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
#cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f

# Pacific walker circulation
#Hind[10,which(analysis.abs$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
#                                                    which(analysis.abs$names=="omega500")] 
pwc_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PWC1.calc"),]
pwc_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PWC2.calc"),]
pwc <- pwc_reg1 - pwc_reg2 
# cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
# cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f

# Dynamic Indian Monsoon index
dimi_reg1 <- aind_tmp$bem[which(aind_tmp$names=="DIMI1.calc"),]
dimi_reg2 <- aind_tmp$bem[which(aind_tmp$names=="DIMI2.calc"),]
dimi <- dimi_reg1 - dimi_reg2 
#cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
#cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f


# NAO based on anomalies
nao_reg1 <- aind_tmp$bem[which(aind_tmp$names=="NAO1.calc"),]
nao_reg2 <- aind_tmp$bem[which(aind_tmp$names=="NAO2.calc"),]
nao <- nao_reg1 - nao_reg2 

# PNA based on anomalies
pna_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PNA1.calc"),]
pna_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PNA2.calc"),]
pna_reg3 <- aind_tmp$bem[which(aind_tmp$names=="PNA3.calc"),]
pna_reg4 <- aind_tmp$bem[which(aind_tmp$names=="PNA4.calc"),]
pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
# PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 

aind$bem <- matrix(0,nrow=length(indices),ncol=2)
#aind <- list(bem=matrix(0,nrow=length(indices),ncol=2), names=indices)
aind$bem[which(indices=='ENH.temp2'),] <- aind_tmp$bem[which(aind$names=='ENH.temp2'),]
aind$bem[which(indices=='NAM.temp2'),] <- aind_tmp$bem[which(aind$names=='NAM.temp2'),]
aind$bem[which(indices=='SAM.temp2'),] <- aind_tmp$bem[which(aind$names=='SAM.temp2'),]
aind$bem[which(indices=='AFR.temp2'),] <- aind_tmp$bem[which(aind$names=='AFR.temp2'),]
aind$bem[which(indices=='ASI.temp2'),] <- aind_tmp$bem[which(aind$names=='ASI.temp2'),]
aind$bem[which(indices=='AUS.temp2'),] <- aind_tmp$bem[which(aind$names=='AUS.temp2'),]
aind$bem[which(indices=='ARC.temp2'),] <- aind_tmp$bem[which(aind$names=='ARC.temp2'),]
aind$bem[which(indices=='ANT.temp2'),] <- aind_tmp$bem[which(aind$names=='ANT.temp2'),]
aind$bem[which(indices=='NEU.temp2'),] <- aind_tmp$bem[which(aind$names=='NEU.temp2'),]
aind$bem[which(indices=='MED.temp2'),] <- aind_tmp$bem[which(aind$names=='MED.temp2'),]
aind$bem[which(indices=='GLO.temp2'),] <- aind_tmp$bem[which(aind$names=='GLO.temp2'),]
aind$bem[which(indices=='NAM.precip'),] <- aind_tmp$bem[which(aind$names=='NAM.precip'),]
aind$bem[which(indices=='SAM.precip'),] <- aind_tmp$bem[which(aind$names=='SAM.precip'),]
aind$bem[which(indices=='AFR.precip'),] <- aind_tmp$bem[which(aind$names=='AFR.precip'),]
aind$bem[which(indices=='ASI.precip'),] <- aind_tmp$bem[which(aind$names=='ASI.precip'),]
aind$bem[which(indices=='AUS.precip'),] <- aind_tmp$bem[which(aind$names=='AUS.precip'),]
aind$bem[which(indices=='ARC.precip'),] <- aind_tmp$bem[which(aind$names=='ARC.precip'),]
aind$bem[which(indices=='ANT.precip'),] <- aind_tmp$bem[which(aind$names=='ANT.precip'),]
aind$bem[which(indices=='NEU.precip'),] <- aind_tmp$bem[which(aind$names=='NEU.precip'),]
aind$bem[which(indices=='MED.precip'),] <- aind_tmp$bem[which(aind$names=='MED.precip'),]
aind$bem[which(indices=='SH.temp2'),] <- aind_tmp$bem[which(aind$names=='SH.temp2'),]
aind$bem[which(indices=='NAM.slp'),] <- aind_tmp$bem[which(aind$names=='NAM.slp'),]
aind$bem[which(indices=='SAM.slp'),] <- aind_tmp$bem[which(aind$names=='SAM.slp'),]
aind$bem[which(indices=='AFR.slp'),] <- aind_tmp$bem[which(aind$names=='AFR.slp'),]
aind$bem[which(indices=='ASI.slp'),] <- aind_tmp$bem[which(aind$names=='ASI.slp'),]
aind$bem[which(indices=='AUS.slp'),] <- aind_tmp$bem[which(aind$names=='AUS.slp'),]
aind$bem[which(indices=='ARC.slp'),] <- aind_tmp$bem[which(aind$names=='ARC.slp'),]
aind$bem[which(indices=='ANT.slp'),] <- aind_tmp$bem[which(aind$names=='ANT.slp'),]
aind$bem[which(indices=='NEU.slp'),] <- aind_tmp$bem[which(aind$names=='NEU.slp'),]
aind$bem[which(indices=='MED.slp'),] <- aind_tmp$bem[which(aind$names=='MED.slp'),]
aind$bem[which(indices=='NH.temp2'),] <- aind_tmp$bem[which(aind$names=='NH.temp2'),]
aind$bem[which(indices=='EU.temp2'),] <- aind_tmp$bem[which(aind$names=='EU.temp2'),]
aind$bem[which(indices=='EU.precip'),] <- aind_tmp$bem[which(aind$names=='EU.precip'),]
aind$bem[which(indices=='EU.slp'),] <- aind_tmp$bem[which(aind$names=='EU.slp'),]
aind$bem[which(indices=='HC.calc'),] <- hc #aind_tmp$bem[which(aind$names=='HC.calc'),]
aind$bem[which(indices=='ITCZ.calc'),] <- itcz #aind_tmp$bem[which(aind$names=='ITCZ.calc'),]
aind$bem[which(indices=='SJ_u200.calc'),] <- sj_u200 #aind_tmp$bem[which(aind$names=='SJ.calc'),]
aind$bem[which(indices=='SJ_slp.calc'),] <- sj_slp
aind$bem[which(indices=='PV.calc'),] <- pv #aind_tmp$bem[which(aind$names=='PV.calc'),]
aind$bem[which(indices=='PWC.calc'),] <- pwc #aind_tmp$bem[which(aind$names=='PWC.calc'),]
aind$bem[which(indices=='DIMI.calc'),] <- dimi #aind_tmp$bem[which(aind$names=='DIMI.calc'),]
aind$bem[which(indices=='NAO.calc'),] <- nao #aind_tmp$bem[which(aind$names=='DIMI.calc'),] #c(0,0)
aind$bem[which(indices=='PNA.calc'),] <- pna #aind_tmp$bem[which(aind$names=='DIMI.calc'),] #c(0,0)
aind$bem[which(indices=='DIMI'),] <- ana_ind$bem[which(ana_ind$names=='DIMI'),]
aind$bem[which(indices=='HC'),] <- ana_ind$bem[which(ana_ind$names=='HC'),]
aind$bem[which(indices=='SJ'),] <- ana_ind$bem[which(ana_ind$names=='SJ'),]
aind$bem[which(indices=='Z100'),] <- ana_ind$bem[which(ana_ind$names=='z100'),]
aind$bem[which(indices=='Z300'),] <- ana_ind$bem[which(ana_ind$names=='z300'),]
aind$bem[which(indices=='PWC'),] <- ana_ind$bem[which(ana_ind$names=='PWC'),]


    
    
    
    # data for ensemble member
    #    aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
    #                       length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
    #                       dim(analysis$data)[2:3]))
    aind_tmp <- list(data=array(Hind %*% array(analysis.abs$data, c(nrow(analysis.abs$data), 
                  length(analysis.abs$data)/nrow(analysis.abs$data))), c(nrow(Hind), 
                  dim(analysis.abs$data)[2:3])), names=indices_tmp)
    # Hadley cell strength
    hc <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
      for (l in unique(analysis.abs$lat)) {
        pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(analysis.abs$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis.abs$data[pos,,m])==2) { zmean[c,] <- analysis.abs$data[pos,,m] }
        if (length(analysis.abs$data[pos,,m])>2) { zmean[c,] <- apply(analysis.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                                      which(unique(analysis.abs$lat) < 30))]
      pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      hc[m,] <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
    }
    # Hadley cell poleward extend
    
    # ITCZ location
    itcz <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
      for (l in unique(analysis.abs$lat)) {
        pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(analysis.abs$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis.abs$data[pos,,m])==2) { zmean[c,] <- analysis.abs$data[pos,,m] }
        if (length(analysis.abs$data[pos,,m])>2) { zmean[c,] <- apply(analysis.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > -30) %in% 
                                                      which(unique(analysis.abs$lat) < 30))]
      pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      itcz_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(itcz_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(itcz_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                    (2*(abs(itcz_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                    (2*(abs(itcz_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz[m,] <- c(latmin1,latmin2)
    }
    
    # subtropical jet (u200)
    sj_u200 <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
      for (l in unique(analysis.abs$lat)) {
        pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(analysis.abs$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(analysis.abs$data[pos,,m])==2) { zmean[c,] <- analysis.abs$data[pos,,m] }
        if (length(analysis.abs$data[pos,,m])>2) { zmean[c,] <- apply(analysis.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% 
                                                      which(unique(analysis.abs$lat) < 57))]
      pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj_tmp <- sj_tmp * -1 
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                   (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                   (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200[m,] <- c(latmin1,latmin2)
    }
    
    # subtropical jet (slp)
    sj_slp <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis.abs$lat)),ncol=2)
      for (l in unique(analysis.abs$lat)) {
        pos1 <- which(trunc(analysis.abs$lat,digits=5)==trunc(l,digits=5))
        pos2 <- which(analysis.abs$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(analysis.abs$data[pos,,m])==2) { zmean[c,] <- analysis.abs$data[pos,,m] }
        if (length(analysis.abs$data[pos,,m])>2) { zmean[c,] <- apply(analysis.abs$data[pos,,m],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis.abs$lat)[which(which(unique(analysis.abs$lat) > 0) %in% 
                                                  which(unique(analysis.abs$lat) < 57))]
      pos3 <- match(round(latlist,digits=5),round(unique(analysis.abs$lat),digits=5)) 
      zmean <- zmean[pos3,]
      sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=5)==round(sj_tmp[1],digits=5)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=5)==round(sj_tmp[2],digits=5)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj_tmp <- sj_tmp * -1 
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp[m,] <- c(latmin1,latmin2)
    }
    
    # stratospheric polar vortex
    pv <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pv_reg1 <- aind_tmp$data[which(aind_tmp$names=="PV1.calc"),,m]
      pv_reg2 <- aind_tmp$data[which(aind_tmp$names=="PV2.calc"),,m]
      pv[m,] <- pv_reg1 - pv_reg2 
    }
    
    # Pacific walker circulation
    pwc <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pwc_reg1 <- aind_tmp$data[which(aind_tmp$names=="PWC1.calc"),,m]
      pwc_reg2 <- aind_tmp$data[which(aind_tmp$names=="PWC2.calc"),,m]
      pwc[m,] <- pwc_reg1 - pwc_reg2 
    }
    
    # Dynamic Indian Monsoon index
    dimi <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      dimi_reg1 <- aind_tmp$data[which(aind_tmp$names=="DIMI1.calc"),,m]
      dimi_reg2 <- aind_tmp$data[which(aind_tmp$names=="DIMI2.calc"),,m]
      dimi[m,] <- dimi_reg1 - dimi_reg2 
    }
    
    # NAO based on anomalies
    nao <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      nao_reg1 <- aind_tmp$data[which(aind_tmp$names=="NAO1.calc"),,m]
      nao_reg2 <- aind_tmp$data[which(aind_tmp$names=="NAO2.calc"),,m]
      nao[m,] <- nao_reg1 - nao_reg2 
    }
    
    # PNA based on anomalies
    pna <- matrix(0,nrow=nmem,ncol=2)
    for (m in 1:nmem) {
      pna_reg1 <- aind_tmp$data[which(aind_tmp$names=="PNA1.calc"),,m]
      pna_reg2 <- aind_tmp$data[which(aind_tmp$names=="PNA2.calc"),,m]
      pna_reg3 <- aind_tmp$data[which(aind_tmp$names=="PNA3.calc"),,m]
      pna_reg4 <- aind_tmp$data[which(aind_tmp$names=="PNA4.calc"),,m]
      pna[m,] <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
    }
    
    aind$data <- array(0,c(length(indices),2,nmem))
    for (m in 1:nmem) {
      aind$data[which(indices=='ENH.temp2'),,m] <- aind_tmp$data[which(aind$names=='ENH.temp2'),,m]
      aind$data[which(indices=='NAM.temp2'),,m] <- aind_tmp$data[which(aind$names=='NAM.temp2'),,m]
      aind$data[which(indices=='SAM.temp2'),,m] <- aind_tmp$data[which(aind$names=='SAM.temp2'),,m]
      aind$data[which(indices=='AFR.temp2'),,m] <- aind_tmp$data[which(aind$names=='AFR.temp2'),,m]
      aind$data[which(indices=='ASI.temp2'),,m] <- aind_tmp$data[which(aind$names=='ASI.temp2'),,m]
      aind$data[which(indices=='AUS.temp2'),,m] <- aind_tmp$data[which(aind$names=='AUS.temp2'),,m]
      aind$data[which(indices=='ARC.temp2'),,m] <- aind_tmp$data[which(aind$names=='ARC.temp2'),,m]
      aind$data[which(indices=='ANT.temp2'),,m] <- aind_tmp$data[which(aind$names=='ANT.temp2'),,m]
      aind$data[which(indices=='NEU.temp2'),,m] <- aind_tmp$data[which(aind$names=='NEU.temp2'),,m]
      aind$data[which(indices=='MED.temp2'),,m] <- aind_tmp$data[which(aind$names=='MED.temp2'),,m]
      aind$data[which(indices=='GLO.temp2'),,m] <- aind_tmp$data[which(aind$names=='GLO.temp2'),,m]
      aind$data[which(indices=='NAM.precip'),,m] <- aind_tmp$data[which(aind$names=='NAM.precip'),,m]
      aind$data[which(indices=='SAM.precip'),,m] <- aind_tmp$data[which(aind$names=='SAM.precip'),,m]
      aind$data[which(indices=='AFR.precip'),,m] <- aind_tmp$data[which(aind$names=='AFR.precip'),,m]
      aind$data[which(indices=='ASI.precip'),,m] <- aind_tmp$data[which(aind$names=='ASI.precip'),,m]
      aind$data[which(indices=='AUS.precip'),,m] <- aind_tmp$data[which(aind$names=='AUS.precip'),,m]
      aind$data[which(indices=='ARC.precip'),,m] <- aind_tmp$data[which(aind$names=='ARC.precip'),,m]
      aind$data[which(indices=='ANT.precip'),,m] <- aind_tmp$data[which(aind$names=='ANT.precip'),,m]
      aind$data[which(indices=='NEU.precip'),,m] <- aind_tmp$data[which(aind$names=='NEU.precip'),,m]
      aind$data[which(indices=='MED.precip'),,m] <- aind_tmp$data[which(aind$names=='MED.precip'),,m]
      aind$data[which(indices=='SH.temp2'),,m] <- aind_tmp$data[which(aind$names=='SH.temp2'),,m]
      aind$data[which(indices=='NAM.slp'),,m] <- aind_tmp$data[which(aind$names=='NAM.slp'),,m]
      aind$data[which(indices=='SAM.slp'),,m] <- aind_tmp$data[which(aind$names=='SAM.slp'),,m]
      aind$data[which(indices=='AFR.slp'),,m] <- aind_tmp$data[which(aind$names=='AFR.slp'),,m]
      aind$data[which(indices=='ASI.slp'),,m] <- aind_tmp$data[which(aind$names=='ASI.slp'),,m]
      aind$data[which(indices=='AUS.slp'),,m] <- aind_tmp$data[which(aind$names=='AUS.slp'),,m]
      aind$data[which(indices=='ARC.slp'),,m] <- aind_tmp$data[which(aind$names=='ARC.slp'),,m]
      aind$data[which(indices=='ANT.slp'),,m] <- aind_tmp$data[which(aind$names=='ANT.slp'),,m]
      aind$data[which(indices=='NEU.slp'),,m] <- aind_tmp$data[which(aind$names=='NEU.slp'),,m]
      aind$data[which(indices=='MED.slp'),,m] <- aind_tmp$data[which(aind$names=='MED.slp'),,m]    
      aind$data[which(indices=='NH.temp2'),,m] <- aind_tmp$data[which(aind$names=='NH.temp2'),,m]
      aind$data[which(indices=='EU.temp2'),,m] <- aind_tmp$data[which(aind$names=='EU.temp2'),,m]
      aind$data[which(indices=='EU.precip'),,m] <- aind_tmp$data[which(aind$names=='EU.precip'),,m]
      aind$data[which(indices=='EU.slp'),,m] <- aind_tmp$data[which(aind$names=='EU.slp'),,m]
      aind$data[which(indices=='HC.calc'),,m] <- hc[m,] #aind_tmp$data[which(aind$names=='HC.calc'),,m]
      aind$data[which(indices=='ITCZ.calc'),,m] <- itcz[m,] #aind_tmp$data[which(aind$names=='ITCZ.calc'),,m]
      aind$data[which(indices=='SJ_u200.calc'),,m] <- sj_u200[m,] #aind_tmp$data[which(aind$names=='SJ.calc'),,m]
      aind$data[which(indices=='SJ_slp.calc'),,m] <- sj_slp[m,]
      aind$data[which(indices=='PV.calc'),,m] <- pv[m,] #aind_tmp$data[which(aind$names=='PV.calc'),,m]
      aind$data[which(indices=='PWC.calc'),,m] <- pwc[m,] #aind_tmp$data[which(aind$names=='PWC.calc'),,m]
      aind$data[which(indices=='DIMI.calc'),,m] <- dimi[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m]
      aind$data[which(indices=='NAO.calc'),,m] <- nao[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m] #c(0,0)
      aind$data[which(indices=='PNA.calc'),,m] <- pna[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m] #c(0,0)
      aind$data[which(indices=='DIMI'),,m] <- ana_ind$data[which(ana_ind$names=='DIMI'),,m]
      aind$data[which(indices=='HC'),,m] <- ana_ind$data[which(ana_ind$names=='HC'),,m]
      aind$data[which(indices=='SJ'),,m] <- ana_ind$data[which(ana_ind$names=='SJ'),,m]
      aind$data[which(indices=='Z100'),,m] <- ana_ind$data[which(ana_ind$names=='z100'),,m]
      aind$data[which(indices=='Z300'),,m] <- ana_ind$data[which(ana_ind$names=='z300'),,m]
      aind$data[which(indices=='PWC'),,m] <- ana_ind$data[which(ana_ind$names=='PWC'),,m]
    }
    
    aind$data[which(indices=='DIMI'),,] <- ana_ind$data[which(ana_ind$names=='DIMI'),,]
    aind$data[which(indices=='HC'),,] <- ana_ind$data[which(ana_ind$names=='HC'),,]
    aind$data[which(indices=='SJ'),,] <- ana_ind$data[which(ana_ind$names=='SJ'),,]
    aind$data[which(indices=='Z100'),,] <- ana_ind$data[which(ana_ind$names=='z100'),,]
    aind$data[which(indices=='Z300'),,] <- ana_ind$data[which(ana_ind$names=='z300'),,]
    aind$data[which(indices=='PWC'),,] <- ana_ind$data[which(ana_ind$names=='PWC'),,]
    
#     aind <- list(ensmean=Hind%*%analysis$ensmean, names=indices)
#     aind$ensmean[which(indices=='DIMI'),] <- ana_ind$ensmean[which(ana_ind$names=='DIMI'),]
#     aind$ensmean[which(indices=='HC'),] <- ana_ind$ensmean[which(ana_ind$names=='HC'),]
#     aind$ensmean[which(indices=='SJ'),] <- ana_ind$ensmean[which(ana_ind$names=='SJ'),]
#     aind$ensmean[which(indices=='Z100'),] <- ana_ind$ensmean[which(ana_ind$names=='z100'),]
#     aind$ensmean[which(indices=='Z300'),] <- ana_ind$ensmean[which(ana_ind$names=='z300'),]
#     aind$ensmean[which(indices=='PWC'),] <- ana_ind$ensmean[which(ana_ind$names=='PWC'),]  
# # data for ensemble member
#     aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
#                        length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
#                        dim(analysis$data)[2:3]))
#     aind$data[which(indices=='DIMI'),,] <- ana_ind$data[which(ana_ind$names=='DIMI'),,]
#     aind$data[which(indices=='HC'),,] <- ana_ind$data[which(ana_ind$names=='HC'),,]
#     aind$data[which(indices=='SJ'),,] <- ana_ind$data[which(ana_ind$names=='SJ'),,]
#     aind$data[which(indices=='Z100'),,] <- ana_ind$data[which(ana_ind$names=='z100'),,]
#     aind$data[which(indices=='Z300'),,] <- ana_ind$data[which(ana_ind$names=='z300'),,]
#     aind$data[which(indices=='PWC'),,] <- ana_ind$data[which(ana_ind$names=='PWC'),,]







    if (vali) { 
      H.giorgi <- compute_giorgi_H_v2(giorgi, validate) #, numvar=3) # 3 vars temp, precip, slp
      indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                   'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                   'NEU.temp2','MED.temp2',
                   'GLO.temp','NAM.precip','SAM.precip','AFR.precip',
                   'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                   'NEU.precip','MED.precip',
                   'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                   'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                   'NEU.slp','MED.slp',
                   'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
#                   'HC.calc', 'ITCZ.calc', 'SJ.calc', 'PV.calc', 
#                   'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc',
                   'HC', 'SJ', 'Z100', 'Z300', 'PWC')
      Hind <- matrix(0,nrow=length(indices),ncol=nrow(validate$data))
      Hind[1,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                             which(validate$names=="temp2")]
      Hind[2,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                             which(validate$names=="temp2")]
      Hind[3,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                             which(validate$names=="temp2")]
      Hind[4,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                             which(validate$names=="temp2")]
      Hind[5,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                             which(validate$names=="temp2")]
      Hind[6,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                             which(validate$names=="temp2")]
      Hind[7,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                             which(validate$names=="temp2")]
      Hind[8,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                             which(validate$names=="temp2")]
      Hind[9,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                             which(validate$names=="temp2")]
      Hind[10,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                              which(validate$names=="temp2")]
      Hind[11,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                               which(validate$names=="temp2")]
      Hind[12,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                               which(validate$names=="precip")]
      Hind[13,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                               which(validate$names=="precip")]
      Hind[14,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                               which(validate$names=="precip")]
      Hind[15,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                               which(validate$names=="precip")]
      Hind[16,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                               which(validate$names=="precip")]
      Hind[17,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                               which(validate$names=="precip")]
      Hind[18,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                               which(validate$names=="precip")]
      Hind[19,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                               which(validate$names=="precip")]
      Hind[20,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                               which(validate$names=="precip")]
      Hind[21,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                            which(validate$names=="temp2")]
      Hind[22,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                            which(validate$names=="slp")]
      Hind[23,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                            which(validate$names=="slp")]
      Hind[24,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                            which(validate$names=="slp")]
      Hind[25,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                            which(validate$names=="slp")]
      Hind[26,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                            which(validate$names=="slp")]
      Hind[27,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                            which(validate$names=="slp")]
      Hind[28,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                            which(validate$names=="slp")]
      Hind[29,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                            which(validate$names=="slp")]
      Hind[30,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                            which(validate$names=="slp")]       
      Hind[31,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                              which(validate$names=="temp2")]
      Hind[32,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                              which(validate$names=="temp2")]
      Hind[33,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                               which(validate$names=="precip")]
      Hind[34,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                            which(validate$names=="slp")] 
#      # land sea mask to take care of missing validata in the ocean      
      vpos <- which(!is.na(validate$data[,1]))
      Hind2 <- Hind[,vpos]
      validatanona <- validate$data[vpos,]
      vind <- list(data=Hind2 %*% validatanona, names=indices)

      try(vind$data[35,] <- vali_ind$data[which(vali_ind$names=="ind_rec_hc"),])
      try(vind$data[36,] <- vali_ind$data[which(vali_ind$names=="ind_rec_sj"),])
      try(vind$data[37,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z100"),])
      try(vind$data[38,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z300"),])
      try(vind$data[39,] <- vali_ind$data[which(vali_ind$names=="ind_rec_pwc"),])
    }
#  }  # end sixmonstatevector loop

  # save annual file
  if (!land_only) {
    if (vali) {
      save(aind,eind,vind, file=paste0('../data/prepplot_season/indices_',cyr,'.Rdata'))
    } else {
      save(aind,eind, file=paste0('../data/prepplot_season/indices_',cyr,'.Rdata'))
    }
  } else {
    if (vali) {
      save(aind,eind,vind, file=paste0('../data/prepplot_season/indices_landonly_',cyr,'.Rdata'))
    } else {
      save(aind,eind, file=paste0('../data/prepplot_season/indices_landonly_',cyr,'.Rdata'))
    }
  }
}    # end year loop
}    # end save prepplot







if (mergetime) {
  for (cyr in syr:eyr){
    print(cyr)
    ptm1 <- proc.time()
    if ((cyr < 1751) | (cyr > 1979)) {
      vali=F                 # switch off prepplot if no vali data selected
    } else {
      vali=T
    }
    
  if (!land_only) {
    load(file=paste0('../data/prepplot_season/indices_',cyr,'.Rdata'))
  } else {
    load(file=paste0('../data/prepplot_season/indices_landonly_',cyr,'.Rdata'))
  }
  # merge timesteps 
  if (cyr == 1751) {
    vind.allts=vind
    vind.allts$time=cyr
  }
  if (cyr == syr) {
    aind.allts=aind
    eind.allts=eind
    aind.allts$time=cyr
    eind.allts$time=cyr
  } else {
    aind.allts$data=abind(aind.allts$data,aind$data,along=2)
    aind.allts$ensmean=cbind(aind.allts$ensmean,aind$ensmean)
    aind.allts$bem=cbind(aind.allts$bem,aind$bem)
    aind.allts$time=c(aind.allts$time,cyr)
    if (vali) {
      vind.allts$data=cbind(vind.allts$data,vind$data)
      vind.allts$ensmean=cbind(vind.allts$ensmean,vind$ensmean)
      if (cyr != 1751) {
        vind.allts$time=c(vind.allts$time,cyr)
      }
    }
    eind.allts$data=abind(eind.allts$data,eind$data,along=2)
    eind.allts$ensmean=cbind(eind.allts$ensmean,eind$ensmean)
    eind.allts$bem=cbind(eind.allts$bem,eind$bem)
    eind.allts$time=c(eind.allts$time,cyr)
  } 
}

# save file
if (!land_only) {  
#  if (vali) {
    save(aind.allts,eind.allts,vind.allts, file='../data/prepplot_season/indices_allts.Rdata')
#  } else {
#    save(aind.allts,eind.allts, file='../data/prepplot_season/indices_allts.Rdata')
#  } 
} else {
#  if (vali) {
    save(aind.allts,eind.allts,vind.allts, file='../data/prepplot_season/indices_allts_landonly.Rdata')
#  } else {
#    save(aind.allts,eind.allts, file='../data/prepplot_season/indices_allts_landonly.Rdata')
#  }
}    # end land only
} #end mergetime










if (load_prepplot) {
  if (land_only) {
    load(paste0(prepplotpath,'indices_allts_landonly.Rdata'))
  } else {
    load(paste0(prepplotpath,'indices_allts.Rdata'))
  }
  
  # calc global mean temp spread/sd pre post 1800
  e.glomean <- array(eind.allts$data[1,,],c(2,(dim(eind.allts$data)[2]/2),
                                            dim(eind.allts$data)[3]))
  a.glomean <- array(aind.allts$data[1,,],c(2,(dim(aind.allts$data)[2]/2),
    dim(aind.allts$data)[3]))
  if (every2grid) {
    load(file="../data/cru4_ens_sd_2ndgrid.Rdata")
  } else {
    load(file="../data/cru4_ens_sd.Rdata")  
  }
  obs.spread.temp <- cbind(as.vector(cru4_oct_apr$data),as.vector(cru4_may_sep$data))
  apply(obs.spread.temp,2,mean,na.rm=T)
  # winter
  mean(apply(e.glomean[1,,],1,sd)[1:200])
  mean(apply(a.glomean[1,,],1,sd)[1:200])
  mean(apply(a.glomean[1,,],1,sd)[201:400])
  #sommer
  mean(apply(e.glomean[2,,],1,sd)[1:200])
  mean(apply(a.glomean[2,,],1,sd)[1:200])
  mean(apply(a.glomean[2,,],1,sd)[201:400])
  pdf(file='../figures/nat_data_paper/glo_mean_ana_ens_spread.pdf',width=9,height=3.5)
    par(mfrow=c(1,2))
   #winter
    plot(aind.allts$time,apply(e.glomean[1,,],1,sd),col=rgb(0,0,0,10,maxColorValue=10),ty='l',
         ylim=c(0.05,0.5),ylab="Ens. std. dev.",xlab='Year',main='Oct. - Apr.')
    lines(aind.allts$time,apply(a.glomean[1,,],1,sd),ty='l',col=rgb(10,0,0,7,maxColorValue=10))
  # summer  
    plot(aind.allts$time,apply(e.glomean[1,,],1,sd),col=rgb(0,0,0,10,maxColorValue=10),ty='l',
         ylim=c(0.05,0.5),ylab="Ens. std. dev.",xlab='Year',main='May - Sep.')
    lines(aind.allts$time,apply(a.glomean[1,,],1,sd),ty='l',col=rgb(10,0,0,7,maxColorValue=10))
  dev.off()
  
  # plot global mean temp
  pdf(file='../figures/nat_data_paper/eNH_summer_temp_recon4.pdf',width=10,height=4)
  wlen <- 1
  linew=1
  linew1=1
  linew2=2
  par(mfrow=c(1,1))
  par(mai=c(0.5,0.8,0.4,0.2))
  fs <- 1.0
  scal=F
  plbem=F
  plmem=F
  attributes(giorgi)
  anomper=(1901:1980)

  noaa <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',skip=3,header=F,sep='\t',
                     stringsAsFactors=FALSE)
  colnames(noaa) <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',nrow=1,sep='\t',
                               stringsAsFactors=FALSE)
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
                     lw=linew,sca=scal,xa='s',plotbem=F,plotmem=T,plotech=F,
                     title='Extratropical Northern Hemisphere Summer Temperature Anomalies')
  ts_e <- ts(data$allind[,"eindmean"],start=data$period[1])
  ts_a <- ts(data$allind[,"aindmean"],start=data$period[1])

  # northern hemisphere extra-tropics (>20N)
  nc=open.nc('../comparison_data/CRUTEM4_A-Smean_ENHmean.nc', write=F)
  print.nc(nc)
  tmp1 <- var.get.nc(nc, "TEMPERATURE_ANOMALY") # for CRU temp
  tlist <- var.get.nc(nc, "TIME")
  unitstr <- "days since 1850-01-15 00:00:00"
  yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
  pos <- yrv %in% data[['period']]
  tmp1 <- tmp1[pos]
  yrv <- yrv[pos]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp1[anopos])
  scalefaccru <- sd(tmp1[anopos])
  ts_i <- ts(tmp1,yrv[1])
  if (wlen > 1) {
    tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
  } else {
    tmp1 <- scale(tmp1,center=centerfac,scale=F)
  }
  lines(yrv[30:length(yrv)],tmp1[30:length(tmp1)],ty='l',lwd=linew,lty=1,
        col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  
  cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
  tmp2 <- cr[,2]
  yrv <- cr[,1]
  pos <- yrv %in% data[['period']]
  tmp2 <- tmp2[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata <- data[['allind']][pos,]
  data[['allind']][,1]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp2[anopos])
  ts_crow <- ts(tmp2,yrv[1])
  if (wlen > 1) {
    tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
  } else {
    tmp2 <- scale(tmp2,center=centerfac,scale=F)
  }
  lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), ylab="",xlab="")

  tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
  yrv <- noaa[,1]
  pos <- yrv %in% data[['period']]
  tmp4 <- tmp4[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata4 <- data[['allind']][pos,]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp4[anopos])
  ts_darr <- ts(tmp4,yrv[1])
  if (wlen > 1) {
    tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
  } else {
    tmp4 <- scale(tmp4,center=centerfac,scale=F)
  }
  lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")
  legend('topleft', c(
         #'ECHAM ens. mean','ECHAM spread',
         'Analysis ens. mean','Analysis spread',
         'CRUTEM4','Crowley et al. 2014','D Arrigo et al. 2006'),
         lwd=c(
          # 2,12,
          2,12,2,2,2),lty=c(1,1,1,1,1,1,1), col=c(
          #'black','lightgrey',
         'red','indianred1','blue','darkcyan','cyan'),cex=fs, bty='o', bg='white',
         box.col='white')
  
  # mark warmest/coldest decade before 1800
#  polygon(c(1790,1800,1800,1790),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(5,0,0,3,maxColorValue=10))
#  polygon(c(1830,1840,1840,1830),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(0,0,5,3,maxColorValue=10))
#  polygon(c(1640,1650,1650,1640),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(0,0,5,3,maxColorValue=10))
# calc corr, RE between series
  ts_all <- ts.union(ts_e,ts_a,ts_i,ts_crow,ts_darr)
  cor(ts_all,use="pairwise.complete.obs",method='spearman')
  cor(ts_all,use="pairwise.complete.obs",method='pearson')
#  ts_all_detr <- detrend(na.omit(ts_all))) #,bp=150)
  # look at high freq. changes only:
  ts_all_detr <- ts_all-(filter(ts_all,rep((1/11),11),method="convolution",sides=2))
  cor(ts_all_detr,use="pairwise.complete.obs",method='spearman')
  cor(ts_all_detr,use="pairwise.complete.obs",method='pearson')
#  1- sum((ts_all[,2]-ts_all[,1])**2,na.rm=T)  / sum((ts_all[,1]-mean(ts_all[,1]))**2,na.rm=T)
  library("dcv")
  nona_i <- window(ts_all,1850,2004)
  nona_p <- window(ts_all,1782,1984)
  nona_pi <- window(ts_all,1850,1984)
  # D'Arrigo 1608-1990
#  test.RE(nona_i[,2],nona_i[,1])
#  test.RE(nona_p[,3],nona_p[,1])
#  test.RE(nona_pi[,2],nona_pi[,3])
  1-((sum(nona_i[,2]-nona_i[,3])^2)/(sum(nona_i[,1]-nona_i[,3])^2)) # RE CRU
  1-((sum(nona_p[,2]-nona_p[,4])^2)/(sum(nona_p[,1]-nona_p[,4])^2)) # RE Crowley
  1-((sum(nona_pi[,2]-nona_pi[,4])^2)/(sum(nona_pi[,1]-nona_pi[,4])^2)) # RE Crowley
  1-((sum(nona_p[,2]-nona_p[,5])^2)/(sum(nona_p[,1]-nona_p[,4])^2)) # RE D'Arrigo
  1-((sum(nona_pi[,2]-nona_pi[,5])^2)/(sum(nona_pi[,1]-nona_pi[,4])^2)) # RE D'Arrigo
dev.off() # end Fig. 1 current version with 1 panel
    
     
  
  
  
  
  
  
  
  
     
# plot difference between seasons, ENH, NH, Glob
  reg3 <- which(eind.allts$names=="GLO.temp2")
  reg2 <- which(eind.allts$names=="NH.temp2")
  reg <- which(eind.allts$names=="ENH.temp2")
  period <- eind.allts$time
  anomper=(1901:1980)
# annual mean
  years <- rep(eind.allts$time,each=2)  
  eindmeanann <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
  aindmeanann <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
  eindmeanannNH <- aggregate(eind.allts$ensmean[reg2,],list(years),mean,na.rm=T)[,2]
  aindmeanannNH <- aggregate(aind.allts$ensmean[reg2,],list(years),mean,na.rm=T)[,2]
  eindmeanannGL <- aggregate(eind.allts$ensmean[reg3,],list(years),mean,na.rm=T)[,2]
  aindmeanannGL <- aggregate(aind.allts$ensmean[reg3,],list(years),mean,na.rm=T)[,2]
  # summer
  somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
  somcolv <- is.even(seq(1,ncol(vind.allts$data)))    
  eindmeansum <- eind.allts$ensmean[reg,somcol]  
  aindmeansum <- aind.allts$ensmean[reg,somcol]  
# winter
  wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
  wincolv <- is.odd(seq(1,ncol(vind.allts$data)))
  eindmeanwin <- eind.allts$ensmean[reg,wincol]  
  aindmeanwin <- aind.allts$ensmean[reg,wincol]  
# join all indices
  allind <- cbind(eindmeanannGL,aindmeanannGL,eindmeanannNH,aindmeanannNH,
              eindmeanann,aindmeanann,eindmeansum,aindmeansum,eindmeanwin,aindmeanwin)
  anopos <- period %in% anomper
  centerfac <- apply(allind[anopos,],2,mean)
  centerfac[2] <- centerfac[1] 
  centerfac[3] <- centerfac[4]
  centerfac[5] <- centerfac[6]
  allind <- scale(allind,center=centerfac,scale=F)  
  plot(period,allind[,2],ty='l',col=rgb(0,0,0,0,maxColorValue=10), ylab="Temperature Anomaly",
       ylim=c(min(allind[,2],na.rm=T)-0.5,max(allind[,2],na.rm=T)+0.5),
       main='',xaxt='s',bty='n')
  lines(period,allind[,5],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
  lines(period,allind[,7],ty='l',lwd=1,lty=1,col=rgb(10,0,0,8,maxColorValue=10), ylab="",xlab="")
  lines(period,allind[,9],ty='l',lwd=1,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
  cor(allind[,c(6,8,10)])
  mean(allind[1:100,7]-allind[1:100,9]) # sum vs win diff
  mean(allind[1:100,5]-allind[1:100,7]) # ann vs sum
#  mean(allind[1:100,8]-allind[1:100,10]) # sum vs win diff
#  mean(allind[1:100,6]-allind[1:100,8]) # ann vs sum
  legend('topleft', c('ECHAM ens. mean annual mean','ECHAM ens. mean Winter',
    'ECHAM ens. mean Summer'),lwd=c(2,2,2),lty=c(1,1,1), 
    col=c('black','blue','red'),cex=fs, bty='o', bg='white',box.col='white')
  
# find and mark coldest/warmest decade
  d <- aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
                     (is.even(seq(1,ncol(aind.allts$ensmean))))]
  dec <- rep(NA,250)
  dectime <- rep(NA,250)
  for (t in 1:length(dec)) {
    if (t==1) {
      dec[t] <- mean(d[1:5])
      dectime[t] <- 1603
      i=2
    } else {
      dec[t] <- mean(d[i:(i+4)])
      dectime[t] <- dectime[t-1] + 1 
      i=i+1
    }
  } 
#  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec,dectime)
  tmp[order(dec),]
  
# mark warmest/coldest decade before 1800
  polygon(c(1790,1800,1800,1790),c(0.9,0.9,1,1),density=NA,col=rgb(5,0,0,3,maxColorValue=10))
  polygon(c(1640,1650,1650,1640),c(0.9,0.9,1,1),density=NA,col=rgb(0,0,5,3,maxColorValue=10))

#   runm <- runmean(aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
#                      (is.even(seq(1,ncol(aind.allts$ensmean))))],11)
#   sortrunm <- sort(runm)
#   aind.allts$time[which(runm==sortrunm[1])]
#   aind.allts$time[which(runm==sortrunm[2])]
#   aind.allts$time[which(runm==sortrunm[3])]
 
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)])]
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)-1])]
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)-2])]
  dev.off()

  
  
#   
#   # ADD NAO plots + reconstructions + skill (corr and RE, ...?)
#   load('../data/prepplot_season/indices_allts.Rdata') # NAO and PNA only in seasonal! why?
#   # NAO
#   data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),
#                      lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,pag=F,plotech=T,title='NAO')
#   # some ocean based indices only NOT in landonly version
#   stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#                         header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
#                         seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#   head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#                      nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
#                                 seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#   newhead <- rep(NA,ncol(head))
#   for (i in 1:ncol(head)) {
#     newhead[i] <- paste(head[1,i],head[2,i],sep='_')
#   }
#   colnames(stefind) <- newhead
#   vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
#   yrv <- stefind[,1]
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(0,0,5,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
#   ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
#                         use="pairwise.complete.obs"),digits=2) 
#   acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
#                         use="pairwise.complete.obs"),digits=2) 
#   ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],
#                        use="pairwise.complete.obs"),digits=2) 
#   acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],
#                        use="pairwise.complete.obs"),digits=2)
#   legend('bottomleft', paste('ECHAM(1,2)/EKF(3,4) ens. mean - 20CR/NNR:',
#     ecormean[1],ecormean[3],acormean[1],acormean[3]),
# #    paste('ECHAM/EKF BEM - 20CR:',ecorbem[1],acorbem[1])),
#     cex=fs, bty='o', bg='white', box.col='white')
#   
# #  legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',
# #    ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',
# #    acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',
# #    ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',
# #    acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) 
#   nao_trou <- read.table("/Users/joerg/Documents/climdata/indices/www/nao-trouet2009.txt",skip=93,header=T)
#   nao_cook <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_cook2002.txt",skip=50,header=F)
#   nao_lut <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_sea_luterbacher_2002.txt",
#                         skip=19,header=T)
#   nao_lut2 <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_mon_luterbacher_2002.txt",
#                         skip=27,header=T)
#   pna_trou <- read.table("/Users/joerg/Documents/climdata/indices/www/pna-trouet2009.txt",skip=75,header=F)
#   
#   # data has to be loaded just before "plot_pages"! why?
#   naol <- nao_lut[nao_lut[,2]=="wi",]
#   naol <- ts(naol[naol[,1]>=1600,3],start=1600)
#   agglist <- c(rep(1659,3),rep(seq(1660,2001),each=3))
#   naol2 <- nao_lut2[(nao_lut2[,2]=="Dec" | nao_lut2[,2]=="Jan" | nao_lut2[,2]=="Feb"),]
# #  naol2 <- nao_lut2[(nao_lut2[,2]=="Oct" | nao_lut2[,2]=="Nov" | nao_lut2[,2]=="Dec" | 
# #                     nao_lut2[,2]=="Jan" | nao_lut2[,2]=="Feb" | nao_lut2[,2]=="Mar"),]
# #  agglist <- c(rep(1659,4),rep(seq(1660,2001),each=6))
#   naol2agg <- aggregate(naol2,list(agglist),mean)
#   naol2 <- ts(naol2agg[,4],start=1659,freq=1)
#   naol3 <- ts(c(naol,naol2),start=1600)
#   naot <- ts(nao_trou[nao_trou[,1]>=1600,2],start=1600)
#   naoc <- ts(nao_cook[nao_cook[,1]>=1600,2],start=1600)
#   vdata <- ts.union(naol3,naot,naoc)
#   yrv <- seq(1600,1600+nrow(vdata)-1)
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
# #  ecormean <- round(cor(data[['allind']][,1],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  acormean <- round(cor(data[['allind']][,2],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  ecorbem <- round(cor(data[['allind']][,3],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  acorbem <- round(cor(data[['allind']][,4],vdata[,],use="pairwise.complete.obs"),digits=2)
#   round(cor(vdata,use="pairwise.complete.obs"),digits=2)
#   round(cor(data[['allind']][,c(1,2)],vdata[,],use="pairwise.complete.obs"),digits=2)
# #  legend('topleft', c(paste(''Analysis ens. mean - Luterb., Trouet, Cook:',
# #    ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens.mean - Lut.,Trouet,Cook:',
# #    acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM BEM - Lut.,Trouet,Cook:',
# #    ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - Lut.,Trouet,Cook:',
# #    acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   
#   # PNA
#   data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#   #abline(v=vf,col='black',lty=2)
#   vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
#   yrv <- stefind[,1]
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
#   ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
#   legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   dev.off()
  
  
  
  

  
  
  
  
  
#  } else {
#    legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread'),
#           lwd=c(2,12,2,12),lty=c(1,1,1,1), col=c('black','lightgrey',
#           'red','indianred1'),cex=fs, bty='o',bg='white', box.col='white')  
#  }
  #data <- plot_pages(wl=wlen,reg=1,seas='sum')
  #data <- plot_pages(wl=wlen,reg=1,seas='win')
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
                     lw=linew,sca=scal,xa='s',plotbem=F,plotmem=T,plotech=F,yl='')  
#  par(new=T)
#  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
#                   lw=linew,sca=scal,xa='n',plotbem=F,plotmem=F,plotech=F,yl='')  
#  par(new=T)  
#  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='win',anomper=(1901:1980),
#                   lw=linew,sca=scal,xa='n',plotbem=F,plotmem=F,plotech=F,yl='')  
  nc=open.nc('../comparison_data/cru3_tmp_25_annmean_enhmean.nc', write=F)
  print.nc(nc)
  tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
  tlist <- var.get.nc(nc, "TIME")
  unitstr <- "months since 1901-01-15 00:00:00"
  yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
  pos <- yrv %in% data[['period']]
  tmp1 <- tmp1[pos]
  yrv <- yrv[pos]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp1[anopos])
  scalefaccru <- sd(tmp1[anopos])
  if (wlen > 1) {
    tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
  } else {
    tmp1 <- scale(tmp1,center=centerfac,scale=F)
  }
  lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  
  cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
  tmp2 <- cr[,2]
  yrv <- cr[,1]
  pos <- yrv %in% data[['period']]
  tmp2 <- tmp2[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata <- data[['allind']][pos,]
  data[['allind']][,1]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp2[anopos])
  if (wlen > 1) {
    tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
  } else {
    tmp2 <- scale(tmp2,center=centerfac,scale=F)
  }
  lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
  
#   tmp3 <- noaa[,which(colnames(noaa)=="esper2002")]
#   yrv <- noaa[,1]
#   pos <- yrv %in% data[['period']]
#   tmp3 <- tmp3[pos]
#   yrv <- yrv[pos]
#   pos <- data[['period']] %in% yrv 
#   cordata3 <- data[['allind']][pos,]
#   anopos <- yrv %in% anomper
#   centerfac <- mean(tmp3[anopos])
#   if (wlen > 1) {
#     tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
#   } else {
#     tmp3 <- scale(tmp3,center=centerfac,scale=(scalefaccru*2))
#   }
#   lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")
  
  tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
  yrv <- noaa[,1]
  pos <- yrv %in% data[['period']]
  tmp4 <- tmp4[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata4 <- data[['allind']][pos,]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp4[anopos])
  if (wlen > 1) {
    tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
  } else {
    tmp4 <- scale(tmp4,center=centerfac,scale=F)
  }
  lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(0,2,0,5,maxColorValue=10), ylab="",xlab="")
  
  # tmp5 <- noaa[,which(colnames(noaa)=="wilson2007")]
  # yrv <- noaa[,1]
  # pos <- yrv %in% data[['period']]
  # tmp5 <- tmp5[pos]
  # yrv <- yrv[pos]
  # anopos <- yrv %in% anomper
  # centerfac <- mean(tmp5[anopos])
  # if (wlen > 1) {
  #   tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
  # } else {
  #   tmp5 <- scale(tmp5,center=centerfac,scale=F)
  # }
  # lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")
  
  #ecormean <- round(cor(data[['allind']][,1],tmp2,use="pairwise.complete.obs"),digits=2) 
  #acormean <- round(cor(data[['allind']][,2],tmp2,use="pairwise.complete.obs"),digits=2) 
  ecormean <- round(cor(cordata[,1],tmp2,use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(cordata[,2],tmp2,use="pairwise.complete.obs"),digits=2) 
  #ecormean3 <- round(cor(cordata3[,1],tmp3,use="pairwise.complete.obs"),digits=2) 
  #acormean3 <- round(cor(cordata3[,2],tmp3,use="pairwise.complete.obs"),digits=2) 
  ecormean4 <- round(cor(cordata4[,1],tmp4,use="pairwise.complete.obs"),digits=2) 
  acormean4 <- round(cor(cordata4[,2],tmp4,use="pairwise.complete.obs"),digits=2) 
  ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
  acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - Crowley 2014:',ecormean,'/',acormean), 
#         paste('ECHAM/Analysis ens. mean - Esper 2002:',ecormean3,'/',acormean3), 
         paste('ECHAM/Analysis ens. mean - DArrigo 2007:',ecormean4,'/',acormean4), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
  legend('bottomright', c('Analysis ens. mean','Analysis spread','CRU TS3','Crowley 2014',
         'D Arrigo 2006'),lwd=c(1,12,2,2,2,2),lty=c(1,1,1,1,1,1), 
         col=c('red','indianred1','blue','cyan','darkgreen'),cex=fs, bty='o', bg='white', box.col='white')  
#legend('bottomright', c('Analysis ens. mean','Analysis spread','CRU TS3','Crowley 2014',
#                        'Esper 2002','D Arrigo 2006'),lwd=c(1,12,2,2,2,2),lty=c(1,1,1,1,1,1), 
#       col=c('red','indianred1','blue','cyan','violet','orange'),cex=fs, bty='o', bg='white',
#       box.col='white')  


  
  # Europe
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='EU.temp2'),lw=linew,sca=scal, yl='Temperature',xa='s',
                     plotbem=plbem,plotmem=plmem,plotech=F,seas='sum')
  reg <- which(vind.allts$names=='EU.temp2')
  years <- rep(vind.allts$time,each=2) 
  years <- years[(length(years)-155):length(years)]
  vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
  vindmean <- aggregate(vdat,list(years),mean)[,2]
  mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
  ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
  polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
          c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
          density=NA,col=rgb(0,5,5,3,maxColorValue=10))
  #lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
  #      ylab="",xlab="")
  lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
        ylab="",xlab="")
  lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
        lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  lines(data[['periodv']],data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,8,maxColorValue=10),
        ylab="",xlab="")
  ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
  ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
  ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - PAGES recon.:',ecormean,'/',acormean), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
#  legend('bottomright', c('Analysis ens. mean','PAGES','PAGES spread','CRU TS3'),
#         lwd=c(2,2,12,2,2),lty=c(1,1,1,1,1), 
#         col=c('red','lightcyan','darkcyan','blue'),cex=fs, bty='o', bg='white', box.col='white')  
#
#    legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
#                        paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
#           cex=fs, bty='o', bg='white', box.col='white')   
  
  # North America
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAM.temp2'),lw=linew,sca=scal, yl='',xa='s',
                     plotbem=plbem,plotmem=plmem,plotech=F,seas='sum')
  reg <- which(vind.allts$names=='NAM.temp2')
  years <- rep(vind.allts$time,each=2) 
  years <- years[(length(years)-155):length(years)]
  vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
  vindmean <- aggregate(vdat,list(years),mean)[,2]
  mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
  ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
  polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
          c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
          density=NA,col=rgb(0,5,5,3,maxColorValue=10))
  #lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
  #      ylab="",xlab="")
  lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
        ylab="",xlab="")
  lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
        lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  lines(data[['periodv']],data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,8,maxColorValue=10), 
        ylab="",xlab="")
  ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
  ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
  ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - PAGES recon.:',ecormean,'/',acormean), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
  legend('bottomright', c('Analysis ens. mean','PAGES','PAGES spread','CRU TS3'),
         lwd=c(2,2,12,2,2),lty=c(1,1,1,1,1), 
         col=c('red','lightcyan3','darkcyan','blue'),cex=fs, bty='o', bg='white', box.col='white') 
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# continental temp comparison
pdf(file='../figures/pages_cont_temp_comp.pdf',width=10,height=10)
wlen <- 1
linew=1
linew1=1
linew2=2
par(mfrow=c(4,2))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
scal=F
plbem=F
plmem=F
attributes(giorgi)



# Arctic
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ARC.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='ARC.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,22],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,23],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,23],use="pairwise.complete.obs"),digits=2)
acormean <- round(cor(data[['allind']][,2],data[['pages']][,23],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,23],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,23],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Europe
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='EU.temp2'),lw=linew,sca=scal, yl='',xa='n',
                   plotbem=plbem,plotmem=plmem,plotech=F,seas='yrmean')
reg <- which(vind.allts$names=='EU.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
        c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
        density=NA,col=rgb(3,3,3,2,maxColorValue=10))
#lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
#      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(3,3,3,7,maxColorValue=10), 
      ylab="",xlab="")
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(period,data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,5,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# North America
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAM.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='NAM.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,6],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,7],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,7],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Asia
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ASI.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='ASI.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,26],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,27],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,27],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}
  
# Africa
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='AFR.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='AFR.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - CRU TS3:',ecorval),
                    paste('Analysis ens. mean - CRU TS3:',acorval),
                    paste('ECHAM BEM - CRU TS3:',ecorvalbem),
                    paste('Analysis BEM - CRU TS3:',acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# South America
data <- plot_pages(wl=wlen,reg=which(attributes(giorgi)$names=='SAM'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='SAM.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,18],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,19],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,19],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Australia
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='AUS.temp2'),lw=linew,sca=scal, yl='',xa='s',plotbem=plbem)
reg <- which(vind.allts$names=='AUS.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,10],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,12],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,12],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Antarctica
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ANT.temp2'),box='n',lw=linew,sca=scal, yl='',xa='s',plotbem=plbem)
#reg <- which(vind.allts$names=='ANT.temp2')
#years <- rep(vind.allts$time,each=2) 
#years <- years[(length(years)-155):length(years)]
#vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
#vindmean <- aggregate(vdat,list(years),mean)[,2]
#lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
#lwd=linew, lty=1, col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,15],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,15],use="pairwise.complete.obs"),digits=2)
#ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
#acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
#ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
#acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon:',ecormean),
                    paste('Analysis ens. mean - PAGES recon:',acormean),
                    paste('ECHAM BEM - PAGES recon:',ecorbem),
                    paste('Analysis BEM - PAGES recon:',acorbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread',
#'Best ECHAM ensemble member','Best Analysis ensemble member',
'CRU TS3','PAGES reconstruction','Other reconstruction'),lwd=c(2,12,2,12,2,2,2),lty=c(1,1,1,1,1,2,2), col=c('black','lightgrey','red','indianred1',
#'grey','darkorange',                                                                                                                                                                                                                                                                              
'blue','cyan','violet'),cex=fs, bty='o', bg='white', box.col='white')

dev.off()












# plot global and hemisperic averages
pdf(file='../figures/glo_hemis_temp_comp.pdf',width=5,height=10)
wlen <- 1
linew <- 1
par(mfrow=c(4,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
scal=F
plbem=F
#eind.allts$names
noaa <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',skip=3,header=F,sep='\t',stringsAsFactors=FALSE)
colnames(noaa) <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',nrow=1,sep='\t',stringsAsFactors=FALSE)

# global mean land only
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='GLO.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_glomean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012global.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008e")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="oerlemans2005")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp3,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp3,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Mann 2008, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Mann 2008, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Leclercq 2012','Mann 2008','Oerlemans 2005'),lwd=c(2,2,2),lty=c(1,1,1), col=c('cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')





# northern hemisphere extra-tropics (>20N)
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem)
#data <- plot_pages(wl=wlen,reg=1,seas='sum')
#data <- plot_pages(wl=wlen,reg=1,seas='win')
nc=open.nc('../comparison_data/cru3_tmp_25_annmean_enhmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
tmp2 <- cr[,2]
yrv <- cr[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
pos <- data[['period']] %in%yrv 
cordata <- data[['allind']][pos,]
data[['allind']][,1]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="esper2002")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
cordata3 <- data[['allind']][pos,]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
cordata4 <- data[['allind']][pos,]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

# tmp5 <- noaa[,which(colnames(noaa)=="wilson2007")]
# yrv <- noaa[,1]
# pos <- yrv %in% data[['period']]
# tmp5 <- tmp5[pos]
# yrv <- yrv[pos]
# anopos <- yrv %in% anomper
# centerfac <- mean(tmp5[anopos])
# if (wlen > 1) {
#   tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
# } else {
#   tmp5 <- scale(tmp5,center=centerfac,scale=F)
# }
# lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")

#ecormean <- round(cor(data[['allind']][,1],tmp2,use="pairwise.complete.obs"),digits=2) 
#acormean <- round(cor(data[['allind']][,2],tmp2,use="pairwise.complete.obs"),digits=2) 
ecormean <- round(cor(cordata[,1],tmp2,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(cordata[,2],tmp2,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Crowley 2014, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Crowley 2014, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Crowley 2014','Esper 2002','D Arrigo 2006'),lwd=c(2,2,2),lty=c(1,1,1), col=c('cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')




# northern hemisphere complete land
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem,use18=F,pag=F)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_nhmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012nh.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008g")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="ammann2007")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

tmp5 <- noaa[,which(colnames(noaa)=="moberg2005")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp5 <- tmp5[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp5[anopos],na.rm=T)
if (wlen > 1) {
  tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
} else {
  tmp5 <- scale(tmp5,center=centerfac,scale=F)
}
lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp4,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp4,use="pairwise.complete.obs"),digits=2)
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Ammann 2007, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Ammann 2007, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Leclercq 2012','Mann 2008', 'Ammann 2007', 'Moberg 2005'),lwd=c(2,2,2,2),lty=c(1,1,1,1), col=c('cyan','violet','yellow','green'),cex=fs, bty='o', bg='white', box.col='white')




# southern hemisphere extra-tropical land ( < -20N )
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='s',plotbem=plbem, pag=F, use18=F, plotmem=T)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_shmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012sh.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008i")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

nk <- read.table('../comparison_data/SH_temp_recon_neukom2014.txt',header=T,skip=106,stringsAsFactors=FALSE)
tmp4 <- nk[,2]
yrv <- nk[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp4,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp4,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Neukom 2014, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Neukom 2014, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')
if (plbem) {
  legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread','Best ECHAM ensemble member','Best Analysis ensemble member','CRU TS3','Leclercq 2012','Mann 2008','Neukom 2014'),lwd=c(2,12,2,12,2,2,2,2,2,2),lty=c(1,1,1,1,2,2,1,1,1,1), col=c('black','lightgrey','red','indianred1','grey','darkorange','blue','cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')
} else {
  legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread','CRU TS3','Leclercq 2012','Mann 2008','Neukom 2014'),lwd=c(2,12,2,12,2,2,2,2),lty=c(1,1,1,1,1,1,1,1), col=c('black','lightgrey','red','indianred1','blue','cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')
}  

# load instr. and recon data for comparison
#system("cdo -fldmean -sellonlatbox,-180,180,20,90 /Users/joerg/Documents/climdata/cru/HadCRUT.4.2.0.0.median.nc ../comparison_data/crutem4_enh_mean.nc")
 
dev.off()







# PLOT atm indices
if (land_only) {
  load(paste0(prepplotpath,'indices_allts_landonly.Rdata'))
} else {
  load(paste0(prepplotpath,'indices_allts.Rdata'))
}
#load(paste0(prepplotpath,'indices_allts.Rdata')) 
# some ocean based indices only NOT in landonly version
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
# (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead
#load('../data/indices_recon.Rdata')
vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]

if (land_only) {
  pdf(file='../figures/indices/atm_index_comp_landonly.pdf',width=20,height=22)
} else {
  pdf(file='../figures/indices/atm_index_comp.pdf',width=20,height=22)
}  
wlen <- 1
linew <- 1
par(mfrow=c(9,2))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
#eind.allts$names

# HC
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='HC.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata1 <- stefind['HC_20C'] # [,'SJ_20C'] 
vdata2 <- stefind['HC_REC'] # [,'SJ_REC'] 
vdata3 <- stefind['HC_NCEP']
vdata4 <- stefind['HC_ERA-40']
vdata <- cbind(vdata1,vdata2,vdata3,vdata4)
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
#centerfac2 <- mean(vdata2,na.rm=T)
#scalefac2 <- sd(vdata2,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#  vdata1 <- runmean(scale(vdata1,center=centerfac1,scale=scalefac1),wlen)
#  vdata2 <- runmean(scale(vdata2,center=centerfac2,scale=scalefac2),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
#  vdata1 <- scale(vdata1,center=centerfac1,scale=scalefac1)
#  vdata2 <- scale(vdata2,center=centerfac2,scale=scalefac2)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='HC'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# SJ_u200
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ_u200.calc'),seas='win',anom=T,anomper=(1901:1980),lw=linew,sca=F,yl='',use18=T)
vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# SJ_slp
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ_slp.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='',use18=T)
vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# PWC
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PWC.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['PWC_20C'],stefind['PWC_REC'],stefind['PWC_NCEP'],stefind['PWC_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PWC'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# DIMI
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='DIMI.calc'),seas='sum',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['DIMI_20C'],stefind['DIMI_REC'],stefind['DIMI_NCEP'],stefind['DIMI_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='DIMI'),seas='sum',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# PV
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PV.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['Z100_20C'],stefind['Z100_REC'],stefind['Z100_NCEP'],stefind['Z100_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z100'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
#abline(v=vf,col='black',lty=2)
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# z300
#data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='z300.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T)
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z300'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['Z300_20C'],stefind['Z300_REC'],stefind['Z300_NCEP'],stefind['Z300_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# ITCZ
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ITCZ.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')


# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

dev.off()










#rm(list=ls())
#setwd('~/unibe/projects/EnSRF/r/src')
#source('EnSRF_functions.R')
# PLOT atm indices
if (land_only) {
  load(paste0(prepplotpath,'indices_allts_landonly.Rdata'))
} else {
  load(paste0(prepplotpath,'indices_allts.Rdata'))
}
#load(paste0(prepplotpath,'indices_allts.Rdata')) # some ocean based indices only NOT in landonly version
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
# (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead
#load('../data/indices_recon.Rdata')
vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]



if (land_only) {
  pdf(file='../figures/indices/atm_index_comp_short_landonly.pdf',width=10,height=10)
} else {
  pdf(file='../figures/indices/atm_index_comp_short.pdf',width=10,height=10)
}
wlen <- 1
linew <- 1
par(mfrow=c(3,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
#eind.allts$names

# PV
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PV.calc'),seas='win',anomper=(1901:1980),
                   lw=linew,sca=T,yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['Z100_20C'],stefind['Z100_REC'],stefind['Z100_NCEP'],stefind['Z100_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])
  cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread',
       '20th cent. reanal.','Brönnimann 2009, NCEP reanal.','ERA40'),
       lwd=c(2,12,2,12,2,2,2,2),lty=c(1,1,1,1,1,1,1,1), 
       col=c('black','lightgrey','red','indianred1','blue','magenta','yellow','cyan'),
       cex=fs, bty='o', bg='white', box.col='white') 

# data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z100'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
# abline(v=vf,col='black',lty=2)
# abline(h=0,col='black',lty=2)
# ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
# legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],',',ecorbem[2],',',ecorbem[3],',',ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],',',acorbem[2],',',acorbem[3],',',acorbem[4])),
  cex=fs, bty='o', bg='white', box.col='white')


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],',',ecorbem[2],',',ecorbem[3],',',ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],',',acorbem[2],',',acorbem[3],',',acorbem[4])),
  cex=fs, bty='o', bg='white', box.col='white')

dev.off()

# Why is recon corr with ECHAM so high (NH, SH, GLO mean temp.)?

# DO Jonas plots for general skill!


} # end load_prepplot











 
# PLOT atm indices
load(paste0(prepplotpath,'indices_allts.Rdata')) # some ocean based indices only NOT in landonly version
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table(paste0(datadir,"indices/stefan/stefan_seasonal_indices.txt"),header=T,
                      skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
                      seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
header <- read.table(paste0(datadir,"indices/stefan/stefan_seasonal_indices.txt"),
                   nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
                   seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
                   # (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(header))
for (i in 1:ncol(header)) {
 newhead[i] <- paste(header[1,i],header[2,i],sep='_')
}
colnames(stefind) <- newhead
# #load('../data/indices_recon.Rdata')
# vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]
#
nao_trou <- read.table(paste0(datadir,"indices/www/nao-trouet2009.txt"),skip=93,header=T)
nao_cook <- read.table(paste0(datadir,"indices/www/nao_cook2002.txt"),skip=50,header=F)
nao_lut_seas <- read.table(paste0(datadir,"indices/www/nao_sea_luterbacher_2002.txt"),skip=19,header=T)
nao_lut_mon <- read.table(paste0(datadir,"indices/www/nao_mon_luterbacher_2002.txt"),skip=27,header=T)
pna_trou <- read.table(paste0(datadir,"indices/www/pna-trouet2009.txt"),skip=75,header=F)

pdf(file='../figures/nat_data_paper/nao_pna_index.pdf',width=10,height=10) 
wlen <- 1
linew <- 1
par(mfrow=c(3,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
#eind.allts$names

# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,plotts=c(1,300),title='NAO')
#abline(v=vf,col='black',lty=2)
naols <- nao_lut_seas[nao_lut_seas[,2]=="wi",] # seas. averages until 1659
naols <- ts(naols[naols[,1]>=1600,3],start=1600)
pos <- nao_lut_mon[,2]=="Dec" | nao_lut_mon[,2]=="Jan" | 
       nao_lut_mon[,2]=="Feb" | nao_lut_mon[,2]=="Mar"
naolm <- nao_lut_mon[pos,]
pos2 <- which(naolm[,2]=="Dec")
naolm[pos2,1] <- (naolm[pos2,1]+1) # move Dec to next year to make winter season annual avg
naolm <- aggregate(naolm[,3],by=list(naolm[,1]),mean) # monthly data since 1659
naolm <- ts(naolm[,2],start=naolm[1,1],freq=1)
naol <- ts(c(naols,naolm),start=tsp(naols)[1]) # annual winter NAO
naot <- ts(nao_trou[nao_trou[,1]>=1600,2],start=1600)
naoc <- ts(nao_cook[nao_cook[,1]>=1600,2],start=1600)
vdata=NULL
vdata <- ts(cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40']),
           start=stefind[1,1],freq=1)
vdata <- ts.union(naol,naot,naoc,vdata)
yrv <- seq(1600,1600+nrow(vdata)-1)
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
vdata <- scale(vdata,center=centerfac,scale=scalefac)
if (wlen > 1) { vdata <- runmean(vdata,wlen) }
lines(vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1602 & data[['period']]<1901)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata,1603,1900),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata,1603,1900),
                      use="pairwise.complete.obs"),digits=2) 

#ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('topleft', c('EKF400','Luterbacher','Trouet','Cook'),
       lty=rep(1,5),lwd=c(2,1,1,1,1),
       col=c(rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
             rgb(0,10,10,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')
legend('bottomleft',c(paste('CCC400 ens. mean - Lut., Trou., Cook:',
       ecormean[1],ecormean[2],ecormean[3]),
       paste('EKF400 ens. mean - Lut., Trou., Cook:',acormean[1],acormean[2],acormean[3])),
       cex=fs, bty='o', bg='white', box.col='white')
#,paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4]))

data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),
                   lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,
                   plotts=c(299,402),title='NAO 20th century')
lines(vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,5],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,6],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,7],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1900 & data[['period']]<2001)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata,1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata,1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
legend('topleft', c('EKF400','20CR','REC','NCEP','ERA40'),
       lty=rep(1,5),lwd=c(2,1,1,1,1),
       col=c(rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
             rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')
legend('bottomleft', 
  c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[4],ecormean[5],ecormean[6],ecormean[7]),
  paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40:',acormean[4],acormean[5],acormean[6],acormean[7])),
  cex=fs, bty='o', bg='white', box.col='white')
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[4],ecorbem[5],ecorbem[3],ecorbem[6]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[4],acorbem[5],acorbem[6],acorbem[7])


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,title='PNA')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
vdata <- ts(vdata,start=yrv[1],freq=1)
pnat <-ts(pna_trou[,2],start=pna_trou[1,1],freq=1)
vdata <- ts.union(vdata,pnat)
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(vdata[,5],ty='l',lwd=linew,lty=1,col=rgb(0,10,0,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1900 & data[['period']]<2001)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata[,1:5],1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata[,1:5],1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
pos2 <- which(data[['period']]>1724 & data[['period']]<2001)
ecormean[5] <- round(cor(data[['allind']][pos2,1],window(vdata[,5],1725,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean[5] <- round(cor(data[['allind']][pos2,2],window(vdata[,5],1725,2000),
                      use="pairwise.complete.obs"),digits=2) 
#ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40, Trouet:',
  ecormean[1],ecormean[2],ecormean[3],ecormean[4],ecormean[5]),paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40, Trouet:',
  acormean[1],acormean[2],acormean[3],acormean[4],acormean[5])),cex=fs, bty='o', bg='white', box.col='white')
# ,paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),
# paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4]))
legend('topleft', c('EKF400','20CR','REC','NCEP','ERA40','Trouet'),
       lty=rep(1,6),lwd=c(2,1,1,1,1,1),
       col=c(rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
             rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10),
             rgb(0,10,0,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')
dev.off()








# PLOT atm indices new short version
load(paste0(prepplotpath,'indices_allts.Rdata')) # some ocean based indices only NOT in landonly version
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt"
                      ,header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),
                      seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),
                      seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
                   nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),
                   seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead

pdf(file='../figures/indices/nao_pna_index_comp3.pdf',width=10,height=10) 
wlen <- 1
linew <- 2
par(mfrow=c(2,1))
par(mai=c(0.5,0.5,0.1,0.1))
fs <- 1.0
anomper=(1901:1980)

# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),
                   lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotts=c(299,402),title=' ')
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
legend('bottomleft', c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             ecormean[1],ecormean[2],ecormean[3],ecormean[4]),
                       paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             acormean[1],acormean[2],acormean[3],acormean[4])),
                       cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('CCC400 NAO','EKF400 NAO','20CR NAO','REC NAO','NCEP NAO','ERA40 NAO'),
       lty=rep(1,6),lwd=c(2,2,1,1,1,1),
       col=c(rgb(0,0,0,10,maxColorValue=10),rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
             rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')

# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),
                   lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotts=c(299,402),title=' ')
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
legend('bottomleft', c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             ecormean[1],ecormean[2],ecormean[3],ecormean[4]),
                       paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             acormean[1],acormean[2],acormean[3],acormean[4])),
                       cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('CCC400 PNA','EKF400 PNA','20CR PNA','REC PNA','NCEP PNA','ERA40 PNA'),
              lty=rep(1,6),lwd=c(2,2,1,1,1,1),
              col=c(rgb(0,0,0,10,maxColorValue=10),rgb(10,0,0,10,maxColorValue=10),
                    rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
                    rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10)),
              cex=fs, bty='o', bg='white', box.col='white')
dev.off()



} # end timeseriesplots

























if (drought_1790) {
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1760,1830)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==1760){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c
    
  # Asia PDSI comparison years 1790, 1792-96
#  for (cyr in c(1790, 1792, 1793, 1794, 1795, 1796)) {
# define drought period
  syr <- 1790
  eyr <- 1799

  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    echam.abs <- echam
    analysis.abs <- analysis
    bemlist <- read.table(file='../data/bem/bem.txt',header=F)
    bem <- bemlist[which(bemlist[,1]==cyr),2]
    echam.abs$bem <- echam.abs$data[,,bem]
    echam.anom$bem <- echam.anom$data[,,bem]
    analysis.abs$bem <- analysis.abs$data[,,bem]
    analysis.anom$bem <- analysis.anom$data[,,bem]
    ech_ind$bem <- ech_ind$data[,,bem]
    ana_ind$bem <- ana_ind$data[,,bem]
    
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      a2tmp <- analysis.abs$ensmean
      e2tmp <- echam.abs$ensmean
      a3tmp <- analysis.anom$bem
      e3tmp <- echam.anom$bem
      vtmp <- validate$data
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      a2tmp <- a2tmp+analysis.abs$ensmean
      e2tmp <- e2tmp+echam.abs$ensmean
      a3tmp <- a3tmp+analysis.anom$bem
      e3tmp <- e3tmp+echam.anom$bem
      vtmp <- vtmp+validate$data
      c=c+1
    }
  }
  anameananom <- atmp/c
  echmeananom <- etmp/c
  eadiff <- anameananom-echmeananom
  anameanabs <- a2tmp/c
  echmeanabs <- e2tmp/c
  eadiffabs <- anameanabs-echmeanabs
  anameananombem <- a3tmp/c
  echmeananombem <- e3tmp/c
  eadiffbem <- anameananombem-echmeananombem
  valmean <- vtmp/c
  valmeananom <- valmean-valclim

  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)



  # open cook pdsi recons
  # first interpolate to echam global grid:
  # cdo -m 'NaNf' -setmissval,9e33 MADApdsi.nc madapdsi2.nc
  # cdo -r remapbil,t63grid madapdsi2.nc mada_pdsi_echamgrid.nc
  # cdo -m 'NaNf' -setmissval,9e33 NADAv2-2008.nc nadapdsi2.nc
  # cdo -r remapbil,t63grid nadapdsi2.nc nada_pdsi_echamgrid.nc
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                      small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                  small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  lut <- valmeananom[(which(validate$names=='precip')),2]
  lut[is.na(lut)] <- 0
#  nada$data <- scale(nada$data,scale=sd(lut))
#  mada$data <- scale(mada$data,scale=sd(lut))
  pdsi <- nada
  # scale factor to make precip and pdsi comparable
  f=2
  tmp <- nada$data*f + mada$data*f + lut
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  pdsi$data <- rep(pdsi$data,11)
  pdsi$names <- rep('precip',length(pdsi$names))


  
# Figure 2
  pdata <- echam
  t=2 # summer season as drought recon is JJA
  pdf(paste0('../figures/figure_2_',syr,'-',eyr,'_v2.pdf'), width=6, height=7, paper='special')
#  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
#  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
  pdata$data <- array(c(anameananom[,t],eadiff[,t]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-10,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)'),
             colnames=c("EKF400 anomaly","EKF400 update"),
             statpanel=2, add=T, rownames='temperature and GHP500', main='1790-1799',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  pdata$data <- array(c(anameananom[,t],pdsi$data), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-10,90),
             names=c('c)', 'd)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation/drought',
             colnames=c("EKF400 anomaly","reconstructions"))
  
  dev.off()
#   
#   # read mann et al temp recon and plot 1816
#   mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
#                      na.string="NaN")
#   ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
#                    na.string="NaN")
#   ts <- which(mann[,1]==1816)
#   plotmann <- plotdata
#   mdata <- as.vector(t(mann[ts,-1]))
#   #plotmann$data <- array(rep(mdata,2),c(length(mdata),1,2))
#   plotmann$data <- array(mdata,c(length(mdata),1))
#   #plotmann$data <- array(mdata,,c(length(mdata),1,1))
#   #plotmann$ensmean <- array(mdata,c(length(mdata),1,1))
#   #plotmann$time <- rep(mann[ts,1],2)
#   plotmann$time <- mann[ts,1]
#   plotmann$names <- rep('temp2',dim(plotmann$data)[1])
#   plotmann$lon <- ll[,1] 
#   plotmann$lat <- ll[,2]
#   plotmann$height <- ll[,3]   
#   plotmann$lsm.i <- ll[,3]  
#   layout(matrix(c(1,2,3,3), 2, 1, byrow = TRUE), height=c(3,1))
#   #layout(matrix(c(1,2,2), 2, 2, byrow = TRUE), height=c(3,1))
#   par(oma=c(0,0,0,0))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
#               lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')

  
  
  
  
  
  
  
  
  
  # make plots
  t=2 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,t],anameananom[,t],eadiff[,t],valmeananom[,t]), c(nrow(echmeananom),1,4))
  plotdata$data <- array(c(anameananom[,t],valmeananom[,t]), c(nrow(echmeananom),1,2))  
  pdf(paste0('../figures/1790-1800/drought_',syr,'-',eyr,'.pdf'), width=9, height=9, paper='special')
  layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]],
             lev=levs, st.col=NULL, stations=NULL, add=T) #, addcontours=T, contvarname='gph500')
  
#  add_contour(plotdata, varname = "temp2", ti=1)
  plotdata$data <- array(c(eadiff[,t],pdsi$data), c(nrow(echmeananom),1,2))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
             lev=levs, st.col=NULL, stations=NULL, add=T)
  dev.off()


# search explanation for droughts: do HC index, omega500 and gph500 anomaly indicate subtropical subsidence?
  t=2 # summer season as drought recon is JJA
  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/drought/gph500_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T) #,
#             addcontours=T, contvarname='gph500',conttype='data',contcol='white')
  dev.off()

  plotdata$data <- array(c(echmeananombem[,t],anameananombem[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/drought/gph500_anom_bem_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-12,12,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/drought/temp_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.1), Inf)
  plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T) #, addcontours=T, contvarname='gph500')
  dev.off()

  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/drought/slp_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.2), Inf)
  plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

  plotdata$data <- array(c(eadiff[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/drought/omega500_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-0.001,0.001,0.0001), Inf)
  plot_echam(plotdata, varname='omega500', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

} # end 1790 droughts





#if (cold_1640) {}
#if (cold_1690) {} 

if (warmcold_decades) {
  ptm1 <- proc.time()
#  syr <- 1640
#  eyr <- 1649
  syr <- 1641
  eyr <- 1645
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1640 <- atmp/c
#  echmeananom1640 <- etmp/c
  anameananom1641 <- atmp/c
  echmeananom1641 <- etmp/c

  syr <- 1658
  eyr <- 1662
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0('../data/prepplot_season/analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
  anameananom1658 <- atmp/c
  echmeananom1658 <- etmp/c
    
#  syr <- 1660
#  eyr <- 1669
  syr <- 1695
  eyr <- 1699
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0('../data/prepplot_season/analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1660 <- atmp/c
#  echmeananom1660 <- etmp/c
  anameananom1695 <- atmp/c
  echmeananom1695 <- etmp/c
  
#  syr <- 1690
#  eyr <- 1699
  syr <- 1676
  eyr <- 1680
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1690 <- atmp/c
#  echmeananom1690 <- etmp/c
  anameananom1676 <- atmp/c
  echmeananom1676 <- etmp/c
  
#  syr <- 1790
#  eyr <- 1799
  syr <- 1791
  eyr <- 1795
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1790 <- atmp/c
#  echmeananom1790 <- etmp/c
  anameananom1791 <- atmp/c
  echmeananom1791 <- etmp/c
  
#  syr <- 1800
#  eyr <- 1809
  syr <- 1801
  eyr <- 1805
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1800 <- atmp/c
#  echmeananom1800 <- etmp/c
  anameananom1801 <- atmp/c
  echmeananom1801 <- etmp/c
  
#  syr <- 1810
#  eyr <- 1819
  syr <- 1815
  eyr <- 1819
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1810 <- atmp/c
#  echmeananom1810 <- etmp/c
  anameananom1815 <- atmp/c
  echmeananom1815 <- etmp/c
  
#  syr <- 1820
#  eyr <- 1829
  syr <- 1826
  eyr <- 1830
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1820 <- atmp/c
#  echmeananom1820 <- etmp/c
  anameananom1826 <- atmp/c
  echmeananom1826 <- etmp/c
  
#  syr <- 1830
#  eyr <- 1839
  syr <- 1832
  eyr <- 1836
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1830 <- atmp/c
#  echmeananom1830 <- etmp/c
  anameananom1832 <- atmp/c
  echmeananom1832 <- etmp/c
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
# # convert monthly_out from 6-mon state vector because of now fixed bug in EnSRF_prepplot
#  echmeananom1791 <- array(echmeananom1791,c((dim(echmeananom1791)[1]/6),
#                       dim(echmeananom1791)[2]*6))
#  anameananom1791 <- array(anameananom1791,c((dim(anameananom1791)[1]/6),
#                                             dim(anameananom1791)[2]*6))


  # Figure 2b 5-yr cold period sommer
  pdata <- echam
  t=2 # 
#  pdf('../figures/figure_2b_v3.pdf', width=12, height=7, paper='special')
  pdf('../figures/figure_2b_som_v4.pdf', width=12, height=7, paper='special')
  #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
  layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
  #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
#  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t],anameananom1810[,t],anameananom1830[,t]), c(nrow(echam.anom$ensmean),1,4))
  pdata$data <- array(c(anameananom1641[,t],anameananom1695[,t],anameananom1815[,t],anameananom1832[,t]), c(nrow(echam.anom$ensmean),1,4))
  
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1641-45","EKF400 1695-99","EKF400 1815-19","EKF400 1832-36"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Cold 5-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
#  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-5,5,1), Inf)
#  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()
  
  # Figure 2b  5-yr cold period winter
  pdata <- echam
  t=1 # see volc high lat winter warming
  #  pdf('../figures/figure_2b_v3.pdf', width=12, height=7, paper='special')
  pdf('../figures/figure_2b_win_v4.pdf', width=12, height=7, paper='special')
  #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
  layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
  #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
  #  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t],anameananom1810[,t],anameananom1830[,t]), c(nrow(echam.anom$ensmean),1,4))
  pdata$data <- array(c(anameananom1641[,t],anameananom1695[,t],anameananom1815[,t],anameananom1832[,t]), c(nrow(echam.anom$ensmean),1,4))
  
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1641-45","EKF400 1695-99","EKF400 1815-19","EKF400 1832-36"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Cold 5-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
  #  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-5,5,1), Inf)
  #  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()
  
  
  # Figure 2c  5-yr warm period summer
  pdata <- echam
  t=2 # summer season as drought recon is JJA
#  pdf('../figures/figure_2c_som_v3.pdf', width=12, height=7, paper='special')
  pdf('../figures/figure_2c_som_v4.pdf', width=12, height=7, paper='special')
  #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
  layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
  #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
#  pdata$data <- array(c(anameananom1660[,t],anameananom1790[,t],anameananom1800[,t],anameananom1820[,t]), c(nrow(echam.anom$ensmean),1,4))
  pdata$data <- array(c(anameananom1658[,t],anameananom1791[,t],anameananom1801[,t],anameananom1826[,t]), c(nrow(echam.anom$ensmean),1,4))
  
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1658-62","EKF400 1791-95","EKF400 1801-05","EKF400 1826-30"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
  #  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-5,5,1), Inf)
  #  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()
  
  # Figure 2 5-yr warm period winter
  pdata <- echam
  t=1 
  #  pdf('../figures/figure_2c_v3.pdf', width=12, height=7, paper='special')
  pdf('../figures/figure_2c_win_v4.pdf', width=12, height=7, paper='special')
  #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
  layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
  #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
  #  pdata$data <- array(c(anameananom1660[,t],anameananom1790[,t],anameananom1800[,t],anameananom1820[,t]), c(nrow(echam.anom$ensmean),1,4))
  pdata$data <- array(c(anameananom1658[,t],anameananom1791[,t],anameananom1801[,t],anameananom1826[,t]), c(nrow(echam.anom$ensmean),1,4))
  
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1658-62","EKF400 1791-95","EKF400 1801-05","EKF400 1826-30"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
  #  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-5,5,1), Inf)
  #  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()

#     
#   # Figure 2e SON, DJF, MAM, JJA to see if better Z500 signal in winter
#   pdata <- echam
#   t=2 # look at winters of 5-yr warm periods
#   #  pdf('../figures/figure_2c_v3.pdf', width=12, height=7, paper='special')
#   pdf('../figures/figure_2c_v4.pdf', width=12, height=7, paper='special')
#   #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
#   # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
#   layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
#   #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#   par(oma=c(0.5,4,4,0))
#   levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
#   contlevs <- seq(-4,4,2)
#   #  pdata$data <- array(c(anameananom1660[,t],anameananom1790[,t],anameananom1800[,t],anameananom1820[,t]), c(nrow(echam.anom$ensmean),1,4))
#   pdata$data <- array(c(anameananom1658[,t],anameananom1791[,t],anameananom1801[,t],anameananom1826[,t]), c(nrow(echam.anom$ensmean),1,4))
#   
#   plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
#              lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
#              colnames=c("EKF400 1658-62","EKF400 1791-95","EKF400 1801-05","EKF400 1826-30"),
#              statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
#              addcontours=T, contvarname='gph500',conttype='data',contcol='black',
#              contlev=contlevs)
#   #  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
#   levs <- c(-Inf, seq(-5,5,1), Inf)
#   #  pdata$data <- array(c(anameananom1640[,t],anameananom1690[,t]), c(nrow(echam.anom$ensmean),1,2))
#   plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
#              names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
#              addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
#   
#   dev.off()  
#}
} # end warmcold_decades









if (plots1816) {
  # check yuri's 1816 data compilation
  #load("../comparison_data/yuri_1816_all_station_and_travel_pressure_data.Rdata")
  syr=1816     
  eyr=1817
  
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1781,1851)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==1781){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c

  # European T, P and SLP comparison year 1816
  cyr = 1816
  load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
  validate.anom <- echam
  validate.anom$data <- echam$ensmean
  validate.anom$data[(1:nrow(validate$data)),] <- validate$data-valclim
  validate.anom$data[(nrow(validate$data)+1):nrow(echam$data),] <- NA
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
  # load station location where data was assimilated
  pcoor <- read.table("../data/coor/prox_coor_1816.5.csv")
  scoor <- read.table("../data/coor/stat_coor_1816.5.csv")
  station <- list(lon=c(scoor[,1],pcoor[,1]),lat=c(scoor[,2],pcoor[,2]),
                  data=matrix(c(rep(1,nrow(pcoor)),rep(2,nrow(scoor))),ncol=1))
  
  # make plots
  tim=2 # summer season as we look at year without a summer
  pdata=echam
  pdata$lon <- pdata$lon[1:dim(validate.anom)[1]]
  pdata$lat <- pdata$lat[1:dim(validate.anom)[1]]
  pdata$names <- pdata$names[1:dim(validate.anom)[1]]
  pdf('../figures/nat_data_paper/eu_1816.pdf', width=9, height=9, paper='special')
    layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,1,3,1))
    #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(2,4,4,0))
    levs <- c(-Inf, seq(-1.75,1.75,0.5), Inf)
    contlevs <- seq(-2,2,0.5)
    pdata$data <- array(c(echam.anom$ensmean[,tim],analysis.anom$ensmean[,tim],
                      validate.anom$data[,tim]), 
                      c(nrow(validate.anom$data),1,3))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=5.2, lonlim=c(-20,40),
             latlim=c(30,75), wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)','c)'),
             colnames=c("CCC400 anomaly","EKF400 anomaly","Reconstruction anomaly"),
             statpanel=2, add=T, rownames='Temperature and SLP', main='1816', #no ghp in vali data
             addcontours=T, contvarname='slp',conttype='data',contcol='black',
             contlev=contlevs)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    levs <- c(-Inf, seq(-9,9,2), Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=5.2, lonlim=c(-20,40),
             latlim=c(30,75), wcol='darkgrey',
             lev=levs, st.col=NULL, stations=NULL, names=c('d)', 'e)','f)'),
             colnames='',
             statpanel=NULL, add=T, rownames='Precipitation and 850hPa wind', main='',
             addvectors=T, vecnames=c('u850','v850'), veccol='black', 
             veclen=scfac*0.01, vecscale=scfac*0.9, vecwd=0.95, every_x_vec=1,
             colorbar=T)
  dev.off()
  
  
  
  #  plotdata$data <- array(c(echmeananom[,t],anameananom[,t],eadiff[,t],valmeananom[,t]), c(nrow(echmeananom),1,4))
  #plotdata$data <- array(c(anameananom[,t],valmeananom[,t]), c(nrow(echmeananom),1,2))
#  plotdata$data <- array(c(anameananom[,t],anameananombem[,t],valmeananom[,t]), c(nrow(echmeananom),1,3))  
# 
#   pdf('../figures/1816/1816.pdf', width=13.5, height=9, paper='special')
# #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
#   layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8,9,10,11,12,12,12), 6, 3, byrow = TRUE), height=c(3,1,3,1,3,1))
#   par(oma=c(0,0,0,0))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:2],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   levs <- c(-Inf, seq(-8,8,2), Inf)
#   plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   levs <- c(-Inf, seq(-4,4,1), Inf)
#   plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, names=pnames[5:6],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   dev.off()
  
  
  
  
  
  
# Figure 3: 1816 alternative  
  t=2 # for summer 1816
  pdata<-echam
#   pdata$data <- array(c(echmeananom[,t],anameananom[,t],anameananombem[,t]), c(nrow(echmeananom),1,3))  
#   pdf('../figures/figure_3_1816_alltemp.pdf', width=10, height=6, paper='special')
#   layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,0.5,3,1))
#   par(oma=c(2,2,4,2))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names=c('a)', 'b)', 'c)'),
#               lev=levs, st.col=NULL, stations=calibrate, statpanel=2,
#               add=T ,units='ºC',seas=1,colnames=c("CCC","EKF","BEM"),
#               rownames=c("this study","reconstructions"),
#               latlim=c(0,90),main='1816',colorbar=F)
  pdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echmeananom),1,2))  
  pdf('../figures/figure_3_1816_alltemp_v3.pdf', width=10, height=6, paper='special')
  layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,0.5,3,1))
  par(oma=c(2,2,4,2))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names=c('a)', 'b)'),
                lev=levs, st.col=NULL, stations=calibrate, statpanel=2,
                add=T ,units='ºC',seas=1,colnames=c("CCC","EKF"),
                rownames="this study",latlim=c(0,90),main='1816',colorbar=F)  
  pdata$data <- array(anameananombem[,t], c(nrow(echmeananom),1,1))  
  coor <- read.csv('/Users/joerg/Documents/unibe/projects/2013_small/pages_paper/ccc400_coor_all.csv')
  stat35 <- calibrate
  stat35$lon <- coor$"echam_lon"
  stat35$lon[stat35$lon>180] <- stat35$lon[stat35$lon>180]-360
  stat35$lat <- coor$"echam_lat"
  stat35$names <- rep('prox',length(stat35$lon))
  plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names='c)',
               lev=levs, st.col=NULL, stations=stat35, statpanel=1,
               add=T ,units='ºC',seas=1,colnames="BEM",
               latlim=c(0,90),main='',colorbar=F) #,centercol='lightgrey')    
    
  plot(NA,axes=F,bty='n',ylim=c(1,2))
  
# Luterbacher recon  
  pdata$data <- array(valmeananom[,t], c(nrow(valmeananom),1,1)) 
  pdata$lon <- validate$lon  
  pdata$lat <- validate$lat
  pdata$names <- validate$names
#  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5,latlim=c(0,90),
              seas=1,colorbar=F,rownames='reconstructions',colnames='Luterbacher',
              lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',
              names='d)') #,centercol='lightgrey')
  
# read mann et al temp recon and plot 1816
  mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                     na.string="NaN")
  ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                   na.string="NaN")
  ts <- which(mann[,1]==1816)
  plotmann <- pdata
  mdata <- as.vector(t(mann[ts,-1]))
  mdata2 <- apply(mann[(ts-35):(ts+35),-1],2,mean)
  plotmann$data <- array(mdata-mdata2,c(length(mdata),1))
  plotmann$time <- mann[ts,1]
  plotmann$names <- rep('temp2',dim(plotmann$data)[1])
  plotmann$lon <- ll[,1] 
  plotmann$lat <- ll[,2]
  plotmann$height <- ll[,3]   
  plotmann$lsm.i <- ll[,3]  
#  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(plotmann, varname='temp2', type='data', cex.pt=2.,names='e)',
              latlim=c(0,90),seas=1,colorbar=F,wcol="black",colnames='Mann',
              lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC') #,centercol='lightgrey')
  
  schw <- read.table("/Users/joerg/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_mxdabd_grid.dat",
                     na.string="-9.990")
  lon <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lon.dat")
  lat <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lat.dat")
  yrs <- seq(1400,1994)
  ts <- which(yrs==1816)
  plotschw <- echam
  sdata <- as.vector(t(schw[ts,]))
  sdata2 <- apply(schw[(ts-35):(ts+35),],2,mean)
  plotschw$data <- array(scale(sdata-sdata2,center=F),c(length(sdata),1))
  plotschw$ensmean <- array(sdata,c(length(sdata),1))
  plotschw$time <- schw[ts,1]
  plotschw$names <- rep('temp2',dim(plotschw$data)[1])
  plotschw$lon <- as.vector(t(lon))
  plotschw$lat <- as.vector(t(lat))
  plotschw$height <- rep(NA,length(lon))
  plotschw$lsm.i <- rep(NA,length(lon))
#  levs <- c(-Inf, seq(-3,3,0.5), Inf)
  plot_echam(plotschw, varname='temp2', type='data', cex.pt=2.,names='f)',
             latlim=c(0,90),seas=1,colorbar=T,wcol="black",colnames='Briffa',
             lev=levs, st.col=NULL, stations=NULL, add=T,units='') #,centercol='lightgrey')
  
dev.off()
  
 
  
  
  
  


  
  
plotdata$data <- array(echmeananom[,t], c(nrow(echmeananom),1))  
pdf('../figures/1816/1816_CCC400.pdf', width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='CCC400')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(anameananom[,t], c(nrow(echmeananom),1))  
pdf('../figures/1816/1816_EKF.pdf', width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='EKF')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(eadiff[,t], c(nrow(echmeananom),1))  
pdf('../figures/1816/1816_EKF-CCC400_all.pdf', width=8, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-20,20,5), Inf)
plot_echam2(plotdata, varname='gph500', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='hPa')
levs <- c(-Inf, seq(-1,1,0.25), Inf)
plot_echam2(plotdata, varname='u850', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='m/s')
plot_echam2(plotdata, varname='v850', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='m/s')
levs <- c(-Inf, seq(-0.005,0.005,0.001), Inf)
plot_echam2(plotdata, varname='omega500', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=NULL, add=T ,units='m/s')
dev.off()

plotdata$data <- array(anameananombem[,t], c(nrow(echmeananom),1))  
pdf('../figures/1816/1816_BEM.pdf', width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='BEM')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(valmeananom[,t], c(nrow(valmeananom),1))  
pdf('../figures/1816/1816_LUT.pdf', width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='Reconstructions')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()


  # read mann et al temp recon and plot 1816
  mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                   na.string="NaN")
  ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                 na.string="NaN")
  ts <- which(mann[,1]==1816)
  plotmann <- plotdata
  mdata <- as.vector(t(mann[ts,-1]))
  #plotmann$data <- array(rep(mdata,2),c(length(mdata),1,2))
  plotmann$data <- array(mdata,c(length(mdata),1))
  #plotmann$data <- array(mdata,,c(length(mdata),1,1))
  #plotmann$ensmean <- array(mdata,c(length(mdata),1,1))
  #plotmann$time <- rep(mann[ts,1],2)
  plotmann$time <- mann[ts,1]
  plotmann$names <- rep('temp2',dim(plotmann$data)[1])
  plotmann$lon <- ll[,1] 
  plotmann$lat <- ll[,2]
  plotmann$height <- ll[,3]   
  plotmann$lsm.i <- ll[,3]  
  pdf('../figures/1816/1816mann.pdf', width=8, height=6, paper='special')
    layout(matrix(c(1,2,3,3), 2, 1, byrow = TRUE), height=c(3,1))
  #layout(matrix(c(1,2,2), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    plot_echam2(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
             lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')
  dev.off()
  
  
# schweingruber mxd grid by rutherford and briffa
#  The (text or ASCII format file) has 115 columns and 595 rows. 
#  Each column represents one grid-box time series, 
#  while each row represents one year, starting in 1400 and ending in 1994
  schw <- read.table("/Users/joerg/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_mxdabd_grid.dat",na.string="-9.990")
  lon <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lon.dat")
  lat <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lat.dat")
  yrs <- seq(1400,1994)
  ts <- which(yrs==1816)
  plotschw <- plotdata
  sdata <- as.vector(t(schw[ts,]))
# plotschw$data <- array(rep(sdata,2),c(length(sdata),1,2))
# plotschw$ensmean <- array(rep(sdata,2),c(length(sdata),1,2))
# plotschw$time <- rep(mann[ts,1],2)
  plotschw$data <- array(sdata,c(length(sdata),1))
  plotschw$ensmean <- array(sdata,c(length(sdata),1))
  plotschw$time <- schw[ts,1]
  plotschw$names <- rep('temp2',dim(plotschw$data)[1])
  plotschw$lon <- as.vector(t(lon))
  plotschw$lat <- as.vector(t(lat))
  plotschw$height <- rep(NA,length(lon))
  plotschw$lsm.i <- rep(NA,length(lon))
  pdf('../figures/1816/1816schw.pdf', width=8, height=4, paper='special')
    layout(matrix(c(1,2), 2, 1, byrow = TRUE), height=c(2,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-3,3,0.5), Inf)
    plot_echam2(plotschw, varname='temp2', type='data', cex.pt=3.0, names='',
             lev=levs, st.col=NULL, stations=NULL, add=T, units='TRW')
  dev.off()  


# pdsi 1816
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                  small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                  small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  pdsi <- nada
  tmp <- nada$data + mada$data
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  #pdsi$data <- rep(pdsi$data,2)
  pdsi$names <- rep('precip',length(pdsi$names))

  # make plots
  t=1 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,t],anameananom[,t],eadiff[,t],valmeananom[,t]), c(nrow(echmeananom),1,4))
#  plotdata$data <- array(c(anameananom[,t],valmeananom[,t]), c(nrow(echmeananom),1,2))

  pdf('../figures/1816/1816_pdsi.pdf', width=8, height=6, paper='special')
    layout(matrix(c(1,2), 2, 1, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    plotdata$data <- array(pdsi$data, c(nrow(echmeananom),1))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T, units='PDSI')
  dev.off()



# 20CR 1816 plot, anomalies to 20cr 2c 1851-80, because no good data for earlier period
#cdo merge t_AMJJAS_1816_monmean_ano.nc p_AMJJAS_1816_monmean_ano.nc slp_AMJJAS_1816_monmean_ano.nc tpslp_AMJJAS_1816_monmean_ano.nc
scout <- read_20cr('tpslp_AMJJAS_1816_summermean_ano.nc', path=twentycrpath, 
            xlim=c(-180,180), ylim=c(-90,90), timlim=c(syr, eyr), small=F, landonly=F)

# make plots
t=1 # summer season as drought recon is JJA
plotdata=scout
#  plotdata$data <- array(c(echmeananom[,t],anameananom[,t],eadiff[,t],valmeananom[,t]), c(nrow(echmeananom),1,4))
#plotdata$data <- array(c(scout$data,scout$data), c(length(scout$data),1,2))
#plotdata$data <- array(scout$data, c(length(scout$data),1,1))
plotdata$data <- array(scout$data, c(length(scout$data),1))

pdf('../figures/1816/20cr_1816.pdf', width=4, height=10, paper='special')
#  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
  par(oma=c(0,0,4,0))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, seas=1,
        names='', lev=levs, st.col=NULL, stations=NULL, add=T, 
        main='20CR', units='ºC')
  levs <- c(-Inf, seq(-4,4,1), Inf)
  plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, 
        names='', lev=levs, st.col=NULL, stations=NULL, add=T,
        units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, 
            names='', lev=levs, st.col=NULL, stations=NULL, add=T,
            units='mm')
dev.off()


# 
# # all 1816 plots together
# plotdata$data <- array(c(anameananom[,2],anameananombem[,2],valmeananom[,2]),scout$data, c(nrow(echmeananom),1,4))
# plotdata$lon <- c(rep(echam$lon,3),scout$lon)
# plotdata$lat <- c(rep(echam$lat,3),scout$lat)
# plotdata$names <- c(rep(echam$names,3),scout$names)
# 
# pdf('../figures/1816/1816.pdf', width=13.5, height=9, paper='special')
# #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
# layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8,9,10,11,12,12,12), 6, 3, byrow = TRUE), height=c(3,1,3,1,3,1))
# par(oma=c(0,0,0,0))
# levs <- c(-Inf, seq(-2,2,0.5), Inf)
# plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:2],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# levs <- c(-Inf, seq(-8,8,2), Inf)
# plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# levs <- c(-Inf, seq(-4,4,1), Inf)
# plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, names=pnames[5:6],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# dev.off()

  
  
} # end 1816 plots 




















# find dryest, warmest, ... years in NH and specific regions
if (multiann) {
  # SW USA droughts
  # no index defined yet for that region. First check somewhere else if it works!
  # Not ideal that simulation finish 2005 and we do not have this drought in your analysis
  # northern hemisphere extra-tropics (>20N)
  par(mfrow=c(1,1))
  par(mai=c(1,1,0.5,0.5))
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='MED.precip'),seas='sum',
                     anomper=(1901:1980),lw=2,sca=F,xa='s',plotbem=F,plotmem=F,
                     plotech=T,pag=F,use18=T,yl="Precip.")
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='MED.slp'),seas='sum',
                     anomper=(1901:1980),lw=2,sca=F,xa='s',plotbem=F,plotmem=T,
                     plotech=T,pag=F,use18=T,yl="SLP")
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='ENH.temp2'),seas='sum',
                     anomper=(1901:1980),lw=2,sca=F,xa='s',plotbem=F,plotmem=T,
                     plotech=T,pag=F,use18=T,yl="Temp.")
  plot(runmean(data$allindann[,"aindmean"],31),ty='l')
  datdetr<-data$allindann[,"aindmean"]-runmean(data$allindann[,"aindmean"],31)
  plot(datdetr,ty='l')
  data$period[which(datdetr < quantile(datdetr, probs = 0.1))]
  
  # 1830-1850 glacier advances
  # NAM, NEU and MED winter precip and summer temp
  wlen=1
  pdf(file=paste0('../figures/1830-50/1830-50_winprecip_sumtemp_',wlen,'.pdf'),width=10,height=6)
  par(mfrow=c(2,3))
  par(mai=c(0.3,0.6,0.2,0))
  pdata <- aind.allts$ensmean[which(aind.allts$names=='NEU.precip'),
                          !is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NEU Winter precip.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='NEU.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2)
  pdata <- aind.allts$ensmean[which(aind.allts$names=='MED.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='MED Winter precip.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='MED.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2)
  pdata <- aind.allts$ensmean[which(aind.allts$names=='NAM.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NAM Winter precip.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='NAM.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2)
  
  pdata <- aind.allts$ensmean[which(aind.allts$names=='NEU.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NEU Summer Temp.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='NEU.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2)
  pdata <- aind.allts$ensmean[which(aind.allts$names=='MED.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='MED Summer Temp.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='MED.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2)
  pdata <- aind.allts$ensmean[which(aind.allts$names=='NAM.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NAM Summer Temp.',ylab='')
  pdata <- eind.allts$ensmean[which(aind.allts$names=='NAM.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
  lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
  abline(v=c(1830,1850),lty=2) 
  legend('topleft', c('EKF400','CCC400'),lwd=c(2,2), col=c('red','blue'), bty='o', 
         bg='white', box.col='white')  
  dev.off()
  

  
  
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1800,1880)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    if (cyr==1800){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c
  
  # 1830-50 glacier advances
  #  for (cyr in c(1790, 1792, 1793, 1794, 1795, 1796)) {
  # define period
  syr <- 1830
  eyr <- 1850
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    echam.abs <- echam
    analysis.abs <- analysis
    bemlist <- read.table(file='../data/bem/bem.txt',header=F)
    bem <- bemlist[which(bemlist[,1]==cyr),2]
    echam.abs$bem <- echam.abs$data[,,bem]
    echam.anom$bem <- echam.anom$data[,,bem]
    analysis.abs$bem <- analysis.abs$data[,,bem]
    analysis.anom$bem <- analysis.anom$data[,,bem]
    ech_ind$bem <- ech_ind$data[,,bem]
    ana_ind$bem <- ana_ind$data[,,bem]
    
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      a2tmp <- analysis.abs$ensmean
      e2tmp <- echam.abs$ensmean
      a3tmp <- analysis.anom$bem
      e3tmp <- echam.anom$bem
      vtmp <- validate$data
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      a2tmp <- a2tmp+analysis.abs$ensmean
      e2tmp <- e2tmp+echam.abs$ensmean
      a3tmp <- a3tmp+analysis.anom$bem
      e3tmp <- e3tmp+echam.anom$bem
      vtmp <- vtmp+validate$data
      c=c+1
    }
  }
  anameananom <- atmp/c
  echmeananom <- etmp/c
  eadiff <- anameananom-echmeananom
  anameanabs <- a2tmp/c
  echmeanabs <- e2tmp/c
  eadiffabs <- anameanabs-echmeanabs
  anameananombem <- a3tmp/c
  echmeananombem <- e3tmp/c
  eadiffbem <- anameananombem-echmeananombem
  valmean <- vtmp/c
  valmeananom <- valmean-valclim
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
  
  
  # open cook pdsi recons
  # first interpolate to echam global grid:
  # cdo -m 'NaNf' -setmissval,9e33 MADApdsi.nc madapdsi2.nc
  # cdo -r remapbil,t63grid madapdsi2.nc mada_pdsi_echamgrid.nc
  # cdo -m 'NaNf' -setmissval,9e33 NADAv2-2008.nc nadapdsi2.nc
  # cdo -r remapbil,t63grid nadapdsi2.nc nada_pdsi_echamgrid.nc
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                    path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                    small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                    path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                    small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  pdsi <- nada
  tmp <- nada$data + mada$data
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  pdsi$data <- rep(pdsi$data,2)
  pdsi$names <- rep('precip',length(pdsi$names))
  
  # make plots
  t=2 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,t],anameananom[,t],eadiff[,t],valmeananom[,t]), c(nrow(echmeananom),1,4))
  plotdata$data <- array(c(anameananom[,t],valmeananom[,t]), c(nrow(echmeananom),1,2))
  
  pdf(paste0('../figures/1830-50/drought_',syr,'-',eyr,'.pdf'), width=9, height=9, paper='special')
  layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]],
             lev=levs, st.col=NULL, stations=NULL, add=T)
  
  plotdata$data <- array(c(eadiff[,t],pdsi$data), c(nrow(echmeananom),1,2))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
             lev=levs, st.col=NULL, stations=NULL, add=T)
  dev.off()
  
  
  # search explanation for droughts: do HC index, omega500 and gph500 anomaly indicate subtropical subsidence?
  t=2 # summer season as drought recon is JJA
  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/1830-50/gph500_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananombem[,t],anameananombem[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/1830-50/gph500_anom_bem_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-12,12,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/1830-50/temp_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.1), Inf)
  plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/1830-50/slp_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.2), Inf)
  plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(eadiff[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/1830-50/omega500_anom_',syr,'-',eyr,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-0.001,0.001,0.0001), Inf)
  plot_echam(plotdata, varname='omega500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  
  
  
  
  # Figure 3 1830-1840 eruptions/cold decade (or 1816 for anniversary version?)
  pdata <- echam
#  t=2 # summer season, but we want also winter precip
  pdf(paste0('../figures/figure_3_',syr,'-',eyr,'.pdf'), width=6, height=6, paper='special')
  #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
  #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
  pdata$data <- array(anameananom[,1:2], c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-10,90),wcol='darkgrey',
             names=pnames[1:dim(pdata$data)[3]], lev=levs, st.col=NULL, stations=calibrate, 
             statpanel=c(1,2), add=T,main="Temp. & GPH500 anom.",
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',contlev=contlevs)
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  pdata$data <- array(anameananom[,1:2], c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-10,90),
             names=pnames[3:4], lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey')
  
  dev.off()
  

} # end 1830-50














if (warmcold_us) {
  #  # find and mark coldest/warmest 5yr period in US
  #   d <- aind.allts$ensmean[(which(eind.allts$names=="NAM.temp2")),
  #                           (is.even(seq(1,ncol(aind.allts$ensmean))))]
  #   dec <- rep(NA,250)
  #   dectime <- rep(NA,250)
  #   for (t in 1:length(dec)) {
  #     if (t==1) {
  #       dec[t] <- mean(d[1:5])
  #       dectime[t] <- 1603
  #       i=2
  #     } else {
  #       dec[t] <- mean(d[i:(i+4)])
  #       dectime[t] <- dectime[t-1] + 1 
  #       i=i+1
  #     }
  #   } 
  #   #  plot(dectime,dec,ty='l',col='red')
  #   tmp=cbind(dec,dectime)
  #   tmp[order(dec),]  
  coldper <- c(1641, 1697, 1835, 1667, 1738, 1836)
  #  warmper <- c(1825, 1790, 1842, 1800, 1772, 1733)
  for (syr in coldper) {
    eyr <- syr+4
    # load nada pdsi recon
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)
    #
    #  plotdata <- echam
    #  plotdata$ensmean <- echmeananom # 1641
    #  scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    #  plot_echam3(plotdata, varname='precip', type='ensmean', cex.pt=6, ti=2,
    #              levs=c(50,90,95,100,105,110,200),
    #              latlim=c(20,60),lonlim=c(-140,-40),
    #              main='Precip. & 200hPa wind anomalies', units='%',
    #              addvectors=T, vecnames=c('u200','v200'), veccol='black', 
    #              veclen=scfac*0.01, vecscale=scfac*0.9, vecwd=0.95, every_x_vec=1,
    #              wcol='darkgrey')
    #  
    #   # change to precip in % 
    #   load(file=paste0('../data/prepplot_season/analysis_',cyr,'.Rdata'))
    #   echam.abs <- echam
    #   load(file=paste0('../data/echam_clim/echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))
    #   # seasonal mean
    #   echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),rowMeans(echam_clim$ensmean[,16:21]))
    #   echam$ensmean <- echmeanclim
    #   plot_echam(echam, varname='precip', type='ensmean', cex.pt=2, ti=1)
    #   pos <- which(echam$names=="precip")[1:4608]
    #   echam$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
    ##   quartz()
    #   plot_echam(echam,levs=c(50,85,90,95,100,105,110,115,200), varname='precip', 
    #              type='ensmean', cex.pt=2, ti=1)
    
    # Figure ??  5-yr cold period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    #   t=2 # summer
    pdf(paste0('../figures/figure_cold_us_',s,'_v1.pdf'), width=12, 
          height=14, paper='special')
#      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
      layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
             16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
      par(oma=c(0.5,4,4,0))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      scfac <- max(echam.anom$ensmean[echam$names=='u200'])
#      pdata$data <- array(c(echmeananom[,t],anameananom[,t]), c(nrow(echam.anom$ensmean),1,2))
      pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                            anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
#      pdata$data <- array(c(echmeananom[,t],echmeananom[,t+2],echmeananom[,t+4],
#                            echmeananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))      
#      pdata$dataechmeanclim[pos,]/echam$ensmean[pos,]*100
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                  statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                  main='EKF400 5-yr cold periods in North America', units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
      # diff echam analysis
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                            (anameananom[,t+2]-echmeananom[,t+2]),
                            (anameananom[,t+4]-echmeananom[,t+4]),
                            (anameananom[,t+6]-echmeananom[,t+6])),
                            c(nrow(echam.anom$ensmean),1,4))
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                  statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                  units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
      pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                            anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
#      #    levs <- c(-Inf, seq(-5,5,1), Inf)
      levs=c(50,80,85,90,95,100,105,110,115,120,200)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),units="%",
                  names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
      pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                            pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
      levs <- seq(-4,4,1)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
#  } # end cold periods loop 
  
  warmper <- c(1825, 1790, 1842, 1800, 1772, 1733)
  for (syr in warmper) {
    eyr <- syr+4
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  
  # Figure ??  5-yr warm period summer and winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_warm_us_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr warm periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('e)', 'f)','g)', 'h)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods US












if (warmcold_eu) {
#  # find and mark coldest/warmest 5yr period in EU
#   d <- aind.allts$ensmean[(which(eind.allts$names=="EU.temp2")),
#                           (is.even(seq(1,ncol(aind.allts$ensmean))))]
#   dec <- rep(NA,250)
#   dectime <- rep(NA,250)
#   for (t in 1:length(dec)) {
#     if (t==1) {
#       dec[t] <- mean(d[1:5])
#       dectime[t] <- 1603
#       i=2
#     } else {
#       dec[t] <- mean(d[i:(i+4)])
#       dectime[t] <- dectime[t-1] + 1 
#       i=i+1
#     }
#   } 
#   #  plot(dectime,dec,ty='l',col='red')
#   tmp=cbind(dec,dectime)
#   tmp[order(dec),]  
 


  coldper <- c(1695, 1641, 1674, 1739, 1713, 1701)
  for (syr in coldper) {
#    # calc validate clim (70 yr period) to substract for anomalies
#     if (syr>1645) {asyr<-syr-35} else {asyr=1610}
#     for (cyr in seq(asyr,(syr+39))) {  
#       load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
#       if (cyr==asyr){
#         vtmp <- validate$data
#         c=1
#       } else {
#         vtmp <- vtmp+validate$data
#         c <- c+1
#       }
#     }
#     valclim <- vtmp/c  
    eyr <- syr+4
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
#        vtmp <- validate$data
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
#        vtmp ...        
        c=c+1
      }
    }
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
#      valmean <- vtmp/c
#      valmeananom <- valmean-valclim
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
#      valmean <- vtmp/c
#      valmeananom <- cbind(valmeananom,(valmean-valclim))
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)

  # Figure ??  5-yr cold period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    #   t=2 # summer
    pdf(paste0('../figures/figure_cold_eu_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
#     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
#                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
#     levs <- seq(-4,4,1)
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
#                 latlim=c(30,75),lonlim=c(-20,40),units="K",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Temp. reconstruction',colnames=rep('',4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr cold periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
#     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
#                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
#     levs <- seq(-4,4,1)
#     plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
#                 latlim=c(30,75),lonlim=c(-20,40),units="mm",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Prec. reconstruction',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop

  
  
  
    
  warmper <- c(1794, 1800, 1846, 1826, 1819, 1779)
  for (syr in warmper) {
    eyr <- syr+4
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  
  # Figure ??  5-yr warm period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_warm_eu_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr warm periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods EU











if (warmcold_nh) {
#  # find and mark coldest/warmest 5yr period in NH
#   d <- aind.allts$ensmean[(which(eind.allts$names=="NH.temp2")),
#                           (is.even(seq(1,ncol(aind.allts$ensmean))))]
#   dec <- rep(NA,250)
#   dectime <- rep(NA,250)
#   for (t in 1:length(dec)) {
#     if (t==1) {
#       dec[t] <- mean(d[1:5])
#       dectime[t] <- 1603
#       i=2
#     } else {
#       dec[t] <- mean(d[i:(i+4)])
#       dectime[t] <- dectime[t-1] + 1 
#       i=i+1
#     }
#   } 
#   #  plot(dectime,dec,ty='l',col='red')
#   tmp=cbind(dec,dectime)
#   tmp[order(dec),]  
#   
  
  
#  coldper <- c(1641, 1695, 1676, 1815, 1667, 1701, 1738)
  coldper <- c(1641, 1695, 1815, 1832) #, 1676, 1667)
  for (syr in coldper) {
    eyr <- syr+4
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                       na.string="NaN")
    ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                     na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)
  colnames(mannmeananom) <- coldper
  # Figure ??  5-yr cold period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_cold_nh_',s,'_v5.pdf'), width=16, 
        height=12, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                             mannmeananom[,4]), c(nrow(plotmann$data),1,4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    #    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
    #                lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')
    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main='EKF400 5-yr cold periods', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
  
  
  
  
  
  warmper <- c(1801, 1826, 1791, 1778, 1844, 1779)
  for (syr in warmper) {
    eyr <- syr+4
        # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                       na.string="NaN")
    ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                     na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  colnames(mannmeananom) <- warmper
  # Figure ??  5-yr warm period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_warm_nh_',s,'_v2.pdf'), width=16, 
        height=12, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                          mannmeananom[,4]), c(nrow(plotmann$data),1,4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
#    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
#                lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')
    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main='EKF400 5-yr warm periods', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods NH




if (dryhumid_us) {
#  # find and mark driest/wettest 5yr period in US
#   d <- aind.allts$ensmean[(which(eind.allts$names=="NAM.precip")),
#                           (is.even(seq(1,ncol(aind.allts$ensmean))))]
#   dec <- rep(NA,250)
#   dectime <- rep(NA,250)
#   for (t in 1:length(dec)) {
#     if (t==1) {
#       dec[t] <- mean(d[1:5])
#       dectime[t] <- 1603
#       i=2
#     } else {
#       dec[t] <- mean(d[i:(i+4)])
#       dectime[t] <- dectime[t-1] + 1 
#       i=i+1
#     }
#   } 
#   #  plot(dectime,dec,ty='l',col='red')
#   tmp=cbind(dec,dectime)
#   tmp[order(dec),]  


  dryper <- c(1819, 1843, 1635, 1671, 1760, 1743)
  for (syr in dryper) {
    eyr <- syr+4
    # load nada pdsi recon
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==dryper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end dry periods loop  
  colnames(anameananom) <- rep(dryper,each=2)
  colnames(echmeananom) <- rep(dryper,each=2)

  # Figure ??  5-yr cold period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_dry_us_',s,'_v1.pdf'), width=12, 
        height=14, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=dryper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr dry periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop

  
  humidper <- c(1688, 1729, 1656, 1802, 1698) #, 1603)
  for (syr in humidper) {
    eyr <- syr+4
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==humidper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end humid periods loop  
  colnames(anameananom) <- rep(humidper,each=2)
  colnames(echmeananom) <- rep(humidper,each=2)
  
  # Figure ??  5-yr warm period summer and winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_humid_us_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=humidper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr humid periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('e)', 'f)','g)', 'h)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop  
} # end dry, humid periods US










if (dryhumid_eu) {
#  # find and mark driest/wettest 5yr period in EU
#   d <- aind.allts$ensmean[(which(eind.allts$names=="EU.precip")),
#                           (is.even(seq(1,ncol(aind.allts$ensmean))))]
#   dec <- rep(NA,250)
#   dectime <- rep(NA,250)
#   for (t in 1:length(dec)) {
#     if (t==1) {
#       dec[t] <- mean(d[1:5])
#       dectime[t] <- 1603
#       i=2
#     } else {
#       dec[t] <- mean(d[i:(i+4)])
#       dectime[t] <- dectime[t-1] + 1 
#       i=i+1
#     }
#   } 
#   #  plot(dectime,dec,ty='l',col='red')
#   tmp=cbind(dec,dectime)
#   tmp[order(dec),]  
  
  
  dryper <- c(1846, 1823, 1808, 1837, 1816, 1794)
  for (syr in dryper) {
    #    # calc validate clim (70 yr period) to substract for anomalies
    #     if (syr>1645) {asyr<-syr-35} else {asyr=1610}
    #     for (cyr in seq(asyr,(syr+39))) {  
    #       load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
    #       if (cyr==asyr){
    #         vtmp <- validate$data
    #         c=1
    #       } else {
    #         vtmp <- vtmp+validate$data
    #         c <- c+1
    #       }
    #     }
    #     valclim <- vtmp/c  
    eyr <- syr+4
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        #        vtmp <- validate$data
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        #        vtmp ...        
        c=c+1
      }
    }
    if (syr==dryper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      #      valmean <- vtmp/c
      #      valmeananom <- valmean-valclim
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      #      valmean <- vtmp/c
      #      valmeananom <- cbind(valmeananom,(valmean-valclim))
    }
  } # end dry periods loop  
  colnames(anameananom) <- rep(dryper,each=2)
  colnames(echmeananom) <- rep(dryper,each=2)
  
  # Figure ??  5-yr dry period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    #   t=2 # summer
    pdf(paste0('../figures/figure_dry_eu_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    #     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
    #                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    #     levs <- seq(-4,4,1)
    #     plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
    #                 latlim=c(30,75),lonlim=c(-20,40),units="K",
    #                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
    #                 st.col=NULL, stations=NULL, add=T,addcontours=F,
    #                 wcol='darkgrey',rownames='Temp. reconstruction',colnames=rep('',4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=dryper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr dry periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    #     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
    #                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    #     levs <- seq(-4,4,1)
    #     plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
    #                 latlim=c(30,75),lonlim=c(-20,40),units="mm",
    #                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
    #                 st.col=NULL, stations=NULL, add=T,addcontours=F,
    #                 wcol='darkgrey',rownames='Prec. reconstruction',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
  
  
  
  
  
  humidper <- c(1617, 1664, 1638, 1701, 1713, 1625)
  for (syr in humidper) {
    eyr <- syr+4
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0('../data/echam_clim/echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0('../data/echam_clim/echam_clim_',syr+2,'-',(syr+3),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotpath,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==humidper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
    }
  } # end humid periods loop  
  colnames(anameananom) <- rep(humidper,each=2)
  colnames(echmeananom) <- rep(humidper,each=2)
  
  # Figure ??  5-yr humid period winter
  pdata <- echam
  for (t in 1:2){
    if (t==1){s="win"} else {s="som"}
    pdf(paste0('../figures/figure_humid_eu_',s,'_v1.pdf'), width=12, 
        height=10.5, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=humidper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr humid periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,t]-echmeananom[,t]),
                          (anameananom[,t+2]-echmeananom[,t+2]),
                          (anameananom[,t+4]-echmeananom[,t+4]),
                          (anameananom[,t+6]-echmeananom[,t+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,t],anameananom[,t+2],anameananom[,t+4],
                          anameananom[,t+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end dry, humid periods EU





if (jetshift){
  if (land_only) {
    load(paste0(prepplotpath,'indices_allts_landonly.Rdata'))
  } else {
    load(paste0(prepplotpath,'indices_allts.Rdata'))
  }
  # Stefans ideas: 
  # decadal shifts of subtropical belt. 1945 northern position, 
  # until 1980 southward movement, since 1980 northward trend
  if (land_only) {
    pdf(file='../figures/indices/N_hadley_cell_extend_landonly.pdf',width=10,height=14)
  } else {
    pdf(file='../figures/indices/N_hadley_cell_extend.pdf',width=10,height=14)
  }
  par(oma=c(1.5,1,1,0))
  par(mfrow=c(3,1))  
  wl=11 # 31yr low-pass filter
  ts=1:402
  anomper <- (1901:1980)
  seas1 <- 'yrmean'
  seas2 <- 'win'
  seas3 <- 'sum'
  reg1 <- 37 # for SJ_u200.calc
  reg2 <- 38 # for SJ_slp.calc
  reg3 <- 46 # for SJ index
  period <- eind.allts$time #seq(syr,eyr)
  
  stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
               header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
               seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
  head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
            nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
            seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
  # (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
  newhead <- rep(NA,ncol(head))
  for (i in 1:ncol(head)) {
    newhead[i] <- paste(head[1,i],head[2,i],sep='_')
  }
  colnames(stefind) <- newhead
  #load('../data/indices_recon.Rdata')
  vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]
  
  for (seas in c(seas1,seas2,seas3)) {
    for (reg in c(reg1,reg2,reg3)) {
      if (seas=='yrmean') {
        years <- rep(eind.allts$time,each=2)  
      } else {
        wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
        somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
        wincolv <- is.odd(seq(1,ncol(vind.allts$data)))
        somcolv <- is.even(seq(1,ncol(vind.allts$data)))    
      }
      if (seas=='yrmean') {
        print("annual mean")
        eindmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
        aindmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
        eindbem <- aggregate(eind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
        aindbem <- aggregate(aind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
        eindmin <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
        eindmax <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
        aindmin <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
        aindmax <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
      } else if (seas=='sum'){  
        print("summer season")
        eindmean <- eind.allts$ensmean[reg,somcol]  
        aindmean <- aind.allts$ensmean[reg,somcol]  
        eindbem <- eind.allts$bem[reg,somcol]
        aindbem <- aind.allts$bem[reg,somcol]
        eindmin <- apply(eind.allts$data[reg,somcol,],1,min,na.rm=T)
        eindmax <- apply(eind.allts$data[reg,somcol,],1,max,na.rm=T)
        aindmin <- apply(aind.allts$data[reg,somcol,],1,min,na.rm=T)
        aindmax <- apply(aind.allts$data[reg,somcol,],1,max,na.rm=T)
      } else if (seas=='win'){
        print("winter season")
        eindmean <- eind.allts$ensmean[reg,wincol]  
        aindmean <- aind.allts$ensmean[reg,wincol]  
        eindbem <- eind.allts$bem[reg,wincol]
        aindbem <- aind.allts$bem[reg,wincol]
        eindmin <- apply(eind.allts$data[reg,wincol,],1,min,na.rm=T)
        eindmax <- apply(eind.allts$data[reg,wincol,],1,max,na.rm=T)
        aindmin <- apply(aind.allts$data[reg,wincol,],1,min,na.rm=T)
        aindmax <- apply(aind.allts$data[reg,wincol,],1,max,na.rm=T)
      }
      print("join indices to allind")
      allind <- cbind(eindmean,aindmean,eindbem,aindbem,eindmin,eindmax,aindmin,aindmax)
      anopos <- period %in% anomper
      centerfac <- apply(allind[anopos,],2,mean)
      centerfac[5] <- centerfac[6] <- centerfac[3] <- centerfac[1]
      centerfac[7] <- centerfac[8] <- centerfac[4] <- centerfac[2]
      #  if (sca) {
      #    scalefac <- apply(allind[anopos,],2,sd)
      #    scalefac[5] <- scalefac[6] <- scalefac[3] <- scalefac[1]
      #    scalefac[7] <- scalefac[8] <- scalefac[4] <- scalefac[2]
      #  } else { scalefac <- rep(1,8) }  
      allindann <- allind
      if (wl > 1) {
        for (i in 1:ncol(allind)){
          allind[,i] <- runmean(allind[,i],wl) #scalefac),wl)
        }
        #  allind <- runmean(allind,wl) #scalefac),wl)
      }
      #print("scale allind")
      #if (wl > 1) {
      #  allind <- runmean(scale(allind,center=centerfac,scale=F),wl) #scalefac),wl)
      #} else {
      #  allind <- scale(allind,center=centerfac,scale=F) #scalefac)
      #}
      if (reg==reg1){allind1 <- allind}
      if (reg==reg2){allind2 <- allind}
      if (reg==reg3){allind3 <- allind}
    } # end reg loop
    # plot all
    if (seas==seas1){
      if (land_only){
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="",
             ylim=c(35,48),main='Shifts in the northern Hadley cell extend (annual mean)')  
      } else {
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="",
          ylim=c(35,39),main='Shifts in the northern Hadley cell extend (annual mean)')
      }
      abline(h=c(mean(allind1[ts,1]),mean(allind2[ts,1])),col='lightgrey')
    } else if (seas==seas2){
      if (land_only){
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="",
           ylim=c(29,45),main='Shifts in the northern Hadley cell extend (winter)')
      } else {
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="",
             ylim=c(32,35),main='Shifts in the northern Hadley cell extend (winter)')
      }
      abline(h=c(mean(allind1[ts,1]),mean(allind2[ts,1])),col='lightgrey')
    } else if (seas==seas3){
      if (land_only){
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="year",
           ylim=c(41,52),main='Shifts in the northern Hadley cell extend (summer)')
      } else {
        plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="lat",xlab="year",
             ylim=c(38,43),main='Shifts in the northern Hadley cell extend (summer)')
      }
      abline(h=c(mean(allind1[ts,1]),mean(allind2[ts,1])),col='lightgrey')
    } 
    lines(period[ts],allind1[ts,2],ty='l',lwd=1,lty=1,col=rgb(10,0,0,10,maxColorValue=10), ylab="",xlab="")
    lines(period[ts],allind2[ts,1],ty='l',lwd=1,lty=2,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
    lines(period[ts],allind2[ts,2],ty='l',lwd=1,lty=2,col=rgb(10,0,0,10,maxColorValue=10), ylab="",xlab="")
#    lines(period[ts],allind3[ts,1],ty='l',lwd=1,lty=2,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
#    lines(period[ts],allind3[ts,2],ty='l',lwd=1,lty=2,col=rgb(10,0,0,10,maxColorValue=10), ylab="",xlab="")
    if (seas==seas1){
      legend('topleft', c('echam zonal mean u200 max','anal. u200 max','echam slp max','anal. slp max'), bg='white', 
           box.col='white',lty=c(1,1,2,2),col=c(1,2,1,2))
    }
#     vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
#     yrv <- stefind[,1]
#     centerfac <- apply(vdata,2,mean,na.rm=T)
#     scalefac <- apply(vdata,2,sd,na.rm=T)
#     #if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=F,scale=F),wl) #centerfac,scale=scalefac),wlen)
#     #} else {
#     #  vdata <- scale(vdata,center=centerfac,scale=scalefac)
#     #}
#     lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
    
    # ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
    # acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
    # ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
    # acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
    # legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],
    #   ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],
    #   acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],
    #   ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],
    #   acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
  }
  dev.off()
} # end jet shift








# plot AMO index in warm/cold periods
#load(paste0(prepplotpath,'.Rdata'))
load(paste0(prepplotpath,'indices_allts_landonly.Rdata'))
v <- read.table('/Volumes/DATA/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]
# 95% quantile of largest eruptions
vf <- v[600:999,1][which(v[600:999,2]<quantile(v[600:999,2],0.05))]
nc=open.nc('../data/pdo_amo/joerg/amo_yr.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
close.nc(nc)
amo_j <- cbind(trunc(tim/100),tmp)

nc=open.nc('../data/pdo_amo/joerg_v1/pdo_yr.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
close.nc(nc)
pdo_j <- cbind(trunc(tim/100),tmp)

nc=open.nc('../data/pdo_amo/joerg_v2/eofcoeff_a.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst_res"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- trunc(tlist/100)
close.nc(nc)
pdo_j2 <- cbind(tim,tmp)

tmp <- read.table('../data/pdo_amo/sina_pdo/PDO.txt',header=T)
tmp2 <- aggregate(tmp[,3],list(tmp[,1]),mean)
pdo_s <- tmp2[which(tmp2[,1]>1599),]

pdo_m <- read.table('../data/pdo_amo/clim_expl/ipdo_mann_1600_2000.dat',header=F)

tmp <- as.matrix(read.table('../data/pdo_amo/clim_expl/ipdo_hadsst3.dat'))
pos <- which(tmp < -999)
tmp[pos] <- NA
pdo_h <- cbind(tmp[,1],apply(tmp[,2:13],1,mean))

nc=open.nc('../data/pdo_amo/sina_amo/amo.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
tmp2 <- cbind(trunc(tim/100),tmp)
amo_s <- tmp2[tmp2[,1]>1599,]

amo_m <- read.table('../data/pdo_amo/clim_expl/iamo_mann_1600_2000.dat',header=F)

tmp <- as.matrix(read.table('../data/pdo_amo/clim_expl/iamo_hadsst.dat'))
pos <- which(tmp < -999)
tmp[pos] <- NA
amo_h <- cbind(tmp[,1],apply(tmp[,2:13],1,mean))

pdf(file='../figures/indices/AMO_PDO_JF_Sina_Mann_HadSST.pdf',width=8,height=4)
par(oma=c(1.5,1,1,0.1),mar=c(1.5,1,1,0.1))
par(mfrow=c(2,1))  
plot(pdo_j[,1],scale(runmean(pdo_j[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
     ylim=c(-3,3),xlim=c(1600,2000),cex=0.5,cex.axis=0.5,main='PDO',xaxt='n',
     col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_j2[,1],scale(runmean(-pdo_j2[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(pdo_s[,1],scale(runmean(-pdo_s[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
#      col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_m[,1],scale(runmean(pdo_m[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_h[,1],scale(runmean(pdo_h[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(0,0,0,7,maxColorValue=10), ylab="",xlab="")
cold <- rep(c(1641,1646,1695,1700,1815,1820,1832,1837),each=2) #*100
warm <- rep(c(1778,1783,1791,1796,1801,1806,1826,1831),each=2) #*100
minmax <- c(-3,3,3,-3)
for (i in c(1,5,9,13)) {
  polygon(cold[i:(i+3)],minmax,col=rgb(0,0,10,3,maxColorValue=10),border=NA)
}
for (i in c(1,5,9,13)) {
  polygon(warm[i:(i+3)],minmax,col=rgb(10,0,0,3,maxColorValue=10),border=NA)
}
abline(v=vf,col='black',lty=2)
legend('topleft',c('PDO Joerg','PDO Sina','PDO Mann','PDO HadSST'),lty=1,cex=0.5,
       col=c(rgb(0,10,10,7,maxColorValue=10),rgb(10,0,10,7,maxColorValue=10),
             rgb(0,0,10,7,maxColorValue=10),rgb(0,0,0,7,maxColorValue=10)),
       bg='white', box.col='white')

plot(amo_j[,1],scale(runmean(amo_j[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
     ylim=c(-3,3),xlim=c(1600,2000),cex=0.5,cex.axis=0.5,main='AMO',
     col=rgb(10,0,0,7,maxColorValue=10), ylab="",xlab="")
lines(amo_s[,1],scale(runmean(amo_s[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,7,7,7,maxColorValue=10), ylab="",xlab="")
lines(amo_m[,1],scale(runmean(amo_m[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(amo_h[,1],scale(runmean(amo_h[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#years <- rep(eind.allts$time,each=2)  
##reg=which(aind.allts$names=="PNA.calc")
#reg=which(aind.allts$names=="ENH.temp2")
#eindyrmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
#aindyrmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
##lines(eind.allts$time*100,scale(runmean(eindyrmean,1),center=T,scale=2),
##      lwd=2,lty=1,col=rgb(0,0,0,5,maxColorValue=10))
##lines(aind.allts$time*100,scale(runmean(aindyrmean,1),center=T,scale=2),
##      lwd=2,lty=1,col=rgb(10,0,0,5,maxColorValue=10))
#minmax <- c(-3,3,3,-3)
for (i in c(1,5,9,13)) {
  polygon(cold[i:(i+3)],minmax,col=rgb(0,0,10,3,maxColorValue=10),border=NA)
}
for (i in c(1,5,9,13)) {
  polygon(warm[i:(i+3)],minmax,col=rgb(10,0,0,3,maxColorValue=10),border=NA)
}
#abline(v=vf*100,col='black',lty=2)
abline(v=vf,col='black',lty=2)
#legend("topleft",c("AMO","PDO"),col=c(rgb(0,10,10,10,maxColorValue=10),
#       rgb(10,0,10,10,maxColorValue=10)),lwd=c(2,2),lty=c(1,1), 
#       cex=1, bty='o', bg='white', box.col='white')
legend('topleft',c('AMO Joerg','AMO Sina','AMO Mann','AMO HadSST'),lty=1,cex=0.5,
       col=c(rgb(10,0,0,7,maxColorValue=10),rgb(10,7,7,7,maxColorValue=10),
             rgb(10,10,0,7,maxColorValue=10),rgb(10,0,10,7,maxColorValue=10)),
       bg='white', box.col='white')
dev.off()









###################################################################
# plot ECHAM global mean temp. vs. scaled forcings ################
###################################################################
# load("/Volumes/DATA/climdata/PAGES_DB_Raphi/ProxyDB_1_3_3.RData")
forc <- read.table("../data/forcings/ccc400_forcings.txt",header=T)
pdf(file='../figures/figure2_glo_mean_scal_forc2.pdf',width=5,height=3.5)
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(forc[1:405,1],forc[1:405,7],ty='l',col='red',lwd=2,xlab="Year",
     ylab="Temperature anomaly [ºC]",ylim=c(-0.9,1.1))
lines(forc[1:405,1],forc[1:405,6],col='blue',lwd=2)
#lines(forc[250:405,1],forc[250:405,8],col='black',lwd=2)
colnames(forc)
reg <- lm(forc[250:405,7]~forc[250:405,2]+forc[250:405,3]+forc[250:405,4]+forc[250:405,5])
summary(reg)
ti <- seq(1601,2005)
fit <- reg$coefficients[1]+reg$coefficients[2]*forc[1:405,2]+reg$coefficients[3]*forc[1:405,3]+
       reg$coefficients[4]*forc[1:405,4]+reg$coefficients[5]*forc[1:405,5]
lines(ti,fit,col='black',lwd=2)
legend('topleft', c('CRUTEM4','CCC400','scaled forcings'),
       lwd=c(2,1,1),lty=c(1,1,1), col=c('red','blue','black'),
       cex=1, bty='o', bg='white',box.col='white')
dev.off()








if (landusebugbias){
  figpath="../figures/land_bug/"
  # 1605-1635 and 1971-2000 temp2 bias
  for (yr in c(1605,1971)){
    syr <- yr
    eyr <- yr+29
    print(c(syr,eyr))
    # 1. load echam mem 103 with corrected land use
    i=0
    emean103 <- array(NA,dim=c(4608,30))
    for (cyr in syr:eyr) {
      i=i+1
      load(paste0(datadir,"EnSRF_analysis/echam_103/echam_",(cyr-1),"-",(cyr),"_2ndgrid.Rdata"))
      emean103[,i] <- apply(echam$data[which(echam$names=="temp2"),13:24,1],1,mean)
    } 
    emean103_period <- apply(emean103,1,mean)
    # 2. load analysis (oct syr-1 - sep eyr)
    amean <- array(NA,dim=c(4608,30))
    emean <- array(NA,dim=c(4608,30))
    i=0
    for (cyr in syr:eyr) {
      i=i+1
      load(paste0(datadir,"EnSRF_analysis/prepplot_v3_seasonal/analysis_",cyr,".Rdata"))
      amean[,i] <- apply(analysis$ensmean[which(echam$names=="temp2"),],1,mean)
      emean[,i] <- apply(echam$ensmean[which(echam$names=="temp2"),],1,mean)
    } 
    amean_period <- apply(amean,1,mean)
    emean_period <- apply(emean,1,mean)
    # calc bias
    bias_a_e103_period <- amean_period-(emean103_period-273.15)
    bias_e_e103_period <- emean_period-(emean103_period-273.15)
    # plot bias
    plotbias <- echam
    plotbias$ensmean <- cbind(bias_e_e103_period,bias_a_e103_period) #,dim=c(4608,2,1))
    plotbias$lon <- echam$lon[which(echam$names=="temp2")]
    plotbias$lat <- echam$lat[which(echam$names=="temp2")]
    plotbias$names <- echam$names[which(echam$names=="temp2")]
    pdf(paste(figpath,"bias_anal-ech103_",syr,"-",eyr,".pdf",sep=''), 
        width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,4,0))
    levs <- c(-Inf, -3,-2,-1,-0.5, -0.2, 0.2, 0.5, 1, 2, 3, Inf)
    plot_echam(plotbias, varname='temp2', type='ensmean', lty=3, lev=levs, 
               ti=1:2, colnames=c("echam mean - 103", "analysis mean - 103"), 
               main="Land use bias", add=T)
    dev.off()
  }
  # result: bias does not improve with assimilation because I only work with anomalies 
  # and add the climatology back on
}








