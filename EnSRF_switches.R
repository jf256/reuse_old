#####################################################################################
# general switches 
#####################################################################################
machine="macbook" #"climcal3" # "climpa12"  
if (machine=="macbook") {
  datadir="/Volumes/DATA/climdata/"
  workdir=('~/unibe/projects/EnSRF/src/')
} else {
  datadir="/scratch/joerg/climdata/"
  workdir=('/scratch/${USER}/projects/reuse/src/')
}
setwd(workdir)
#!!!ATTENTION: sixmonstatevector year starts in October of previous year (cyr-1)
sixmonstatevector=T    # 6 months of data in state vector for real proxy multiple 
# regression approach. ATTENTION: monthly has to be TRUE
if (sixmonstatevector) {
  s <- 2 # only 2 seasons to calculate but still monthly results in long state vector
} else {
  s <- 12
}
season=c(3,9)  # 3,9 = apr-sep and oct-mar, num=end month of season; # season=c(2,5,8,11)
nmem=30        # number of ensemble members

#####################################################################################
# EnSRF_data switches
#####################################################################################
# load or generate data from scratch
generate_ECHAM=F       # if TRUE -> orig. echam data is read
generate_ECHAM_1901_70=F # ECHAM data for bias calc with real_prox data 
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
generate_ind_recon=F   # read Stefan's indices 1900-2000 from .txt to .RData
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
generate_DOCUM=F        # if TRUE -> yuri's docu. data collection is read 

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
  docum=F                 # read documentary based data
} else {
  docum=T
}
print(paste("instr:",instrumental, "proxies:",real_proxies, "documentary:",docum))


calc_decorr_dist=F     # calculate decorrelation distance for each variable from ECHAM to set L
l_dist_temp2=1000*1.5  # factor *1.5 after stefans recommendation
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
# l_dist_ind=999999 # precalculated indices should be removed

# how to treat multiple input series in same grid box
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
no_stream=T            # all echam vars but stream function as there is problem with 
                       # 5/9 levels, which are in lat dimension before and after 1880
load_71yr_anom=T       # load 71yr echam anomalies calculated with cdo
anom_reload=F          # reload anom calculated in R (next option)
anom_save=F            # save anom calculated in R to reload next time
if (load_71yr_anom==T) {
  anom_reload=F
  anom_save=F}
check_assimdata=T      # screen assimilation data before using it

if (no_stream & tps_only) {
  tps_only = F
  print('ACHTUNG: tps_only was set to FALSE')
}


# other options
scaleprox=T            # scale standardized docu and prox data the echam variance at location
anomaly_assim=T        # work with anomalies to avoid reg. konst in state vector
nseas <- 12            # year with 12 months
check_dist=F           # test for ideal cut-off distance of spatial correlations
#H_non_lin=F           # new H operator that also allows non-linear functions
ana.enssize=F
NCEP_SOCOL=F

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

#####################################################################################
# prepare plot switches
#####################################################################################
monthly_out=F    # if sixmonstatevector=T output is backtransformed to seasonal 
                 # average or monthly data if monthly_out=T 
if (monthly_out) {
  prepplotdir=paste0(datadir,'ENSRF_analysis/prepplot_monthly/')
} else {
  prepplotdir=paste0(datadir,'/Volumes/DATA/climdata/ENSRF_analysis/prepplot_v3_seasonal/')  
}

save_prepplot=F  # save half year averages calc from monthly data into /prepplot folder
load_prepplot=T  # ATTENTION check if folder prepplot on scratch contains monthly or seasonal data!
load_image=T     # directly load image for syr-eyr period: 1902-2001 or 1651-1750 image
write_coor=T     # write ascii files with assimilated stations and data per ts


#####################################################################################
# plot switches
#####################################################################################
monthly=F
pseudoproxy=F
plot_dweights=F
write_nc=F
recalc <- F
reload <- F
plstat <- NULL #calibrate # NULL or calibrate
#PAGES <- F          # write output for PAGES paper



# -----------------------------------------------------------------------------------------
