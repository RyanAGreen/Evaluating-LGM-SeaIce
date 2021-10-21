import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.path as mpath
import matplotlib.lines as mlines

# load in PMIP3 data
CNRM_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/sic_OIclim_CNRM-CM5_lgm_r1i1p1_180001-199912-climregrid.nc',decode_times=False)
FGOALS_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/sic_OIclim_FGOALS-g2_lgm_r1i1p1_055001-064912-climregrid.nc',decode_times=False)
IPSL_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/sic_OIclim_IPSL-CM5A-LR_lgm_r1i1p1_260101-280012-climregrid.nc',decode_times=False)
MIROC_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/sic_OIclim_MIROC-ESM_lgm_r1i1p1_460001-469912-climregrid.nc',decode_times=False)
MRI_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/sic_OIclim_MRI-CGCM3_lgm_r1i1p1_250101-260012-climregrid.nc',decode_times=False)

LOVE1_sum_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/LOVE1_albq_summer.nc',decode_times=False) # check if this is loveclim1 or 2
feb = LOVE1_sum_3.FEB.mean(dim='AX005')
mar = LOVE1_sum_3.MAR.mean(dim='AX006')
LOVE1_sum_3 = (mar + feb)/2

LOVE1_win_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/V3LNAwSeaIceConcWinter.nc',decode_times=False) # check if this is loveclim1 or 2
jul = LOVE1_win_3.ALJUL.mean(dim='AX007')
aug2 = LOVE1_win_3.ALAUG.mean(dim='AX008')
LOVE1_win_3 = (jul + aug2)/2

#### Summer

# data sets with preset variables
CCSM4_sum_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMinimum-sic-CCSM4_lgm.nc',decode_times=False)
CCSM4_sum_3 = (CCSM4_sum_3.S1A + CCSM4_sum_3.S1B)/2
GISS_sum_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMinimum-sic-GISS_lgm.nc',decode_times=False)
GISS_sum_3 = (GISS_sum_3.S2A + GISS_sum_3.S2B)/2
MPI_sum_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMinimum-sic-MPI_lgm.nc',decode_times=False)
MPI_sum_3 = (MPI_sum_3.S3A + MPI_sum_3.S3B)/2


# original data file that has 12 time points
CNRM_3_sum = xr.concat([CNRM_3.sel(time=(CNRM_3.time[1])), CNRM_3.sel(time=(CNRM_3.time[2]))], dim="time")
CNRM_3_sum = CNRM_3_sum.mean(dim='time') # or whatever the time axis is called for this variable (different for loveclim)
FGOALS_3_sum = xr.concat([FGOALS_3.sel(time=(FGOALS_3.time[2])), FGOALS_3.sel(time=(FGOALS_3.time[3]))], dim="time")
FGOALS_3_sum = FGOALS_3_sum.mean(dim='time')
IPSL_3_sum = xr.concat([IPSL_3.sel(time=(IPSL_3.time[1])), IPSL_3.sel(time=(IPSL_3.time[2]))], dim="time")
IPSL_3_sum = IPSL_3_sum.mean(dim='time')
MIROC_3_sum = xr.concat([MIROC_3.sel(time=(MIROC_3.time[1])), MIROC_3.sel(time=(MIROC_3.time[2]))], dim="time")
MIROC_3_sum = MIROC_3_sum.mean(dim='time')
MRI_3_sum = xr.concat([MRI_3.sel(time=(MRI_3.time[1])), MRI_3.sel(time=(MRI_3.time[2]))], dim="time")
MRI_3_sum = MRI_3_sum.mean(dim='time')
LOVE2_sum_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/V3LNAwSOwSHWwSeaIceConcSUMMER.nc',decode_times=False)  # check if this is loveclim1 or 2
aug = LOVE2_sum_3.ALFEB.mean(dim='AX005')
sep = LOVE2_sum_3.ALMAR.mean(dim='AX006')
LOVE2_sum_3 = (aug + sep)/2



# Multi model mean
MMM_sum = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/AnnualMinimum15%.nc',decode_times=False)

#### Winter

# data sets with preset variables
# data sets with preset variables
CCSM4_win_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMaximum-sic-CCSM4_lgm.nc',decode_times=False)
CCSM4_win_3 = (CCSM4_win_3.S1A + CCSM4_win_3.S1B)/2
GISS_win_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMaximum-sic-GISS_lgm.nc',decode_times=False)
GISS_win_3 = (GISS_win_3.S2A + GISS_win_3.S2B)/2
GISS_win_3 = GISS_win_3.isel(LAT=(GISS_win_3.LAT > -70))
MPI_win_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/2MonthMaximum-sic-MPI_lgm.nc',decode_times=False)
MPI_win_3 = (MPI_win_3.S3A + MPI_win_3.S3B)/2


# original data file that has 12 time points
CNRM_3_win = xr.concat([CNRM_3.sel(time=(CNRM_3.time[8])), CNRM_3.sel(time=(CNRM_3.time[9]))], dim="time")
CNRM_3_win = CNRM_3_win.mean(dim='time') # or whatever the time axis is called for this variable (different for loveclim)
FGOALS_3_win = xr.concat([FGOALS_3.sel(time=(FGOALS_3.time[8])), FGOALS_3.sel(time=(FGOALS_3.time[9]))], dim="time")
FGOALS_3_win = FGOALS_3_win.mean(dim='time')
FGOALS_3_win = FGOALS_3_win.isel(lat=(FGOALS_3_win.lat > -70))
IPSL_3_win = xr.concat([IPSL_3.sel(time=(IPSL_3.time[7])), IPSL_3.sel(time=(IPSL_3.time[8]))], dim="time")
IPSL_3_win = IPSL_3_win.mean(dim='time')
MIROC_3_win = xr.concat([MIROC_3.sel(time=(MIROC_3.time[8])), MIROC_3.sel(time=(MIROC_3.time[9]))], dim="time")
MIROC_3_win = MIROC_3_win.mean(dim='time')
MRI_3_win = xr.concat([MRI_3.sel(time=(MRI_3.time[8])), MRI_3.sel(time=(MRI_3.time[9]))], dim="time")
MRI_3_win = MRI_3_win.mean(dim='time')
MRI_3_win = MRI_3_win.isel(lat=(MRI_3_win.lat > -70))
LOVE2_win_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/V3LNAwSOwSHWwSeaIceConcWINTER.nc',decode_times=False)  # check if this is loveclim1 or 2
aug = LOVE2_win_3.ALAUG.mean(dim='AX005')
sep = LOVE2_win_3.ALSEP.mean(dim='AX006')
LOVE2_win_3 = (aug + sep)/2


# Multi model mean winter
MMM_win = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/AnnualMaximum15%.nc',decode_times=False)

# load in PMIP4 data

iLOVECLIM = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/cresum18250_regrid.nc',decode_times=False)
LOVEsic = iLOVECLIM.albq
feb = LOVEsic[1:2400:12]
mar = LOVEsic[2:2400:12]
aug = LOVEsic[7:2400:12]
sep = LOVEsic[8:2400:12]
mar = mar.mean(dim='time')
feb = feb.mean(dim='time')
sep = sep.mean(dim='time')
aug = aug.mean(dim='time')
LOVEsummer = (mar+feb)/2
LOVEwinter = (sep + aug)/2
LOVEsummer = LOVEsummer * 100
LOVEwinter = LOVEwinter * 100


AWI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/AWI_PMIP4_siconca_regrid.nc',decode_times=False)
MIROC = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/MIROC_PMIP4_siconc_regrid.nc',decode_times=False)
MPI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/MPIPMIP4_siconc_regrid.nc',decode_times=False)
CESM12 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/b.e12.B1850C5.f19_g16.i21ka.03.pop.h.vars.08010900.climo_regrid.nc',decode_times=False)
CCSM4UoT = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/siconc_SImon_UofT-CCSM4_lgm_r1i1p1f1_gn_110101-120012_regrid.nc',decode_times=False)
IPSL = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/IPSLCM5A2_LGM_regrid.nc',decode_times=False)


# datasets that need to be averaged over 2 month period
AWI_winter = xr.concat([AWI.sel(time=(AWI.time[7::12])), AWI.sel(time=(AWI.time[8::12]))], dim="time") # August September
AWI_winter = AWI_winter.mean(dim='time')
AWI_winter = AWI_winter.isel(lat=(AWI_winter.lat > -70))
AWI_summer = xr.concat([AWI.sel(time=(AWI.time[1::12])), AWI.sel(time=(AWI.time[2::12]))], dim="time") # Feb March
AWI_summer = AWI_summer.mean(dim='time')

MIROC_winter = xr.concat([MIROC.sel(time=(MIROC.time[8::12])), MIROC.sel(time=(MIROC.time[9::12]))], dim="time") # Sep october
MIROC_winter = MIROC_winter.mean(dim='time')
MIROC_summer = xr.concat([MIROC.sel(time=(MIROC.time[1::12])), MIROC.sel(time=(MIROC.time[2::12]))], dim="time") # Feb march
MIROC_summer = MIROC_summer.mean(dim='time')

MPI_winter = xr.concat([MPI.sel(time=(MPI.time[8::12])), MPI.sel(time=(MPI.time[9::12]))], dim="time")
MPI_winter = MPI_winter.mean(dim='time')
MPI_summer = xr.concat([MPI.sel(time=(MPI.time[1::12])), MPI.sel(time=(MPI.time[2::12]))], dim="time")
MPI_summer = MPI_summer.mean(dim='time')

CESM12_winter = xr.concat([CESM12.sel(time=(CESM12.time[7])), CESM12.sel(time=(CESM12.time[8]))], dim="time")
CESM12_winter = CESM12_winter.mean(dim='time')
CESM12_summer = xr.concat([CESM12.sel(time=(CESM12.time[1])), CESM12.sel(time=(CESM12.time[2]))], dim="time")
CESM12_summer = CESM12_summer.mean(dim='time')

CCSM4UoT_winter = xr.concat([CCSM4UoT.sel(time=(CCSM4UoT.time[8::12])), CCSM4UoT.sel(time=(CCSM4UoT.time[9::12]))], dim="time")
CCSM4UoT_winter = CCSM4UoT_winter.mean(dim='time')
CCSM4UoT_summer = xr.concat([CCSM4UoT.sel(time=(CCSM4UoT.time[2::12])), CCSM4UoT.sel(time=(CCSM4UoT.time[3::12]))], dim="time")
CCSM4UoT_summer = CCSM4UoT_summer.mean(dim='time')

IPSL_winter = xr.concat([IPSL.sel(time_counter=(IPSL.time_counter[7])), IPSL.sel(time_counter=(IPSL.time_counter[8]))], dim="time_counter") # August September
IPSL_winter = IPSL_winter.mean(dim='time_counter')
IPSL_winter = IPSL_winter.isel(lat=(IPSL_winter.lat > -70))
IPSL_summer = xr.concat([IPSL.sel(time_counter=(IPSL.time_counter[1])), IPSL.sel(time_counter=(IPSL.time_counter[2]))], dim="time_counter") # Feb March
IPSL_summer = IPSL_summer.mean(dim='time_counter')

#####################

#####PMIP3#####
# austral summer MMM
MMM_CCSM4 = (MMM_sum.CCSM4R2 + MMM_sum.CCSM4R1)/2
MMM_GISS = (MMM_sum.GISSP150 + MMM_sum.GISSP151)/2
MMM_MPI = (MMM_sum.MPIP1 + MMM_sum.MPIP2)/2
MMM_sum['CCSM4'] = MMM_CCSM4
MMM_sum['GISS'] = MMM_GISS
MMM_sum['MPI'] = MMM_MPI
MMM_sum = MMM_sum.drop_vars('CCSM4R1')
MMM_sum = MMM_sum.drop_vars('CCSM4R2')
MMM_sum = MMM_sum.drop_vars('GISSP150')
MMM_sum = MMM_sum.drop_vars('GISSP151')
MMM_sum = MMM_sum.drop_vars('MPIP1')
MMM_sum = MMM_sum.drop_vars('MPIP2')
MMM_LOVE1_sum_3 = LOVE1_sum_3*100
MMM_LOVE2_sum_3 = LOVE2_sum_3*100

# multi model mean
MMM_sum_mean = (MMM_sum.CNRM + MMM_sum.FGOALS + MMM_sum.IPSL + MMM_sum.MIROC + MMM_sum.MRI +  MMM_sum.CCSM4 + MMM_sum.MPI + MMM_sum.GISS)/8

# calculating STD
var1=(MMM_sum.CCSM4-MMM_sum_mean)**2
var2=(MMM_sum.GISS-MMM_sum_mean)**2
var3=(MMM_sum.MPI-MMM_sum_mean)**2
var4=(MMM_sum.CNRM-MMM_sum_mean)**2
var5=(MMM_sum.FGOALS-MMM_sum_mean)**2
var6=(MMM_sum.IPSL-MMM_sum_mean)**2
var7=(MMM_sum.MIROC-MMM_sum_mean)**2
var8=(MMM_sum.MRI-MMM_sum_mean)**2
STD=((var1+var2+var3+var4+var5+var6+var7+var8)/8)**0.5
STD2=(MMM_sum_mean+STD)
STD3=(MMM_sum_mean-STD)

MMM_summer_PMIP3 = MMM_sum_mean
STD_MMM_summer_PMIP3 = STD2
MMM_summer_PMIP3_STD = STD3

# austral winter
MMM_CCSM4 = (MMM_win.CCSM4R2 + MMM_win.CCSM4R1)/2
MMM_GISS = (MMM_win.GISSP150 + MMM_win.GISSP151)/2
MMM_MPI = (MMM_win.MPIP1 + MMM_win.MPIP2)/2
MMM_win['CCSM4'] = MMM_CCSM4
MMM_win['GISS'] = MMM_GISS
MMM_win['MPI'] = MMM_MPI
MMM_win = MMM_win.drop_vars('CCSM4R1')
MMM_win = MMM_win.drop_vars('CCSM4R2')
MMM_win = MMM_win.drop_vars('GISSP150')
MMM_win = MMM_win.drop_vars('GISSP151')
MMM_win = MMM_win.drop_vars('MPIP1')
MMM_win = MMM_win.drop_vars('MPIP2')
MMM_LOVE1_win_3 = LOVE1_win_3*100
MMM_LOVE2_win_3 = LOVE2_win_3*100

# multi model mean
MMM_win_mean = (MMM_win.CNRM + MMM_win.FGOALS + MMM_win.IPSL + MMM_win.MIROC + MMM_win.MRI +  MMM_win.CCSM4 + MMM_win.MPI + MMM_win.GISS)/8

# calculating STD
var1=(MMM_win.CCSM4-MMM_win_mean)**2
var2=(MMM_win.GISS-MMM_win_mean)**2
var3=(MMM_win.MPI-MMM_win_mean)**2
var4=(MMM_win.CNRM-MMM_win_mean)**2
var5=(MMM_win.FGOALS-MMM_win_mean)**2
var6=(MMM_win.IPSL-MMM_win_mean)**2
var7=(MMM_win.MIROC-MMM_win_mean)**2
var8=(MMM_win.MRI-MMM_win_mean)**2
STD=((var1+var2+var3+var4+var5+var6+var7+var8)/8)**0.5
STD4=(MMM_win_mean+STD)
STD5=(MMM_win_mean-STD)

MMM_winter_PMIP3 = MMM_win_mean
STD_MMM_winter_PMIP3 = STD4
MMM_winter_PMIP3_STD = STD5

####### LOVECLIM sensitivity  #####

#austral summer

MMM = (MMM_LOVE2_sum_3 + MMM_LOVE1_sum_3)/2
var1 = (MMM_LOVE1_sum_3 - MMM)**2
var2 = (MMM_LOVE2_sum_3 - MMM)**2
STD = ((var1+var2)/2)**0.5
STD1 = (MMM+STD)
STD2 = (MMM-STD)

MMM_summer_LOVE = MMM
STD_MMM_summer_LOVE = STD1
MMM_summer_LOVE_STD = STD2

MMM = (MMM_LOVE1_win_3 + MMM_LOVE2_win_3 )/2
var1 = (MMM_LOVE1_win_3 - MMM)**2
var2 = (MMM_LOVE2_win_3 - MMM)**2
STD = ((var1+var2)/2)**0.5
STD1 = (MMM+STD)
STD2 = (MMM-STD)

MMM_winter_LOVE = MMM
STD_MMM_winter_LOVE = STD1
MMM_winter_LOVE_STD = STD2

#### PMIP4 ########

MMM_summer_mean = (CCSM4UoT_summer.siconc+CESM12_summer.IFRAC*100+LOVEsummer+AWI_summer.siconca+MPI_summer.siconc+MIROC_summer.siconc+IPSL_summer.fract_sic*100)/7
var1=(CCSM4UoT_summer.siconc-MMM_summer_mean)**2
var2=(CESM12_summer.IFRAC*100-MMM_summer_mean)**2
var3=(AWI_summer.siconca-MMM_summer_mean)**2
var4=(MPI_summer.siconc-MMM_summer_mean)**2
var5=(MIROC_summer.siconc-MMM_summer_mean)**2
var6=(IPSL_summer.fract_sic*100-MMM_summer_mean)**2
var7=(LOVEsummer-MMM_summer_mean)**2
STD=((var1+var2+var3+var4+var5+var6+var7)/7)**0.5
STD2=(MMM_summer_mean+STD)
STD3=(MMM_summer_mean-STD)

MMM_summer_PMIP4 = MMM_summer_mean
STD_MMM_summer_PMIP4 = STD2
MMM_summer_PMIP4_STD = STD3

MMM_winter_mean = (CCSM4UoT_winter.siconc+CESM12_winter.IFRAC*100+LOVEwinter+AWI_winter.siconca+MPI_winter.siconc+MIROC_winter.siconc+IPSL_winter.fract_sic*100)/7
var1=(CCSM4UoT_winter.siconc-MMM_winter_mean)**2
var2=(CESM12_winter.IFRAC*100-MMM_winter_mean)**2
var3=(AWI_winter.siconca-MMM_winter_mean)**2
var4=(MPI_winter.siconc-MMM_winter_mean)**2
var5=(MIROC_winter.siconc-MMM_winter_mean)**2
var6=(IPSL_winter.fract_sic*100-MMM_winter_mean)**2
var7=(LOVEwinter-MMM_winter_mean)**2
STD=((var1+var2+var3+var4+var5+var6+var7)/7)**0.5
STD2=(MMM_winter_mean+STD)
STD3=(MMM_winter_mean-STD)

MMM_winter_PMIP4 = MMM_winter_mean
STD_MMM_winter_PMIP4 = STD2
MMM_winter_PMIP4_STD = STD3

# proxy data points for Summer
noSSIlat = [-57.02,-63.08,-60.22,-61.24,-56.9,-57.83,-56.06,-54.91,-62.16,-58.68,-59.97,-59.02,-45.1,-53.93,-52.6,-53.55,-55.55,-55.14,-56.74,-59.79,-59.79,-60.3,-59.05,-58.72,-54.64,-47.9,-50.25,-54.91,-56.84,-55.97,-59.23,-55.19,-57.57,-58.85,-53.05,-46.51,-51.7,-44.15,-50.4,-49.01,-48.23,-37.27,-50.95,-42.88,-53.07,-53.88,-43.17,-53.18,-52.63,-54.91,-53.63,-52.59,-51.98,-51.83,-53.66,-53.17,-53.18,-53.18,-53.19,-50.16,-49.98,-44.96,-46.94,-48.9,-46.6,-42.87,-55.37,-42.52,-43.22,-47.54,-52.02,-45.73,-43.7,-49.1,-43.85,-41.86,-49.57,-56.57,-44.56,-46.14,-58.99,-50.69,-53.23,-51.88,-42.88,-51.91,-55.01,-43.7,-49,-38,-43.49,-54.48,-44.98,-52.2,-38.75,-50.47,-50.32,-48.89,-55.33,-51.06,-50.37,-37.8,-55,-54.91,-38.54,-46.06,-46.02,-55.95,-44.88,-52.9,-51.01,-55.07,-53.04,-53.5,-52.48,-53.44,-45.06,-52.94,-47.77,-45.01,-48.86,-48.03,-54.19,-56.38,-60.39,-56.67]
noSSIlon = [-160.09,-135.12,-127.03,-116.05,-115.24,-115.21,-115.06,-114.7,-109.09,-108.8,-101.32,-99.76,-57.95,-48.04,-46.88,-45.29,-45.02,-44.11,-42.97,-42.68,-39.6,-36.65,-35.61,-33.04,-23.95,-23.7,-23.24,-22.71,-22.32,-22.22,-19.73,-18.61,-17.1,-16.65,-16.45,-15.33,-15.3,-14.23,-14.08,-12.7,-11.04,-10.1,-7.51,-6.02,-4.99,-4.93,-4.06,-0.35,-0.13,3.31,3.86,4.48,4.52,4.81,5.1,5.11,5.13,5.13,5.33,5.72,5.87,5.97,6.26,6.71,7.63,8.92,9.98,11.67,11.74,15.36,20.47,25.65,25.73,27.38,27.6,28.54,30.02,34.18,34.79,35.9,37.63,40.13,40.8,41.65,42.35,42.88,45.01,45.06,45.21,51.18,51.32,53.05,53.28,54.47,59.3,59.58,61.2,61.66,65.47,67.72,68.39,71.53,73.26,73.84,79.87,90.09,96.43,104.95,106.52,109.85,109.99,110.02,110.05,111.33,114.09,114.26,114.37,116.99,123.1,125.98,126.02,126.13,144.79,145.3,157.53,160.23]
SSIlat = [-50.87,-53.80,-54.09,-62.49,-63.65,-64.67]
SSIlon = [-9.87,-8.22,-0.35,95.89,101.15,119.51]
maybeSSIlat = [-59.79,-60.30,-53.88,-53.18,-53.63,-53.17,-52.02]
maybeSSIlon = [-39.60,-36.65,-4.93,-0.35, 3.86, 5.11,20.47]

# proxy data points for winter
noWSIlat = [-58.55,-59.70,-60.87,-59.04,-55.53,-60.22,-54.22,-56.90,-56.06,-55.16,-54.91,-58.68,-52.81,-59.97,-59.02,-45.10,-56.84,-55.97,-46.51,-44.15,-50.40,-49.01,-48.23,-37.27,-42.88,-43.17,-50.16,-44.96,-46.94,-48.90,-46.60,-42.87,-42.52,-43.22,-47.54,-45.73,-43.70,-43.85,-41.86,-44.56,-46.14,-42.88,-43.70,-49.00,-38.00,-43.49,-44.98,-38.75,-48.89,-37.80,-38.54,-46.06,-46.02,-44.88,-45.06,-52.94,-47.77,-45.01,-48.86,-48.03]
noWSIlon = [-172.70,-171.36,-169.55,-158.36,-156.14,-127.03,-125.43,-115.24,-115.06,-114.79,-114.70,-108.80,-107.81,-101.32,-99.76,-57.95,-22.32,-22.22,-15.33,-14.23,-14.08,-12.70,-11.04,-10.10,-6.02,-4.06,5.72,5.97,6.26,6.71,7.63,8.92,11.67,11.74,15.36,25.65,25.73,27.60,28.54,34.79,35.90,42.35,45.06,45.21,51.18,51.32,53.28,59.30,61.66,71.53,79.87,90.09,96.43,106.52,114.37,116.99,123.10,125.98,126.02,126.13]
WSIlat = [-63.69,-61.94,-57.02,-57.20,-57.56,-61.01,-63.08,-62.03,-61.24,-57.83,-56.15,-59.21,-62.16,-53.93,-52.60,-53.55,-55.55,-55.14,-56.74,-59.79,-59.79,-60.30,-59.05,-58.72,-54.64,-47.90,-50.25,-54.91,-59.23,-55.19,-57.57,-58.85,-53.05,-51.70,-50.87,-53.80,-50.95,-53.07,-53.88,-54.09,-53.18,-52.63,-54.91,-53.63,-52.59,-51.98,-51.83,-53.66,-53.17,-53.18,-53.18,-53.19,-49.98,-55.37,-52.02,-49.10,-49.57,-56.57,-58.99,-50.69,-53.23,-51.88,-51.91,-55.01,-54.48,-52.20,-50.47,-50.32,-55.33,-51.06,-50.37,-55.00,-54.91,-62.49,-63.65,-55.95,-52.90,-51.01,-55.07,-53.04,-53.50,-52.48,-53.44,-64.67,-54.19,-56.38,-59.62,-60.39,-56.67]
WSIlon = [-169.07,-160.12,-160.09,-151.61,-151.22,-139.46,-135.12,-116.12,-116.05,-115.21,-115.13,-114.89,-109.09,-48.04,-46.88,-45.29,-45.02,-44.11,-42.97,-42.68,-39.60,-36.65,-35.61,-33.04,-23.95,-23.70,-23.24,-22.71,-19.73,-18.61,-17.10,-16.65,-16.45,-15.30,-9.87,-8.22,-7.51,-4.99,-4.93,-0.35,-0.35,-0.13,3.31,3.86,4.48,4.52,4.81,5.10,5.11,5.13,5.13,5.33,5.87,9.98,20.47,27.38,30.02,34.18,37.63,40.13,40.80,41.65,42.88,45.01,53.05,54.47,59.58,61.20,65.47,67.72,68.39,73.26,73.84,95.89,101.15,104.95,109.85,109.99,110.02,110.05,111.33,114.09,114.26,119.51,144.79,145.30,155.24,157.53,160.23]
maybeWSIlat = [-59.04,-55.53,-56.84,-55.97,-50.40,-50.16,-47.45,-43.70,-49.00]
maybeWSIlon = [-158.36,-156.14,-22.32,-22.22,-14.08,5.72,15.36,45.06,45.21]

# proxy line from
proxyline = pd.read_excel('~/Downloads/Table_LGM_SI_Zenedo (1).xlsx',engine="openpyxl")
WSI = proxyline[["Longitude", "Latitude WSI"]]
SSI = proxyline[["Longitude", "Latitude SSI"]]
WSI = WSI[1:361]
SSI = SSI[1:361]
WSI = WSI.rename(columns={"Latitude WSI": "Latitude"})
SSI = SSI.rename(columns={"Latitude SSI": "Latitude"})



# used to count for sea ice agreement
# if point lands on line then count as inside

extent = [-180, 180, -90, -35]
fig = plt.figure(figsize=[10, 5])

ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.coastlines()
ax.add_feature(cfeature.LAND, facecolor='0.8', edgecolor='k',zorder=0)



extent = [-180, 180, -90, -35]
ax.set_extent(extent, crs=ccrs.PlateCarree())
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.gridlines(linestyle='--',zorder=6)

#test = 'summer'
test = 'winter'

# no contains the maybe points as well ? I am pretty sure
if test == 'summer':
    ax.plot([noSSIlon], [noSSIlat],color='red', marker='o',markersize=2,transform=ccrs.PlateCarree())
    ax.plot([SSIlon], [SSIlat],color='blue', marker='o',markersize=2,transform=ccrs.PlateCarree())

if test == 'winter':
    ax.plot([noWSIlon], [noWSIlat],color='red', marker='o',markersize=1.5,transform=ccrs.PlateCarree())
    ax.plot([WSIlon], [WSIlat],color='blue', marker='o',markersize=1.5,transform=ccrs.PlateCarree())

ax.plot([maybeSSIlon], [maybeSSIlat],color='black', marker='o',markersize=1.5,transform=ccrs.PlateCarree())

ax.contour(MMM_winter_PMIP4.lon,MMM_winter_PMIP4.lat,MMM_winter_PMIP4,colors=['black'],levels=[15], transform=ccrs.PlateCarree(),linewidths=3)
#ax.contour(LOVE2_win_3.LON,LOVE2_win_3.LAT,LOVE2_win_3,colors=['black'],levels=[0.15],linestyles=['solid'], transform=ccrs.PlateCarree(),linewidths=3)
#ax[9].contour(LOVE2_win_3.LON,LOVE2_win_3.LAT,LOVE2_win_3,colors=['#F7DC6F'],levels=[0.15],linestyles=['dotted'], transform=ccrs.PlateCarree(),linewidths=1.5)
#ax.contour(MMM_sum_mean.LON,MMM_sum_mean.LAT,MMM_sum_mean,colors=['black'],levels=[15], transform=ccrs.PlateCarree(),linewidths=3)

plt.savefig('Figures/counting/MMM_PMIP4_{}_withmaybe.pdf'.format(test))
