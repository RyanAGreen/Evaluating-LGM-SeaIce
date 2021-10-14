# two columns PMIP
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

MMM_summer_PMIP3 = MMM_summer_PMIP3.to_dataset(name='sic')
MMM_summer_PMIP3.to_netcdf("MMM_STD/PMIP3_summer_MMM.nc")

STD_MMM_summer_PMIP3 = STD_MMM_summer_PMIP3.to_dataset(name='sic')
STD_MMM_summer_PMIP3.to_netcdf("MMM_STD/PMIP3_summer_STD_MMM.nc")

MMM_summer_PMIP3_STD = MMM_summer_PMIP3_STD.to_dataset(name='sic')
MMM_summer_PMIP3_STD.to_netcdf("MMM_STD/PMIP3_summer_MMM_STD.nc")

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

MMM_winter_PMIP3 = MMM_winter_PMIP3.to_dataset(name='sic')
MMM_winter_PMIP3.to_netcdf("MMM_STD/PMIP3_winter_MMM.nc")

STD_MMM_winter_PMIP3 = STD_MMM_winter_PMIP3.to_dataset(name='sic')
STD_MMM_winter_PMIP3.to_netcdf("MMM_STD/PMIP3_winter_STD_MMM.nc")

MMM_winter_PMIP3_STD = MMM_winter_PMIP3_STD.to_dataset(name='sic')
MMM_winter_PMIP3_STD.to_netcdf("MMM_STD/PMIP3_winter_MMM_STD.nc")

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

MMM_summer_LOVE = MMM_summer_LOVE.to_dataset(name='sic')
MMM_summer_LOVE.to_netcdf("MMM_STD/LOVE_summer_MMM.nc")

STD_MMM_summer_LOVE = STD_MMM_summer_LOVE.to_dataset(name='sic')
STD_MMM_summer_LOVE.to_netcdf("MMM_STD/LOVE_summer_STD_MMM.nc")

MMM_summer_LOVE_STD = MMM_summer_LOVE_STD.to_dataset(name='sic')
MMM_summer_LOVE_STD.to_netcdf("MMM_STD/LOVE_summer_MMM_STD.nc")

MMM = (MMM_LOVE1_win_3 + MMM_LOVE2_win_3 )/2
var1 = (MMM_LOVE1_win_3 - MMM)**2
var2 = (MMM_LOVE2_win_3 - MMM)**2
STD = ((var1+var2)/2)**0.5
STD1 = (MMM+STD)
STD2 = (MMM-STD)

MMM_winter_LOVE = MMM
STD_MMM_winter_LOVE = STD1
MMM_winter_LOVE_STD = STD2

MMM_winter_LOVE = MMM_winter_LOVE.to_dataset(name='sic')
MMM_winter_LOVE.to_netcdf("MMM_STD/LOVE_winter_MMM.nc")

STD_MMM_winter_LOVE = STD_MMM_winter_LOVE.to_dataset(name='sic')
STD_MMM_winter_LOVE.to_netcdf("MMM_STD/LOVE_winter_STD_MMM.nc")

MMM_winter_LOVE_STD = MMM_winter_LOVE_STD.to_dataset(name='sic')
MMM_winter_LOVE_STD.to_netcdf("MMM_STD/LOVE_winter_MMM_STD.nc")

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

MMM_summer_PMIP4 = MMM_summer_PMIP4.to_dataset(name='sic')
MMM_summer_PMIP4.to_netcdf("MMM_STD/PMIP4_summer_MMM.nc")

STD_MMM_summer_PMIP4 = STD_MMM_summer_PMIP4.to_dataset(name='sic')
STD_MMM_summer_PMIP4.to_netcdf("MMM_STD/PMIP4_summer_STD_MMM.nc")

MMM_summer_PMIP4_STD = MMM_summer_PMIP4_STD.to_dataset(name='sic')
MMM_summer_PMIP4_STD.to_netcdf("MMM_STD/PMIP4_summer_MMM_STD.nc")


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

MMM_winter_PMIP4 = MMM_winter_PMIP4.to_dataset(name='sic')
MMM_winter_PMIP4.to_netcdf("MMM_STD/PMIP4_winter_MMM.nc")

STD_MMM_winter_PMIP4 = STD_MMM_winter_PMIP4.to_dataset(name='sic')
STD_MMM_winter_PMIP4.to_netcdf("MMM_STD/PMIP4_winter_STD_MMM.nc")

MMM_winter_PMIP4_STD = MMM_winter_PMIP4_STD.to_dataset(name='sic')
MMM_winter_PMIP4_STD.to_netcdf("MMM_STD/PMIP4_winter_MMM_STD.nc")
