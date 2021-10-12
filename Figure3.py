# Trying to rewrite figure 3 in Python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import xarray as xr
import matplotlib as mpl
import cartopy.crs as ccrs
import cmocean
import cartopy
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

SST = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/MIROC/SST/tos_Omon_MIROC-ES2L_lgm_r1i1p1f2_gn_320001-329912.nc',decode_times=False)
SST = SST.mean(dim='time')

proxy = pd.read_excel('~/Desktop/UNSW/Table_Recap_LGM_Final.xlsx',engine="openpyxl",sheet_name='Data')
data = proxy[['Latitude','Longitude', 'LGM']]
data = data[3:]
data = data.dropna()
data = data.reset_index(drop=True)
# lat is -90:90 and lon is -180:180
lat = data['Latitude']
lon = data['Longitude']
sst = data['LGM']
data['LGM'] = data['LGM'].astype(float)
def LonTo360(dlon):
    # Convert longitudes to 0-360 deg
    dlon = ((360 + (dlon % 360)) % 360)
    return dlon
lonconverted = LonTo360(lon)
#Read the data and convert into numpy array
y = np.array(data.Latitude)
x = np.array(data.Longitude)
z = np.array(data.LGM)

# PMIP3 SST
CNRM_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/tos_Oclim_CNRM-CM5_lgm_r1i1p1_180001-199912-climregrid.nc',decode_times=False)
FGOALS_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/tos_Oclim_FGOALS-g2_lgm_r1i1p1_055001-064912-climregrid.nc',decode_times=False)
IPSL_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/tos_Oclim_IPSL-CM5A-LR_lgm_r1i1p1_260101-280012-climregrid.nc',decode_times=False)
MIROC_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/tos_Oclim_MIROC-ESM_lgm_r1i1p1_460001-469912-climregrid.nc',decode_times=False)
MRI_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/MRISST_winter_PMIP3.nc',decode_times=False)
CCSM4_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/CCSM4SST_winter_PMIP3.nc',decode_times=False)
MPI_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/MPISST_winter_PMIP3.nc',decode_times=False)
GISS_3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP3SST/GISSSST_winter_PMIP3.nc',decode_times=False)
CNRM_3 = xr.concat([CNRM_3.sel(time=(CNRM_3.time[1])), CNRM_3.sel(time=(CNRM_3.time[2]))], dim="time")
CNRM_3 = CNRM_3.mean(dim='time')
FGOALS_3= xr.concat([FGOALS_3.sel(time=(FGOALS_3.time[2])), FGOALS_3.sel(time=(FGOALS_3.time[3]))], dim="time")
FGOALS_3 = FGOALS_3.mean(dim='time')
IPSL_3 = xr.concat([IPSL_3.sel(time=(IPSL_3.time[1])), IPSL_3.sel(time=(IPSL_3.time[2]))], dim="time")
IPSL_3 = IPSL_3.mean(dim='time')
MIROC_3 = xr.concat([MIROC_3.sel(time=(MIROC_3.time[1])), MIROC_3.sel(time=(MIROC_3.time[2]))], dim="time")
MIROC_3 = MIROC_3.mean(dim='time')
MRI_3 = MRI_3.rename({"LON": "lon","LAT": "lat","TIME8": "time","TOSMRI":"tos"})
MRI_3 = xr.concat([MRI_3.sel(time=(MRI_3.time[1])), MRI_3.sel(time=(MRI_3.time[2]))], dim="time")
MRI_3 = MRI_3.mean(dim='time')
CCSM4_3 = CCSM4_3.rename({"LON": "lon","LAT": "lat","CCSM4":"tos"})
GISS_3 = GISS_3.rename({"LON": "lon","LAT": "lat","GISS":"tos"})
MPI_3 = MPI_3.rename({"LON": "lon","LAT": "lat","MPI":"tos"})
CNRM_3 -= 273.15
GISS_3 -= 273.15
IPSL_3 -= 273.15
MIROC_3 -= 273.15
MPI_3 -= 273.15
MRI_3 -= 273.15
FGOALS_3 -= 273.15
CCSM4_3 -= 273.15
PMIP3names = ['CNRM','GISS-E2-R','IPSL-CM5A-LR','MIROC-ESM-P','MPI-ESM-P','MRI-CGCM3','FGOALS-G2','CCSM4']
PMIP3models = [CNRM_3,GISS_3,IPSL_3,MIROC_3,MPI_3,MRI_3,FGOALS_3,CCSM4_3]

# PMIP4 SST
AWI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/AWI_SST_lonconverted.nc',decode_times=False)
CESM12 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/b.e12.B1850C5.f19_g16.i21ka.03.pop.h.vars.08010900.climo_regrid.nc',decode_times=False)
LOVECLIM = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4seaice/cresum18250_regrid.nc',decode_times=False)
MIROC_4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/MIROC_PMIP4_tos_regrid.nc',decode_times=False)
IPSL_4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/IPSL_SSTfixed.nc',decode_times=False)
IPSL_4 = IPSL_4.rename({"LON": "lon","LAT": "lat","TIME_COUNTER":"time","SA":"tos"})
CCSM4UoT = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/CCSM4-UoT_sst.nc',decode_times=False)
MPI_4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/MPI_PMIP4_tos_regrid.nc',decode_times=False)

# getting summer months
AWI = xr.concat([AWI.sel(time=(AWI.time[1::12])), AWI.sel(time=(AWI.time[2::12]))], dim="time") # Feb March
AWI = AWI.mean(dim='time')
AWI = AWI.tos
MIROC_4 = xr.concat([MIROC_4.sel(time=(MIROC_4.time[1::12])), MIROC_4.sel(time=(MIROC_4.time[2::12]))], dim="time") # Feb march
MIROC_4 = MIROC_4.mean(dim='time')
MIROC_4 = MIROC_4.tos
CESM12 = xr.concat([CESM12.sel(time=(CESM12.time[1])), CESM12.sel(time=(CESM12.time[2]))], dim="time")
CESM12 = CESM12.mean(dim='time')
CESM12 = CESM12.TEMP
CESM12 = CESM12[0,:,:]
CCSM4UoT = xr.concat([CCSM4UoT.sel(time=(CCSM4UoT.time[2::12])), CCSM4UoT.sel(time=(CCSM4UoT.time[3::12]))], dim="time")
CCSM4UoT = CCSM4UoT.mean(dim='time')
CCSM4UoT = CCSM4UoT.tos
IPSL_4 = xr.concat([IPSL_4.sel(time=(IPSL_4.time[1])), IPSL_4.sel(time=(IPSL_4.time[2]))], dim="time") # Feb March
IPSL_4 = IPSL_4.mean(dim='time')
IPSL_4 = IPSL_4.tos-273.15
MPI_4 = xr.concat([MPI_4.sel(time=(MPI_4.time[1::12])), MPI_4.sel(time=(MPI_4.time[2::12]))], dim="time")
MPI_4 = MPI_4.mean(dim='time')
MPI_4 = MPI_4.tos
LOVEsic = LOVECLIM.sst
feb = LOVEsic[1:2400:12]
mar = LOVEsic[2:2400:12]
mar = mar.mean(dim='time')
feb = feb.mean(dim='time')
LOVECLIM = (mar+feb)/2
PMIP4models = [MIROC_4, IPSL_4, MPI_4,AWI,LOVECLIM,CESM12,CCSM4UoT]
PMIP4names= ['MIROC-ES2L','IPSL-CM5A2','MPI-ESM1-2','AWI-ESM-1','LOVECLIM','CESM1.2','UoT-CCSM4']


# LOVECLIM sensitivity models SST
weakNA = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/LOVESST/V3LNAw-SSTREGRID.nc')
weakNA_AB = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/LOVESST/weakNA_ABSST_sum.nc')
weakNA = weakNA.rename({"TIME":"time","SST":"tos"})
weakNA_AB = weakNA_AB.rename({"AX006":"time","SSTNHWSOW":"tos","LON":"lon","LAT":"lat"})

weakNA = weakNA.tos
feb = weakNA[1::12,:,:]
feb = feb.mean(dim='time')
mar = weakNA[2::12,:,:]
mar = mar.mean(dim='time')
weakNA_sum = (feb+mar)/2
weakNA_AB = weakNA_AB.tos
feb1 = weakNA_AB[1::12,:,:]
feb1 = feb1.mean(dim='time')
mar1 = weakNA_AB[2::12,:,:]
mar1 = mar1.mean(dim='time')
weakNA_AB_sum = (feb1+mar1)/2
LOVEmodels = [weakNA_sum,weakNA_AB_sum]
LOVEnames = ['weakNA','weakNA_AB']

# legends

proxy = mlines.Line2D([], [], color='black', linestyle ='solid',label = 'Proxy 1˚ isoline')
model = mlines.Line2D([], [], color='black', linestyle ='dashed',label = 'Model 1˚ isoline')


nrows=5
ncols=4
extent = [-180, 180, -90, -35]

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,
                        subplot_kw={'projection': ccrs.SouthPolarStereo()},
                        figsize=(16,14))
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.weight"] = "bold"
ax=ax.flatten()
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

cm = plt.cm.get_cmap('RdYlBu_r')
cmap = mpl.cm.RdYlBu_r
bounds = [-2,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,10,12]

norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')


# plt.
# 60 has been good
# 21 levels brought me down to 3.5 MB
levels= 21

for n in range(15):
    ax[n].set_boundary(circle, transform=ax[n].transAxes)
    ax[n].gridlines(linestyle='--',zorder=4)
    ax[n].set_extent(extent, crs=ccrs.PlateCarree())
    #gl = ax[n].gridlines(crs=ccrs.PlateCarree(),draw_labels=True, dms=True, x_inline=True, y_inline=True,linestyle='--')
    ax[n].coastlines(zorder=4)
    ax[n].add_feature(cfeature.LAND, facecolor='1', edgecolor='k',zorder=3)
    ax[n].set_aspect('equal')
for n in range(16,18):
    ax[n].set_boundary(circle, transform=ax[n].transAxes)
    ax[n].gridlines(linestyle='--',zorder=4)
    ax[n].set_extent(extent, crs=ccrs.PlateCarree())
    #gl = ax[n].gridlines(crs=ccrs.PlateCarree(),draw_labels=True, dms=True, x_inline=True, y_inline=True,linestyle='--')
    ax[n].coastlines(zorder=4)
    ax[n].add_feature(cfeature.LAND, facecolor='1', edgecolor='k',zorder=3)
    ax[n].set_aspect('equal')

for i in range(len(PMIP3models)):
    ax[i].contourf(PMIP3models[i].lon,PMIP3models[i].lat,PMIP3models[i].tos,norm=norm,cmap=cm,levels = levels, transform=ccrs.PlateCarree())
    ax[i].contour(PMIP3models[i].lon,PMIP3models[i].lat,PMIP3models[i].tos,colors='black',linestyles='dashed',levels = [1], transform=ccrs.PlateCarree())
    ax[i].scatter(data.Longitude,data.Latitude,c=data.LGM, norm=norm,s=27,cmap=cm,edgecolors='k', transform=ccrs.PlateCarree(),zorder=5)
    one = ax[i].tricontour(lonconverted,y,z,colors=['black'],levels=[1], transform=ccrs.PlateCarree(),zorder=4)
    ax[i].text(160,-28,PMIP3names[i],transform=ccrs.PlateCarree())

for j in range(7):
    ax[j+8].contourf(PMIP4models[j].lon,PMIP4models[j].lat,PMIP4models[j],norm=norm,cmap=cm,levels = levels, transform=ccrs.PlateCarree())
    ax[j+8].contour(PMIP4models[j].lon,PMIP4models[j].lat,PMIP4models[j],colors='black',linestyles='dashed',levels = [1], transform=ccrs.PlateCarree())
    ax[j+8].scatter(data.Longitude,data.Latitude,c=data.LGM, norm=norm,s=27,cmap=cm,edgecolors='k', transform=ccrs.PlateCarree(),zorder=5)
    one = ax[j+8].tricontour(lonconverted,y,z,colors=['black'],levels=[1], transform=ccrs.PlateCarree(),zorder=4)
    ax[j+8].text(160,-28,PMIP4names[j],transform=ccrs.PlateCarree())

for l in range(2):
    ax[l+16].contourf(LOVEmodels[l].lon,LOVEmodels[l].lat,LOVEmodels[l],norm=norm,cmap=cm,levels = levels, transform=ccrs.PlateCarree())
    ax[l+16].contour(LOVEmodels[l].lon,LOVEmodels[l].lat,LOVEmodels[l],colors='black',linestyles='dashed',levels = [1],transform=ccrs.PlateCarree())
    ax[l+16].scatter(data.Longitude,data.Latitude,c=data.LGM, norm=norm,s=27,cmap=cm,edgecolors='k', transform=ccrs.PlateCarree(),zorder=5)
    one = ax[l+16].tricontour(lonconverted,y,z,colors=['black'],levels=[1], transform=ccrs.PlateCarree(),zorder=4)
    ax[l+16].text(160,-28,LOVEnames[l],transform=ccrs.PlateCarree())

# plt.subplots_adjust(hspace=0.07,wspace=0)

# cb_ax =
# cbar =
# cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),cax=cb_ax)
plt.tight_layout()
cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax.ravel().tolist(), shrink=0.95)
cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
cb.ax.tick_params(labelsize='large')
plt.legend(handles=[proxy,model],frameon=False,fontsize=15)
ax[15].axis('off')
ax[18].axis('off')
ax[19].axis('off')
plt.show()
#plt.savefig('Figures/Figure3_9.27.21.eps',dpi=200)
