import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# PMIP3
CCSM4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/CCSM4ptemp_PMIP3.nc')
CNRM = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/CNRMptemp_PMIP3.nc')
FGOALS = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/FGOALSptemp_PMIP3.nc')
GISS = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/GISSptemp_PMIP3.nc')
IPSL = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/IPSLptemp_PMIP3.nc')
MIROC = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/MIROCptemp_PMIP3.nc')
MPI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/MPIptemp_PMIP3.nc')
MRI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/MRIptemp_PMIP3.nc')

# PMIP4
CESM12 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/PMIP4SST/b.e12.B1850C5.f19_g16.i21ka.03.pop.h.vars.08010900.climo_regrid.nc',decode_times=False)
AWI = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/AWIptemp_PMIP4.nc',decode_times=False)
MPI4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/MPIptemp_PMIP4.nc',decode_times=False)
LOVE = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/LOVECLIM_ptemp_ATL.nc',decode_times=False) # annual mean instead of summer
MIROC_4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/MIROC_ptemp_PMIP4_atl.nc',decode_times=False)
# UoTCCSM4 only sst data
# IPSL only sst data


# LOVECLIM sensitivies
weakNA = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/weakNA_ABptemp_PMIP3.nc',decode_times=False)
weakNA_AB = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/ptemp/weakNAptemp_PMIP3.nc',decode_times=False)

# cleaning up PMIP3
CCSM4 = CCSM4.PTEMP
CNRM = CNRM.PTEMP
FGOALS = FGOALS.PTEMP
GISS = GISS.PTEMP
IPSL = IPSL.PTEMP
MIROC = MIROC.PTEMP
MPI = MPI.PTEMP
MRI = MRI.PTEMP

# cleaning up PMIP4
CESM12 = xr.concat([CESM12.sel(time=(CESM12.time[1])), CESM12.sel(time=(CESM12.time[2]))], dim="time")
CESM12 = CESM12.mean(dim='time')
CESM12 = CESM12.TEMP
CESMatl = xr.concat([CESM12.sel(lon=(CESM12.lon[300:360])), CESM12.sel(lon=(CESM12.lon[0:20]))], dim="lon")
CESMatl = CESMatl.mean(dim='lon')
CESMatl['z_t'] = CESMatl['z_t']/100
CESM12 = CESMatl
LOVE = LOVE.PT
LOVE = LOVE.mean(dim='TIME')
AWI = AWI.PTEMP
AWI = AWI.mean(dim='AX007')
MPI4 = MPI4.PTEMP
MPI4 = MPI4.mean(dim='AX007')
MIROC_4 = MIROC_4.PTEMP
# MIROC_4 = MIROC_4.mean(dim='time')

# cleaning up LOVE
weakNA = weakNA.PTEMP
weakNA_AB = weakNA_AB.PTEMP


# putting into lists
PMIP3models = [CNRM,GISS,IPSL,MIROC,MPI,MRI,FGOALS,CCSM4]
PMIP3names = ['CNRM','GISS-E2-R','IPSL-CM5A-LR','MIROC-ESM-P','MPI-ESM-P','MRI-CGCM3','FGOALS-G2','CCSM4']

#PMIP4models = [MIROC_4, IPSL_4, MPI_4,AWI,LOVECLIM,CESM12,CCSM4UoT]
#PMIP4names= ['MIROC-ES2L','IPSL-CM5A2','MPI-ESM1-2','AWI-ESM-1','LOVECLIM','CESM1.2','UoT-CCSM4']
PMIP4names= ['MIROC-ES2L','MPI-ESM1-2','AWI-ESM-1','LOVECLIM','CESM1.2']

LOVEmodels = [weakNA,weakNA_AB]
LOVEnames = ['weakNA','weakNA_AB']

# renaming
PMIP3models[0] = PMIP3models[0].rename({"LEV1": "lev","LAT": "lat"})
PMIP3models[1] = PMIP3models[1].rename({"LEV3": "lev","LAT": "lat"})
PMIP3models[2] = PMIP3models[2].rename({"LEV4": "lev","LAT": "lat"})
PMIP3models[3] = PMIP3models[3].rename({"LEV5": "lev","LAT": "lat"})
PMIP3models[4] = PMIP3models[4].rename({"LEV6": "lev","LAT": "lat"})
PMIP3models[5] = PMIP3models[5].rename({"LEV7": "lev","LAT": "lat"})
PMIP3models[6] = PMIP3models[6].rename({"LEV2": "lev","LAT": "lat"})
PMIP3models[7] = PMIP3models[7].rename({"LEV": "lev","LAT": "lat"})

# setting up desired levels
first = np.arange(-2,2.2,0.2)
second = np.arange(2,8.5,.5)
third = np.arange(8,30,3)
bounds = []
for numbers in first:
    bounds.append(numbers)
for numbers in second:
    bounds.append(numbers)
for numbers in third:
    bounds.append(numbers)
nrows=4
ncols=4

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(16,10),sharey='row',sharex='col')

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.weight"] = "bold"
ax=ax.flatten()

cm = plt.cm.get_cmap('RdYlBu_r')
cmap = mpl.cm.RdYlBu_r

norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
levels = 75 #75

for n in range(15):
    for axis in ['top','bottom','left','right']:
        ax[n].spines[axis].set_linewidth(3)

 # PMIP 3
for i in range(len(PMIP3names)):
    ax[i].contourf(PMIP3models[i].lat,PMIP3models[i].lev,PMIP3models[i],norm=norm,cmap=cm,levels = levels)
    ax[i].text(15,4800,PMIP3names[i])
    ax[i].grid(ls=':')
    #ax[i].invert_yaxis()
    ax[i].set_xlim(-70,65)
    ax[i].set_ylim(5000,0) #by 500 is how i did it before
    ax[i].tick_params(axis="both", direction="out", length=5, width=3, color="black")


# until I get all the models going
ax[8].contourf(MIROC_4.LAT,MIROC_4.LEV,MIROC_4,norm=norm,cmap=cm,levels = levels)
ax[9].contourf(MPI4.LAT,MPI4.LEV,MPI4,norm=norm,cmap=cm,levels = levels)
ax[10].contourf(AWI.LAT,AWI.DEPTH,AWI,norm=norm,cmap=cm,levels = levels)
ax[11].contourf(LOVE.LAT,LOVE.Z,LOVE-273.15,norm=norm,cmap=cm,levels = 1000)
ax[12].contourf(CESM12.lat,CESM12.z_t,CESM12,norm=norm,cmap=cm,levels = levels)

# PMIP4
for j in range(5):
   # ax[j+8].contourf(PMIP4models[j].lon,PMIP4models[j].lat,PMIP4models[j],norm=norm,cmap=cm,levels = levels)
    ax[j+8].text(15,4800,PMIP4names[j])
    ax[j+8].grid(ls=':')
    #ax[i].invert_yaxis() CESM needs inversion
    ax[j+8].set_xlim(-70,65)
    ax[j+8].set_ylim(5000,0) #by 500 is how i did it before
    ax[j+8].tick_params(axis="both", direction="out", length=5, width=3, color="black")



# LOVE sensitivity
for l in range(2):
    ax[l+13].contourf(LOVEmodels[l].LAT,LOVEmodels[l].Z,LOVEmodels[l],norm=norm,cmap=cm,levels = levels)
    ax[l+13].text(15,4800,LOVEnames[l])
    ax[l+13].grid(ls=':')
    ax[l+13].invert_yaxis()
    ax[l+13].set_xlim(-70,65)
    ax[l+13].set_ylim(5000,0) #by 500 is how i did it before
    ax[l+13].tick_params(axis="both", direction="out", length=5, width=3, color="black")

ax[15].axis('off')
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel("Depth (m)",fontsize=15,fontweight='bold')
plt.xlabel("Latitude (˚)",fontsize=15,fontweight='bold')

plt.tight_layout()
cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax.ravel().tolist(), shrink=0.975)
cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
cb.ax.tick_params(labelsize='large')

plt.show()
#plt.savefig('Figures/Figure4_9.27.21.eps')
