import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import xarray as xr
import matplotlib as mpl

# most area to least area finished with MMM

# PMIP3
# CCSM4, FGOALS, MRI, MPI, MIROC, IPSL, GISS, CNRM  (MMM is commented after because removed later MMM)

seaiceedge_PMIP3 = [-55.5,-61.5,-62.5,-65,-66.5,-70,-65.5,-75.5] #,-62.5]
SSTtest_PMIP3 = [-0.5208,1.998,2.885,2.692,2.911,2.884,3.708,5.582] #,2.95]
seaicearea_PMIP3 =  [ 27.46,13.62,12.54,5.187,3.530,2.414,2.391,.06047] #,9.34]
colorsPMIP3 = ['#ff7f0e','#e377c2','#1f77b4','#bcbd22','#9467bd','#2ca02c','#17becf','black']

# PMIP4
#UoTCCSM4, CESM1.2, LOVECLIM, AWI, MPI, IPSL, MIROC (MMM is commented after because removed later MMM)
seaiceedge_PMIP4 = [-53,-57.5,-59.5,-62,-65,-70,-75.5] #,-59]
SSTtest_PMIP4 = [-1.017,0.161,2.218,0.957,2.451,3.404,5.583] #,1.774]
seaicearea_PMIP4 =  [ 33.15,23.75,17.55,14.73,5.18,2.46,0.36] #,19.08]
# need to decide color for MMM
colorsPMIP4 = ['#ff7f0e','#8c564b','#F7DC6F','#B8255F','#bcbd22','#2ca02c','#9467bd']


# LOVECLIM
# weakNA_AB (LOVE2), weak_AB (LOVE1) (MMM is commented after because removed later MMM)
seaiceedge_LOVE = [-58.5, -59.5] #,-59]
SSTtest_LOVE = [1.61,2.15] #, 1.875]
seaicearea_LOVE = [20.27, 15.73] #, 18.47]
colorsLOVE = ['#F7DC6F','#F7DC6F']


# proxy
proxyseaicearea = [15.9]
proxyseaiceedge = [-61]
proxySST = [1.52]
colorsproxy=['white']
facecolorsproxy=['black']

# total (excluding MMM and proxy)
# PMIP3, PMIP4, LOVECLIM PMIP3
#seaiceedge_total = np.array([-55.5,-61.5,-62.5,-65,-66.5,-70,-65.5,-75.5,-53,-57.5,-59.5,-62,-65,-70,-75.5])
seaiceedge_total = seaiceedge_PMIP3 + seaiceedge_PMIP4 + seaiceedge_LOVE
#SSTtest_total = [-0.5208,1.998,2.885,2.692,2.911,2.884,3.708,5.582,-1.017,0.161,2.218,0.957,2.451,3.404,5.583]
SSTtest_total = SSTtest_PMIP3 + SSTtest_PMIP4 + SSTtest_LOVE
#seaicearea_total = [27.46,13.62,12.54,5.187,3.530,2.414,2.391,.06047,33.15,23.75,17.55,14.73,5.18,2.46,0.36]
seaicearea_total = seaicearea_PMIP3 + seaicearea_PMIP4 + seaicearea_LOVE

# linear fit for sea ice edge vs SST
fit = np.polyfit(SSTtest_total,seaiceedge_total,1)
ang_coeff = fit[0]
intercept = fit[1]
fit_sie = ang_coeff*np.asarray(SSTtest_total) + intercept

# linear fit for sea ice area vs SST
fit = np.polyfit(SSTtest_total,seaicearea_total,1)
ang_coeff = fit[0]
intercept = fit[1]
fit_sia = ang_coeff*np.asarray(SSTtest_total) + intercept

# proxy data
proxy = pd.read_excel('~/Desktop/UNSW/Table_Recap_LGM_Final.xlsx',engine="openpyxl",sheet_name='Data')
data = proxy[['Latitude', 'LGM']]
data = data[3:]
data = data.dropna()
data = data.reset_index(drop=True)

MIROC4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/MIROCsstzonallyavg.nc')
IPSL4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/IPSLsstzonallyavg.nc')
MPI4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/MPIsstzonallyavg.nc')
AWI4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/AWIsstzonallyavg.nc')
LOVE4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/LOVEsstzonallyavg.nc')
CESM4 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/CESMsstzonallyavg.nc')
CCSM44 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/CCSM4sstzonallyavg.nc')
PMIP3 = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/TOSSummerAllModelsnew.nc')
LOVEsens = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/LOVE_summer_zonnalyaveSST.nc')

# PMIP3
CNRM3 = PMIP3.TOS2 - 273.15
GISS3 = PMIP3.TOS4 - 273.15
IPSL3 = PMIP3.TOS5 - 273.15
MIROC3 = PMIP3.TOS6 - 273.15
MPI3 = PMIP3.TOS7 - 273.15
MRI3 = PMIP3.TOS8
FGOALS3 = PMIP3.TOS3 - 273.15
CCSM43 = PMIP3.TOS1 - 273.15
PMIP3_zonal = [CCSM43,FGOALS3,MRI3,MPI3,MIROC3,IPSL3,GISS3,CNRM3]
PMIP3_zonal_colors = ['#ff7f0e','#e377c2','#1f77b4','#bcbd22','#9467bd','#2ca02c','#17becf','black']

# PMIP4
MIROC4 = MIROC4.TOS
IPSL4 = IPSL4.TOS
MPI4 = MPI4.TOS
AWI4 = AWI4.TOS
LOVE4 = LOVE4.TOS
CESM4 = CESM4.TOS
CCSM44 = CCSM44.TOS
PMIP4_zonal = [CCSM44,CESM4,LOVE4,AWI4,MPI4,IPSL4,MIROC4]
PMIP4_zonal_colors = ['#ff7f0e','#8c564b','#F7DC6F','#B8255F','#bcbd22','#2ca02c','#9467bd']

# LOVECLIM
weakNA = LOVEsens.WEAKNA_SUM_SST
weakNA_AB = LOVEsens.WEAKNA_AB_SUM_SST
LOVE_zonal = [weakNA,weakNA_AB]
LOVE_zonal_colors = ['#F7DC6F','#F7DC6F']
LOVE_linestyle = ['dotted','dashdot']

# set up x axis
latitude = np.arange(-89.5,90.5,1)

# setting up the legends

# legend will be least to greatest sea ice for all models
CNRM_legend = mpatches.Patch(color='black', label='CNRM')
MIROC_legend = mpatches.Patch(color='#9467bd', label='MIROC-ESM-P/MIROC-ES2L')
GISS_legend = mpatches.Patch(color='#17becf', label='GISS-E2-R')
IPSL_legend = mpatches.Patch(color='#2ca02c', label='IPSL-CM5A-LR/IPSL-CM5A2')
MPI_legend = mpatches.Patch(color='#bcbd22', label='MPI-ESM-P/MPI-ESM1-2')
MRI_legend = mpatches.Patch(color='#1f77b4', label='MRI-CGCM3')
FGOALS_legend = mpatches.Patch(color='#e377c2', label='FGOALS-G2')
AWI_legend = mpatches.Patch(color='#B8255F', label='AWI-ESM-1')
LOVE_legend = mpatches.Patch(color='#F7DC6F', label='LOVECLIM')
CESM_legend = mpatches.Patch(color='#8c564b', label='CESM1.2')
CCSM4_legend = mpatches.Patch(color='#ff7f0e', label='CCSM4/UoTCCSM4')
#MMM_legend = mpatches.Patch(color='black', label='Multi-model means')
proxy_legend = mpatches.Patch(color='grey', label='Proxy data')


PMIP3_legend = mlines.Line2D([], [], color='black', marker='^', linestyle='None',markersize=8,markerfacecolor='white', label='PMIP3 models')
PMIP4_legend = mlines.Line2D([], [], color='black', marker='s', linestyle='None',markersize=8,markerfacecolor='white', label='PMIP4 models')
weakNA_legend = mlines.Line2D([], [], color='black', marker='P', linestyle='None',markersize=8,markerfacecolor='white', label='LOVECLIM-weakNA')
weakNA_AB_legend = mlines.Line2D([], [], color='black', marker='X', linestyle='None',markersize=8,markerfacecolor='white', label='LOVECLIM-WeakNA_AB')
Proxy_legend = mlines.Line2D([], [], color='black', marker='o', linestyle='None',markersize=8,markerfacecolor='white', label='Proxy estimate')
Proxy_uncertainty = mlines.Line2D([], [], color='black', marker='o', linestyle='None',markersize=8,markerfacecolor='white',alpha=0.3, label='Proxy uncertainty')
#MMM = mlines.Line2D([], [], color='black', marker='o', linestyle='None',markersize=8, label='Multi-model mean')
densly = (0, (3, 1, 1, 1))
LOVE1_legend = mlines.Line2D([], [], color='#F7DC6F', linestyle ='dotted',label = 'LOVECLIM-weakNA')
LOVE2_legend = mlines.Line2D([], [], color='#F7DC6F', linestyle ='dashdot',label = 'LOVECLIM-weakNA_AB')
PMIP3_line = mlines.Line2D([], [], color='black', linestyle ='solid',label = 'PMIP3 model')
PMIP4_line = mlines.Line2D([], [], color='black', linestyle =densly,label = 'PMIP4 model')
degree = mlines.Line2D([], [], color='grey', linestyle ='dashed',label = '-2 ˚C')

uncertainty_sie = plt.Circle((proxySST, proxyseaiceedge), 0.67, color='k', fill=False)
uncertainty_sia = plt.Circle((proxySST, proxyseaicearea), 0.67, color='k', fill=False)

plt.figure(figsize=(16, 16))
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.weight"] = "bold"
mpl.rcParams['axes.linewidth'] = 3 #set the value globally

# plotting sea ice edge

plt.subplot(3, 2, 1)
#plt.gca().add_patch(uncertainty_sie)
plt.plot(proxySST,proxyseaiceedge,'ko',fillstyle='none',markersize=69,alpha=0.3)
for i in range(len(seaiceedge_PMIP3)):
    plt.plot(SSTtest_PMIP3[i],seaiceedge_PMIP3[i],marker='^',markersize=8,markeredgecolor='k',color=colorsPMIP3[i],zorder=4)
for i in range(len(seaiceedge_PMIP4)):
    plt.plot(SSTtest_PMIP4[i],seaiceedge_PMIP4[i],marker='s',markersize=8,markeredgecolor='k',color=colorsPMIP4[i])
# no for loop for love because two different symbols
plt.plot(SSTtest_LOVE[0],seaiceedge_LOVE[0],marker='X',markersize=8,color=colorsLOVE[0],markeredgecolor='k')
plt.plot(SSTtest_LOVE[1],seaiceedge_LOVE[1],marker='P',markersize=8,color=colorsLOVE[0],markeredgecolor='k')
plt.plot(proxySST,proxyseaiceedge,marker='o',markersize=8,color='white',markeredgecolor='k',zorder=4)
# not sure what to do about this symbol
#plt.plot(SSTtest_LOVE[2],seaiceedge_LOVE[2],marker='X',markersize=8,color='black')
# no proxy sea ice edge yet
#linear fit
plt.plot(SSTtest_total,fit_sie,'grey')

# labeling
plt.grid()
plt.xlabel('SST (˚ C)',fontweight='bold',fontsize=15)
plt.ylabel('Sea ice edge (˚ S)',fontweight='bold',fontsize=15)
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', direction="in", length=7, width=3, color="black")



# plotting sea ice area

plt.subplot(3,2,2)
plt.tick_params(bottom=False, top=True, left=False, right=False)
#plt.gca().add_patch(uncertainty_sia)
plt.plot(proxySST,proxyseaicearea,'ko',fillstyle='none',markersize=68,alpha=0.3)
for i in range(len(seaiceedge_PMIP3)):
    plt.plot(SSTtest_PMIP3[i],seaicearea_PMIP3[i],marker='^',markersize=8,color=colorsPMIP3[i],markeredgecolor='k',zorder=4)
for i in range(len(seaiceedge_PMIP4)):
    plt.plot(SSTtest_PMIP4[i],seaicearea_PMIP4[i],marker='s',markersize=8,color=colorsPMIP4[i],markeredgecolor='k')
# no for loop for love because two different symbols
plt.plot(SSTtest_LOVE[0],seaicearea_LOVE[0],marker='P',markersize=8,color=colorsLOVE[0],markeredgecolor='k')
plt.plot(SSTtest_LOVE[1],seaicearea_LOVE[1],marker='X',markersize=8,color=colorsLOVE[0],markeredgecolor='k')
# not sure what to do about this symbol
#plt.plot(SSTtest_LOVE[2],seaicearea_LOVE[2],marker='X',markersize=8,color='black')
# proxy
plt.plot(proxySST,proxyseaicearea,marker='o',markersize=8,color='white',markeredgecolor='k')
#linear fit
plt.plot(SSTtest_total,fit_sia,'grey',zorder=0)
# labeling
plt.grid()
plt.xlabel('SST (˚ C)',fontweight='bold',fontsize=15)
plt.ylabel('Sea ice extent (10$^6$ km$_2$) ',fontweight='bold',fontsize=15)
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', direction="in", length=7, width=3, color="black")
plt.legend(handles=[PMIP3_legend,PMIP4_legend,weakNA_legend,weakNA_AB_legend,Proxy_legend,Proxy_uncertainty],ncol=2)


# plotting zonal averages

plt.subplot(3,2,(3,4))
# two different ways of doing the same thing
for i in range(len(PMIP3_zonal)):
    plt.plot(latitude,PMIP3_zonal[i],color=PMIP3_zonal_colors[i])
for y,c,l in zip(LOVE_zonal, LOVE_zonal_colors,LOVE_linestyle):
    plt.plot(latitude,y,color=c,linestyle=l)

proxy = pd.read_excel('~/Desktop/UNSW/Table_Recap_LGM_Final.xlsx',engine="openpyxl",sheet_name='Data')
data = proxy[['Latitude', 'LGM']]
data = data[3:]
data = data.dropna()
data = data.reset_index(drop=True)

plt.scatter(data.Latitude,data.LGM,color='darkgray')
ave = data
ave.Latitude = ave.Latitude.round()
ave = ave.astype(float)
proxyavg = ave.groupby('Latitude').mean().reset_index()
plt.plot(proxyavg.Latitude,proxyavg.LGM,color='darkgray')
plt.hlines(y=-2,xmin=-60,xmax=-35,linestyles='dashed',color='grey')
plt.xlim(-75,-35)
plt.ylim(-5,25)
plt.ylabel('SST (˚ C) ',fontweight='bold',fontsize=15)
#plt.xlabel('Latitude (˚ S)',fontweight='bold',fontsize=15)
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', direction="in", length=7, width=3, color="black")
plt.grid()
plt.legend(handles=[CNRM_legend,MIROC_legend,GISS_legend,IPSL_legend,MPI_legend,MRI_legend,FGOALS_legend,AWI_legend,LOVE_legend,CESM_legend,CCSM4_legend,proxy_legend,LOVE1_legend,LOVE2_legend,degree],ncol=2)
plt.title('PMIP3 models and LOVECLIM sensitivity runs',fontweight='bold',fontsize=15)

# PMIP4

plt.subplot(3,2,(5,6))
for y,c in zip(PMIP4_zonal, PMIP4_zonal_colors):
    plt.plot(latitude,y,color=c)

proxy = pd.read_excel('~/Desktop/UNSW/Table_Recap_LGM_Final.xlsx',engine="openpyxl",sheet_name='Data')
data = proxy[['Latitude', 'LGM']]
data = data[3:]
data = data.dropna()
data = data.reset_index(drop=True)
plt.scatter(data.Latitude,data.LGM,color='darkgray')


# for some reason was changing the OG data frame so doing changes here
# cleaning up proxy data for zonally averaged line
# cleaning up proxy data for zonally averaged line
ave = data
ave.Latitude = ave.Latitude.round()
ave = ave.astype(float)
proxyavg = ave.groupby('Latitude').mean().reset_index()
plt.plot(proxyavg.Latitude,proxyavg.LGM,color='darkgray')
plt.hlines(y=-2,xmin=-60,xmax=-35,linestyles='dashed',color='grey')
plt.xlim(-75,-35)
plt.ylim(-5,25)
plt.ylabel('SST (˚ C) ',fontweight='bold',fontsize=15)
plt.xlabel('Latitude (˚ S)',fontweight='bold',fontsize=15)
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', direction="in", length=7, width=3, color="black")
plt.grid()
plt.title('PMIP4 models',fontweight='bold',fontsize=15)
plt.show()
plt.tight_layout()
#plt.savefig('Figures/Figure2.pdf')
