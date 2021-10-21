import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib as mpl

curl = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/windstresscurlPMIP4.nc')
PMIP3_bad = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/windstresscurlbadLOVE1.nc')
PMIP3_good = xr.open_dataset('~/Desktop/UNSW/PMIP4Models/windstresscurlgoodLOVE1.nc')

curl['CESM'] = curl.CESM*0.1

LOVE = [PMIP3_good.CURLNAW,PMIP3_good.CURLNAWSOW]
PMIP3 = [PMIP3_bad.CURLCNRM,PMIP3_bad.CURLGISS,PMIP3_bad.CURLIPSL,PMIP3_good.CURLMIROC,PMIP3_good.CURLMPI,PMIP3_good.CURLMRI,PMIP3_good.CURLFGOALS,PMIP3_bad.CURLCCSM4]
PMIP4 = [curl.MIROC,curl.IPSL,curl.MPI,curl.AWI,curl.LOVE,curl.CESM,curl.CCSM4]

Thermodynamically_PMIP4 = [PMIP4[0],PMIP4[1],PMIP4[5],PMIP4[6]]
Dynamically_PMIP4 = [PMIP4[2],PMIP4[3],PMIP4[4]]

Thermodynamically_PMIP3 = [PMIP3[0],PMIP3[1],PMIP3[2],PMIP3[3],PMIP3[7]]
Dynamically_PMIP3 = [PMIP3[4],PMIP3[5],PMIP3[6]]

# these are exact but I rounded
# PMIP3_sie = [-75.5,-65.5,-70,-66.5,-65,-62.5,-61.5,-55.5]
# PMIP4_sie = [-75.5,-70,-65,-62,-59.5,-57.5,-53]
# LOVE_sie = [-59.5,58.5]

# rounded down
PMIP3_sie = [-75.5,-65.5,-69.5,-66.5,-64.5,-62.5,-61.5,-55.5]
PMIP4_sie = [-75.5,-69.5,-64.5,-61.5,-59.5,-57.5,-52.5]
LOVE_sie = [-59.5,-58.5]
Thermodynamically_PMIP4_sie = [PMIP4_sie[0],PMIP4_sie[1],PMIP4_sie[5],PMIP4_sie[6]]
Dynamically_PMIP4_sie = [PMIP4_sie[2],PMIP4_sie[3],PMIP4_sie[4]]

Thermodynamically_PMIP3_sie = [PMIP3_sie[0],PMIP3_sie[1],PMIP3_sie[2],PMIP3_sie[3],PMIP3_sie[7]]
Dynamically_PMIP3_sie = [PMIP3_sie[4],PMIP3_sie[5],PMIP3_sie[6]]


# rounded up
# PMIP3_sie = [-75.5,-65.5,-70.5,-66.5,-65.5,-62.5,-61.5,-55.5]
# PMIP4_sie = [-75.5,-71.5,-65.5,-62.5,-59.5,-57.5,-53.5]
# LOVE_sie = [-59.5,-58.5]

PMIP3_colors = ['#ff7f0e','#e377c2','#1f77b4','#bcbd22','#9467bd','#2ca02c','#17becf','black']
PMIP3_colors.reverse()

PMIP4_colors = ['#ff7f0e','#8c564b','#F7DC6F','#B8255F','#bcbd22','#2ca02c','#9467bd']
PMIP4_colors.reverse()
LOVE_colors = ['#F7DC6F','#F7DC6F']
LOVE_linestyle = ['dotted','dashdot']



# set up legends

PMIP3names = ['CNRM','GISS-E2-R','IPSL-CM5A-LR','MIROC-ESM-P','MPI-ESM-P','MRI-CGCM3','FGOALS-G2','CCSM4']
PMIP4names= ['MIROC-ES2L','IPSL-CM5A2','MPI-ESM1-2','AWI-ESM-1','LOVECLIM','CESM1.2','UoT-CCSM4']
LOVEnames = ['weakNA','weakNA_AB']
Thermodynamically_PMIP4_names = [PMIP4names[0],PMIP4names[1],PMIP4names[5],PMIP4names[6]]
Dynamically_PMIP4_names = [PMIP4names[2],PMIP4names[3],PMIP4names[4]]

Thermodynamically_PMIP3_names= [PMIP3names[0],PMIP3names[1],PMIP3names[2],PMIP3names[3],PMIP3names[7]]
Dynamically_PMIP3_names = [PMIP3names[4],PMIP3names[5],PMIP3names[6]]


# colors are set from most sea ice to least
PMIP3_colors = ['#ff7f0e','#e377c2','#1f77b4','#bcbd22','#9467bd','#2ca02c','#17becf','black']
PMIP3_colors.reverse()
PMIP4_colors = ['#ff7f0e','#8c564b','#F7DC6F','#B8255F','#bcbd22','#2ca02c','#9467bd']
PMIP4_colors.reverse()
LOVE_colors = ['#F7DC6F','#F7DC6F']

Thermodynamically_PMIP4_colors = [PMIP4_colors[0],PMIP4_colors[1],PMIP4_colors[5],PMIP4_colors[6]]
Dynamically_PMIP4_colors = [PMIP4_colors[2],PMIP4_colors[3],PMIP4_colors[4]]

Thermodynamically_PMIP3_colors= [PMIP3_colors[0],PMIP3_colors[1],PMIP3_colors[2],PMIP3_colors[3],PMIP3_colors[7]]
Dynamically_PMIP3_colors = [PMIP3_colors[4],PMIP3_colors[5],PMIP3_colors[6]]


#PMIP3
CNRM_leg = mlines.Line2D([], [], color=PMIP3_colors[0], linestyle ='solid',label = PMIP3names[0])
GISS3_leg = mlines.Line2D([], [], color=PMIP3_colors[1], linestyle ='solid',label = PMIP3names[1])
IPSL3_leg = mlines.Line2D([], [], color=PMIP3_colors[2], linestyle ='solid',label = PMIP3names[2])
MIROC3_leg = mlines.Line2D([], [], color=PMIP3_colors[3], linestyle ='solid',label = PMIP3names[3])
MPI3_leg = mlines.Line2D([], [], color=PMIP3_colors[4], linestyle ='solid',label = PMIP3names[4])
MRI_leg = mlines.Line2D([], [], color=PMIP3_colors[5], linestyle ='solid',label = PMIP3names[5])
FGOALS_leg = mlines.Line2D([], [], color=PMIP3_colors[6], linestyle ='solid',label = PMIP3names[6])
CCSM4_leg = mlines.Line2D([], [], color=PMIP3_colors[7], linestyle ='solid',label = PMIP3names[7])

#PMIP4
MIROC4_leg = mlines.Line2D([], [], color=PMIP4_colors[0], linestyle ='solid',label = PMIP4names[0])
IPSL4_leg = mlines.Line2D([], [], color=PMIP4_colors[1], linestyle ='solid',label = PMIP4names[1])
MPI4_leg = mlines.Line2D([], [], color=PMIP4_colors[2], linestyle ='solid',label = PMIP4names[2])
AWI_leg = mlines.Line2D([], [], color=PMIP4_colors[3], linestyle ='solid',label = PMIP4names[3])
LOVE_leg = mlines.Line2D([], [], color=PMIP4_colors[4], linestyle ='solid',label = PMIP4names[4])
CESM_leg = mlines.Line2D([], [], color=PMIP4_colors[5], linestyle ='solid',label = PMIP4names[5])
CCSM4UoT_leg = mlines.Line2D([], [], color=PMIP4_colors[6], linestyle ='solid',label = PMIP4names[6])

#LOVE
weakNA_leg = mlines.Line2D([], [], color=LOVE_colors[0], linestyle ='dotted',label = LOVEnames[0])
weakNA_AB_leg = mlines.Line2D([], [], color=LOVE_colors[1], linestyle ='dashdot',label = LOVEnames[1])

# straight from Figure 2
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
PMIP4legend = mlines.Line2D([], [], color='black', linestyle ='solid',label = 'PMIP4 models')
PMIP3legend = mlines.Line2D([], [], color='black', linestyle ='dashed',label = 'PMIP3 models')

PMIP3_legend = mlines.Line2D([], [], color='black', marker='^', linestyle='None',markersize=8,markerfacecolor='white', label='PMIP3 models')
PMIP4_legend = mlines.Line2D([], [], color='black', marker='s', linestyle='None',markersize=8,markerfacecolor='white', label='PMIP4 models')
weakNA_legend = mlines.Line2D([], [], color='black', marker='P', linestyle='None',markersize=8,markerfacecolor='white', label='LOVECLIM-weakNA')
weakNA_AB_legend = mlines.Line2D([], [], color='black', marker='X', linestyle='None',markersize=8,markerfacecolor='white', label='LOVECLIM-WeakNA_AB')
Proxy_legend = mlines.Line2D([], [], color='black', marker='o', linestyle='None',markersize=8,markerfacecolor='white', label='Estimate from \nproxy SST')
#MMM = mlines.Line2D([], [], color='black', marker='o', linestyle='None',markersize=8, label='Multi-model mean')
densly = (0, (3, 1, 1, 1))
LOVE1_legend = mlines.Line2D([], [], color='#F7DC6F', linestyle ='dotted',label = 'LOVECLIM-weakNA')
LOVE2_legend = mlines.Line2D([], [], color='#F7DC6F', linestyle ='dashdot',label = 'LOVECLIM-weakNA_AB')
PMIP3_line = mlines.Line2D([], [], color='black', linestyle ='solid',label = 'PMIP3 model')
PMIP4_line = mlines.Line2D([], [], color='black', linestyle =densly,label = 'PMIP4 model')

fig,ax = plt.subplots(2,figsize=(10,10),sharex=True)
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.weight"] = "bold"
ymin = -2.5
ymax = 2.5
for i in range(2):
    ax[i].set_xlim(-80,-30)
    ax[i].set_ylim(ymin,ymax)
    ax[i].vlines(-75.5,ymin=ymin,ymax=ymax,colors='k',linestyles='dashdot',label='Antarctic Coast')
    ax[i].grid(ls=':')
    ax[i].tick_params(axis="both", direction="out", length=5, width=3, color="black")
    for axis in ['top','bottom','left','right']:
        ax[i].spines[axis].set_linewidth(3)

for i in range(len(Thermodynamically_PMIP3)):
    ax[0].plot(PMIP3_bad.LAT,Thermodynamically_PMIP3[i],color=Thermodynamically_PMIP3_colors[i],ls = '--')
    ax[0].plot(Thermodynamically_PMIP3_sie[i],Thermodynamically_PMIP3[i].where(curl.LAT==Thermodynamically_PMIP3_sie[i],drop=True),marker='o',color=Thermodynamically_PMIP3_colors[i],markeredgecolor='k')
for i in range(len(Thermodynamically_PMIP4)):
    ax[0].plot(PMIP3_bad.LAT,Thermodynamically_PMIP4[i],color=Thermodynamically_PMIP4_colors[i])
    ax[0].plot(Thermodynamically_PMIP4_sie[i],Thermodynamically_PMIP4[i].where(curl.LAT==Thermodynamically_PMIP4_sie[i],drop=True),marker='o',color=Thermodynamically_PMIP4_colors[i],markeredgecolor='k')

#ax[0].legend(handles=[CNRM_leg,GISS3_leg,IPSL3_leg,MIROC3_leg,MPI3_leg,MRI_leg,FGOALS_leg,CCSM4_leg])

for i in range(len(Dynamically_PMIP3)):
    ax[1].plot(PMIP3_bad.LAT,Dynamically_PMIP3[i],color=Dynamically_PMIP3_colors[i],ls = '--')
    ax[1].plot(Dynamically_PMIP3_sie[i],Dynamically_PMIP3[i].where(curl.LAT==Dynamically_PMIP3_sie[i],drop=True),marker='o',color=Dynamically_PMIP3_colors[i],markeredgecolor='k')
for i in range(len(Dynamically_PMIP4)):
    ax[1].plot(PMIP3_bad.LAT,Dynamically_PMIP4[i],color=Dynamically_PMIP4_colors[i])
    ax[1].plot(Dynamically_PMIP4_sie[i],Dynamically_PMIP4[i].where(curl.LAT==Dynamically_PMIP4_sie[i],drop=True),marker='o',color=Dynamically_PMIP4_colors[i],markeredgecolor='k')
#ax[1].legend(handles=[MIROC4_leg,IPSL4_leg,MPI4_leg,AWI_leg,LOVE_leg,CESM_leg,CCSM4UoT_leg])

for i in range(len(LOVE)):
    ax[1].plot(PMIP3_bad.LAT,LOVE[i],color=LOVE_colors[i],ls = LOVE_linestyle[i])
    ax[1].plot(LOVE_sie[i],LOVE[i].where(curl.LAT==LOVE_sie[i],drop=True),marker='o',color=LOVE_colors[i],markeredgecolor='k')

# ax[2].legend(handles=[weakNA_leg,weakNA_AB_leg])

# for i in range(3):
#     ax[i].text(-76.5,0.5,'Antarctic coast',rotation='vertical',fontsize=10,fontweight='bold')

ax[1].text(-76.5,0.5,'Antarctic coast',rotation='vertical',fontsize=10,fontweight='bold')
ax[0].text(-76.5,0.5,'Antarctic coast',rotation='vertical',fontsize=10,fontweight='bold')
ax[0].text(-79,2.2,'a)',fontsize=10,fontweight='bold')
ax[1].text(-79,2.2,'b)',fontsize=10,fontweight='bold')
#ax[2].text(-79,2.2,'c)',fontsize=10,fontweight='bold')
ax[0].set_title('Thermodynamically driven models',fontsize=15,fontweight='bold')
ax[0].legend(handles=[PMIP4legend,PMIP3legend,CNRM_legend,MIROC_legend,GISS_legend,IPSL_legend,CESM_legend,CCSM4_legend],ncol=3)

ax[1].set_title('Dynamically driven models',fontsize=15,fontweight='bold')
#ax[2].set_title('LOVECLIM sensitivity runs',fontsize=15,fontweight='bold')
ax[1].legend(handles=[PMIP4legend,PMIP3legend,MPI_legend,MRI_legend,FGOALS_legend,AWI_legend,LOVE_legend,LOVE1_legend,LOVE2_legend],ncol=3)
plt.tight_layout(pad=2.2)

# add labels
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel("Wind Stress Curl (m/s)",fontsize=15,fontweight='bold')
plt.xlabel("Latitude (˚ S)",fontsize=15,fontweight='bold')

plt.show()
#plt.savefig('Figures/Figure5_separate.pdf')
