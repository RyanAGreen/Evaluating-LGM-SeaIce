import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import xarray as xr
import matplotlib as mpl
from scipy import stats

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
proxyseaicearea = [15.35]
proxyseaiceedge = [-61.5]
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


correlation_matrix_sie = np.corrcoef(SSTtest_total, seaiceedge_total)
correlation_xy_sie = correlation_matrix_sie[0,1]
r_squared_sie = correlation_xy_sie**2
print('Rsquared for sea ice edge vs SST is ',r_squared_sie)

correlation_matrix_sia = np.corrcoef(SSTtest_total, seaicearea_total)
correlation_xy_sia = correlation_matrix_sia[0,1]
r_squared_sia = correlation_xy_sia**2
print('Rsquared for sea ice area vs SST is ',r_squared_sia)

# checking answer
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)
     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    correlation = np.corrcoef(x, y)[0,1]

     # r
    results['correlation'] = correlation
     # r-squared
    results['determination'] = correlation**2

    return results
print(polyfit(SSTtest_total,seaicearea_total,1))


from sklearn.metrics import r2_score

R_square = r2_score(SSTtest_total, seaicearea_total)
print('Coefficient of Determination', R_square)

slope, intercept, r_value, p_value, std_err = stats.linregress(SSTtest_total, seaicearea_total)
print("r-squared:", r_value**2)

slope, intercept, r_value, p_value, std_err = stats.linregress(SSTtest_total, seaiceedge_total)
print("r-squared:", r_value**2)
