import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


proxy = pd.read_excel('~/Desktop/UNSW/Table_Recap_LGM_Final.xlsx',engine="openpyxl",sheet_name='Data')
data = proxy[['Latitude', 'LGM']]
data = data[3:]
data = data.dropna()
data = data.reset_index(drop=True)
column = data['Latitude']

df = data[data.Latitude > -52]
plt.scatter(data.Latitude,data.LGM,color='darkgray')
ave = data
ave.Latitude = ave.Latitude.round()
ave = ave.astype(float)

ave2 = df
ave2.Latitude = ave2.Latitude.round()
ave2 = ave2.astype(float)
column = ave['Latitude']
proxyavg = ave2.groupby('Latitude').mean().reset_index()
#plt.plot(proxyavg.Latitude,proxyavg.LGM,color='darkgray')
#plt.show()
slope = np.polyfit(proxyavg.Latitude,proxyavg.LGM,1)[0]
print(slope)
print(column.min())
#print(proxyavg.LGM.mean())
