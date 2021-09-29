import numpy as np 
import pandas as pd
from pandas import DataFrame
import seaborn as sns
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'
data = pd.read_csv(path + '2021-09-26-tube-data-300.csv', header=None)

def plot(data, name='NO-NAME'):
    pressures = [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20]
    qsources = [36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5,
                32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28]
    
    #cmap = sns.color_palette("Pastel1", 4)
    cmap = sns.color_palette("hls", 4)  
    hm = sns.heatmap(data, cmap = cmap, ax = ax,
                     linewidths=1, linecolor='k',
                     xticklabels=pressures, yticklabels=qsources)
    hm.set_yticklabels(hm.get_yticklabels(), rotation=0, fontsize=14)
    hm.set_xticklabels(hm.get_xticklabels(), rotation=0, fontsize=14)
    hm.set_xlabel("Internal pressure [a.u.]", fontsize=20)
    hm.set_ylabel("$Q_{Source}$ [a.u.]", fontsize=20)
    hm.set_title('Tube weld simulation results at ' + name, fontsize=22)
    
    colorbar = ax.collections[0].colorbar 
    r = colorbar.vmax - colorbar.vmin 
    colorbar.set_ticks([colorbar.vmin + r / 4 * (0.5 + i) for i in range(4)])
    colorbar.set_ticklabels(['Lack of\nPenetration',
                             'Burst',
                             'Collapsed',
                             'Successful'])
    colorbar.ax.tick_params(labelsize=14)

fig, ax = plt.subplots(figsize=(9,10))  
data = pd.read_csv(path + '2021-09-26-tube-data-300.csv', header=None)
plot(data, '0.3 [mm]')
ax.figure.savefig(save + '2021-09-26-tubeScore-300.pdf', bbox_inches='tight', dpi=300)
     
fig2, ax = plt.subplots(figsize=(10,10))  
data = pd.read_csv(path + '2021-09-26-tube-data-350.csv', header=None)
plot(data, '0.35 [mm]')
ax.figure.savefig(save + '2021-09-26-tubeScore-350.pdf', bbox_inches='tight', dpi=300)

fig3, ax = plt.subplots(figsize=(10,10))  
data = pd.read_csv(path + '2021-09-26-tube-data-250.csv', header=None)
plot(data, '0.25 [mm]')
ax.figure.savefig(save + '2021-09-26-tubeScore-250.pdf', bbox_inches='tight', dpi=300)


 