import math 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

UTS_el = pd.read_csv(path + '2021-07-07-raw-UTS-data.csv')
sample = UTS_el.iloc[:,0]
UTS = UTS_el.iloc[:,1]
elon = UTS_el.iloc[:,2]

fig1, ax1 = plt.subplots()
HV = pd.read_csv(path + '2021-07-07-raw-HV-data.csv')
distance = HV.iloc[:,0]
HV = HV.iloc[:,1]
ax1.plot(distance,HV,lw=1,marker='.', markersize=10,color='k')
ax1.set(ylabel = 'Hardness Value [HV0.1]', xlabel = 'Distance [mm]')
ax1.axvline(x=-2.5, color = 'k', linestyle = ':')
ax1.axvline(x=-0.8, color = 'k', linestyle = ':')
ax1.axvline(x=1.5, color = 'k', linestyle = ':')
ax1.axvline(x=3.25, color = 'k', linestyle = ':')
ax1.text(-6, 135, 'Base Metal\n    Sleeve', fontsize=12)
ax1.text(-2.5, 160, '  HAZ\nSleeve', fontsize=12)
ax1.text(-0.4, 130, 'Fusion\n Zone\nSleeve', fontsize=12)
ax1.text(1.75, 190, 'HAZ\nTube', fontsize=12)
ax1.text(4, 140, 'Base Metal\n    Tube', fontsize=12)
ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_minor_locator(MultipleLocator(5))
fig1.savefig(save + '2021-07-07-HV.pdf', bbox_inches='tight', dpi=300)


fig2, ax2 = plt.subplots()
ax2.plot(sample,UTS,lw=1,marker='*',markersize=10,color='blue', label = 'UTS')
ax2.tick_params(axis='y', labelcolor='blue', color='blue')
ax2.xaxis.set_major_locator(MultipleLocator(1))
ax2.spines['right'].set_color('red')
ax2.set(ylabel = 'Ultimate Tensile Strength [MPa]', xlabel = 'Specimen Number')

ax3 = ax2.twinx()
ax3.plot(sample,elon,lw=1,marker='.',markersize=10,color='red')
ax3.tick_params(axis='y', labelcolor='red', color='red')
ax3.spines['left'].set_color('blue')
ax3.spines['right'].set_color('red')
ax3.set(ylabel = 'Total Elongation [%]')

fig2.legend(('Ultimate Tensile Strength','Total Elongation'),loc=(0.35,0.175))
fig2.savefig(save + '2021-07-07-UTSel.pdf', bbox_inches='tight', dpi=300)

#plt.plot(HV[0],HV[1])

