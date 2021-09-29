import math 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

current = pd.read_excel(path + '2021-07-10-Shunt-Data-Orig.xlsx', '6A-current2')
current.columns = ['time','current']
current['time100'] = current['time'].values - 98.25
average = current['current'].mean()
print(average)

fig, ax = plt.subplots()
ax.plot(current['time100'],current['current'], color='blue', linewidth=0.5)
ax.hlines(average,0,30,color='red',linestyles='-', alpha=0.85)
ax.axvline(x=5, color = 'k', linestyle = ':', alpha=0.85)
ax.axvline(x=15, color = 'k', linestyle = ':', alpha=0.85)
ax.set_xlim(0,30)
ax.set_ylim(5.8,6.8)

ax.xaxis.set_minor_locator(MultipleLocator(1.25))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))

ax.set_xlabel('Time [ms]')
ax.set_ylabel('Output Current [A]')

fig.savefig(save + '2021-07-10-current-shape.pdf', bbox_inches='tight', dpi=300)
