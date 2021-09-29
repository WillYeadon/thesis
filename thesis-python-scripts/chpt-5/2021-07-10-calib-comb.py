import math 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import scipy
import scipy.stats

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

data = pd.read_csv(path + '2021-07-10-Calib-Shunt-Comb.csv')
data.columns = ['current1','s1','s2','s3','sA','r1','r2','r3','rA','current2',
                'c1','c2','c3','cA','current3','p','h']

x = data['current2'].dropna()
y = data['cA'].dropna()
m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
print('Gradient\t' + str(m))
print('Intercept\t' + str(b))
print('r_sq\t' + str(r_value))
xRange = np.arange(0,10,0.2)


fig, ax = plt.subplots()
ax.plot(data['current1'],data['s1'], marker='.', markersize=5, color='k',linewidth=0)
ax.plot(data['current1'],data['s2'], marker='*', markersize=5, color='k',linewidth=0)
ax.plot(data['current1'],data['s3'], marker='^', markersize=5, color='k',linewidth=0)

ax.plot(data['current1'],data['r1'], marker='.', markersize=5, color='red',linewidth=0)
ax.plot(data['current1'],data['r2'], marker='*', markersize=5, color='red',linewidth=0)
ax.plot(data['current1'],data['r3'], marker='^', markersize=5, color='red',linewidth=0)

ax.plot(data['current2'],data['c1'], marker='.', markersize=5, color='blue',linewidth=0)
ax.plot(data['current2'],data['c2'], marker='*', markersize=5, color='blue',linewidth=0)
ax.plot(data['current2'],data['c3'], marker='^', markersize=5, color='blue',linewidth=0)

ax.plot(data['current1'],data['rA'], marker='.', markersize=0, color='red',linewidth=1, label = 'RAL')
ax.plot(data['current1'],data['sA'], marker='.', markersize=0, color='k',linewidth=1, label = 'Sheffield 1')
ax.plot(data['current2'],data['cA'], marker='.', markersize=0, color='blue',linewidth=1, label = 'Sheffield 2')

ax.plot(xRange, m*xRange + b, color='blue', linestyle='--', label = 'Sheffield 2 fit', alpha=1)

ax.plot(data['current3'],data['p'], linestyle=':', color='magenta',linewidth=2.5, label = 'Idealized')
ax.plot(data['current1'],data['h'], marker='.', markersize=8, color='k',linewidth=0, label = 'Run #1')
ax.plot(data['current1'],data['h'], marker='*', markersize=8, color='k',linewidth=0, label = 'Run #2')
ax.plot(data['current1'],data['h'], marker='^', markersize=8, color='k',linewidth=0, label = 'Run #3')

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25))

ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))

ax.set_xlim(0,9.5)
ax.set_ylim(0,9.5)

ax.legend(loc='best')

ax.set_xlabel('Set Current [A]')
ax.set_ylabel('Measured Output Current [A]')

fig.savefig(save + '2021-07-10-current-comp.pdf', bbox_inches='tight', dpi=300)
