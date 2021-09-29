import math 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import scipy.stats

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

data = pd.read_csv(path + '2021-07-14-gas-pressures.csv')
data.columns = ['init-gas','final-gas']

def dropFails(data):
    if data['final-gas'] < 2.5:
        return np.NAN
    else:
        return data

dataPass = data.apply(dropFails, axis = 1)
dataPass.reset_index(inplace = True, drop = True)
dataPass = dataPass.dropna()

x = dataPass['init-gas']
y = dataPass['final-gas']
m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
print('Gradient\t' + str(m))
print('Intercept\t' + str(b))
print('r_sq\t' + str(r_value))
xRange = np.arange(0,10,0.2)

fig, ax = plt.subplots()

ax.plot(dataPass['init-gas'],dataPass['final-gas'], marker='.', markersize=15, color='k',linewidth=0, label = 'Data', alpha=0.2)
ax.plot(xRange, m*xRange + b, color='red', linestyle='--', label = 'Fit', alpha=1)

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))
ax.set_xlim(3.5,8)
ax.set_ylim(3.5,8)

ax.set_xlabel('Mechanical Joint Pressure [inAq]')
ax.set_ylabel('Welded Joint Pressure [inAq]')

ax.legend(loc='best')

fig.savefig(save + '2021-07-14-gas-pressure.pdf', bbox_inches='tight', dpi=300)

