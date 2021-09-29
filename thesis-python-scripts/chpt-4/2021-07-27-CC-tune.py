import math 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import scipy.stats
from sklearn.linear_model import LinearRegression
plt.rcParams["figure.figsize"] = (6, 4)

path = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-4/images/'

CC_data = pd.read_excel(path + '2021-08-02-CC-data.xlsx', 'data')
CC_data.columns = ['Q','speed','l','omega','CC']
CC_data_X = CC_data.drop(columns='CC').to_numpy()
CC_data_Y = CC_data['CC']
regr = LinearRegression().fit(CC_data_X, CC_data_Y)
regrCoeff = regr.coef_
CC_data_X_trans = np.dot(CC_data_X, regrCoeff) 
print('Coefficients: \n', regr.coef_)

fig, ax = plt.subplots()
ax.plot(CC_data_X_trans, CC_data_Y, linewidth=0, marker='.', color='k', markersize=15)
x = CC_data_X_trans
y = CC_data_Y
xRange = np.arange(80,180,1)
m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
ax.plot(xRange, m*xRange + b, color='red', linestyle='--', label = 'Fit', alpha=1)
ax.plot(xRange, xRange, color='blue', linestyle=':')
plt.xlim(80,170)
plt.ylim(90,200)
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(5))

ax.set_ylabel('Tuned $C_C$ Computational Factor')
ax.set_xlabel('Calculated $C_C$ Computational Factor')

save = 'C:/Users/php17wgy/Documents/thesis-chapter-4/images/'
ax.figure.savefig(save + '2021-07-27-CC-tune.pdf', bbox_inches='tight', dpi=300)