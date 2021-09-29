import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import scipy.stats
from sklearn.linear_model import LinearRegression

def addScore(baseArr, wt, path):
    for i in range(280,365,5):
        for j in range(80, 125, 5):       
            score = 0
            q = i / 10
            p = j / 100
            q2 = '{:.1f}'.format(q)
            p2 = '{:.2f}'.format(p)
            case = q2.strip() + '-' + p2.strip()     
            with open(path + case + '.txt') as var:
                for line in var:
                    line = line.rstrip()
                    if line == 'pass':
                        score = 1
                    else:
                        score = 0 
            newrow = np.array([wt, float(q2), float(p2), score])
            baseArr = np.vstack([baseArr, newrow])
            
    return baseArr
    
newArr = np.array([0,0,0,0])

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults25/'    
newArr = addScore(newArr, 0.25, path)
path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults275/'
newArr = addScore(newArr, 0.275, path)
path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults30/'
newArr = addScore(newArr, 0.30, path)
path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults325/'
newArr = addScore(newArr, 0.325, path)
path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults35/'    
newArr = addScore(newArr, 0.35, path)
data = np.delete(newArr, 0, 0)

data_X = data[:, [0,1,2]]
data_Y = data[:,3]

successful = np.greater(data[:,3], 0.5)
extraction = data[successful]
finalSuccessful = np.delete(extraction, 3, axis=1)

failed = np.less(data[:,3], 0.5)
extraction = data[failed]
finalFailed = np.delete(extraction, 3, axis=1)

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

fig = plt.figure(figsize=(9,8))
ax = Axes3D(fig)

ax.scatter(finalSuccessful[:,0], finalSuccessful[:,1], finalSuccessful[:,2], 
           marker='x', color='red', s=100)
ax.scatter(finalFailed[:,0], finalFailed[:,1], finalFailed[:,2], 
           marker='o', color='blue', s=25, alpha=0.5)

ax.set_xlabel('Wall thickness [mm]', fontsize=20)#, rotation=150)
ax.set_ylabel('$Q_{Source}$ [a.u.]', fontsize=20)
ax.set_zlabel('Internal pressure [a.u.]',fontsize=20)#, rotation=60)

ax.xaxis.labelpad=20
ax.yaxis.labelpad=20
ax.zaxis.labelpad=20

#ax.xaxis.pane.fill = False
#ax.yaxis.pane.fill = False
#ax.zaxis.pane.fill = False

ax.xaxis.pane.set_edgecolor('k')
ax.yaxis.pane.set_edgecolor('k')
ax.zaxis.pane.set_edgecolor('k')

ax.w_xaxis.line.set_color('k')
ax.w_yaxis.line.set_color('k')
ax.w_zaxis.line.set_color('k')


# Axes color is red
ax.set(facecolor='white')
# Figure color is green
fig.set(facecolor='white')

#plt.tight_layout()
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'
plt.savefig(save + '2021-09-27-tubeScore-stacked.pdf', bbox_inches='tight', dpi=300)

plt.show()