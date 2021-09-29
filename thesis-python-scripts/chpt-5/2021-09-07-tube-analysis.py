import math 
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  
import seaborn as sns

from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib
import matplotlib.pyplot as plt
from pylab import rcParams
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

tube = pd.read_csv(path + '2020-12-08-tube-welds-2275-raw-titled.csv')
# normalizing
tmp = tube.values 
min_max_scaler = preprocessing.MinMaxScaler()
tmp_scaled = min_max_scaler.fit_transform(tmp)
tube = pd.DataFrame(tmp_scaled)
tube.columns = ['sample', 'electrode-tool-setting','machine-gas','gas-flow','init-gas-pressure',
                'final-gas-pressure','main-current','background','speed','init-current',
                'electrode-position','warm-up','total-weld-with-electrode','score','binary']
#tube.drop(columns=['machine-gas','speed','init-current','electrode-position','warm-up'], inplace=True)

data_X = tube[['electrode-tool-setting','gas-flow','init-gas-pressure',
          'final-gas-pressure','main-current','background','total-weld-with-electrode']]
#         ,'machine-gas','speed','init-current','electrode-position','warm-up']]
data_check = data_X.copy()

data_X.rename(columns={'electrode-tool-setting' : 'Electrode tool\nsetting',
                       'gas-flow' : 'Gas flow',
                       'init-gas-pressure' : 'Mechanical gas pressure',
                       'final-gas-pressure' : 'Weld gas pressure',
                       'main-current' : 'Main current',
                       'background' : 'Background current',
                       'total-weld-with-electrode' : 'Total welds with electrode'}, inplace=True)

data_Y = tube['score']

x = StandardScaler().fit_transform(data_X)

pca = PCA(n_components=7)
principalComponents = pca.fit_transform(x)
exp_var_cumul = np.cumsum(pca.explained_variance_ratio_)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
princiaplDf = pd.DataFrame(data = principalComponents,
                           columns = ['principal component 1',
                                      'principal component 2'])

finalDf = pd.concat([princiaplDf, data_Y], axis=1)
#print(pca.explained_variance_ratio_)

arrow = pd.DataFrame(np.transpose(pca.components_[0:2, :]))

fig = plt.figure(figsize = (12,16))
ax = fig.add_subplot(2,1,1) 
ax.set_xlabel('Principal component 2', fontsize = 15)
ax.set_ylabel('Principal component 1', fontsize = 15)
ax.set_title('Principal component biplot', fontsize = 20)
targets = [0, 0.5, 1]
colors = ['#1b34a9', '#4CB91D', '#C52B38']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['score'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 2'], 
               finalDf.loc[indicesToKeep, 'principal component 1'], 
               c = color, s = 40, alpha=0.95)#, edgecolors='k')

legend = ['Failed', 'Flawed', 'Successful']
ax.legend(legend, prop={'size': 18})
#ax.grid()

for i in range(7):
        #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
        plt.arrow(0, 0, arrow.iloc[i,0] * 1, arrow.iloc[i,1] * 1, linewidth = 2,
                  color='k', alpha = 0.5, length_includes_head=True, shape='full',
                  head_width=0.08, head_length=0.2)
#        plt.text(arrow.iloc[i,0] * 1.25, arrow.iloc[i,1] * 1.25, data_X.columns[i],
#                 color = 'k', ha = 'center', va = 'center',fontsize=15)
   
ax = fig.add_subplot(2,1,2)
ax.set_xlabel('Principal component 2', fontsize = 15)
ax.set_ylabel('Principal component 1', fontsize = 15)
ax.set_title('Zoomed in principal component biplot', fontsize = 20)
ax.set_ylim(-1,1)
ax.set_xlim(-1,1)
targets = [0, 0.5, 1]
colors = ['#1b34a9', '#4CB91D', '#C52B38']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['score'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 2'], 
               finalDf.loc[indicesToKeep, 'principal component 1'], 
               c = color, s = 40, alpha=0.95)#, edgecolors='k')

for i in range(7):
        #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
        plt.arrow(0, 0, arrow.iloc[i,0] * 1, arrow.iloc[i,1] * 1, linewidth = 1,
                  color='k', alpha = 0.5, length_includes_head=True, shape='full',
                  head_width=0.04, head_length=0.05, label=data_X.columns[i])
        plt.text(arrow.iloc[i,0] * 1.15, arrow.iloc[i,1] * 1.15, data_X.columns[i], # '{x}'.format(x=i),#
                 color = 'k', ha = 'center', va = 'center',fontsize=20, label=data_X.columns[i])

#ax.legend()

ax.figure.savefig(save + '2021-09-12-PCA.pdf', bbox_inches='tight', dpi=300)
     