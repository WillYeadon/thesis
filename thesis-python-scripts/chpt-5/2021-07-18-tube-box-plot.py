import math 
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import preprocessing

import matplotlib
import matplotlib.pyplot as plt
from pylab import rcParams
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

tube = pd.read_csv(path + '2020-12-08-tube-welds-2275-raw-titled.csv')
tmp = tube.values 
min_max_scaler = preprocessing.MinMaxScaler()
tmp_scaled = min_max_scaler.fit_transform(tmp)
tube = pd.DataFrame(tmp_scaled)
tube.columns = ['sample', 'Electrode\ntool\nsetting','Machine\ngas','Gas\nflow','Mechanical\ngas\npressure',
                'Weld\ngas\npressure','Main\ncurrent','Back-\nground\ncurrent','Rotational\nhead\nspeed','Initial\ncurrent',
                'Electrode\nposition','Warmup\nwelds','Total\nwelds with\nelectrode','Score','Binary']
#tube.columns = ['sample', 'electrode-tool-setting','machine-gas','gas-flow','init-gas-pressure',
#                'final-gas-pressure','main-current','background','speed','init-current',
#                'electrode-position','warm-up','total-weld-with-electrode','score']

plotTube = tube[['Electrode\ntool\nsetting','Machine\ngas','Gas\nflow','Mechanical\ngas\npressure',
                'Weld\ngas\npressure','Main\ncurrent','Back-\nground\ncurrent','Rotational\nhead\nspeed','Initial\ncurrent',
                'Electrode\nposition','Warmup\nwelds','Total\nwelds with\nelectrode']]

rcParams['figure.figsize'] = 15, 8
sns.set(font_scale=1.5)
sns.set_style("ticks")
fig = sns.boxplot(data=plotTube)
fig.grid(False)
fig.set_xlabel('Parameter', size = 20)
fig.set_ylabel('Normalized value', size = 20)

fig.figure.savefig(save + '2021-07-20-parameters.pdf', bbox_inches='tight', dpi=300)