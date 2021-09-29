import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/'
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

mapp = pd.read_excel(path + '2021-09-17-tube-welds-mapped.xlsx')

mapp_cut = mapp[['gas-pressure', 'Q_source', 'score']].copy().sort_values(by='Q_source')
mapp_cut_final = mapp_cut.sort_index().iloc[282:]
mapp_cut_final.columns = ['gas','q-source','score']
mapp_cut.columns = ['gas','q-source','score']
mapp_cut_pass = mapp_cut.loc[mapp_cut['score'] == 1]
mapp_cut_okay= mapp_cut.loc[mapp_cut['score'] == 0.5]
mapp_cut_fail = mapp_cut.loc[mapp_cut['score'] == 0]

fig, ax = plt.subplots(figsize = (12,8))

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.set_xlabel('Mechanical gas pressure [Pa]', fontsize = 22.5)
ax.set_ylabel('Average current [A]', fontsize = 22.5)
ax.set_title('Experimental results', fontsize = 24)
ax.tick_params(axis='both', size=7.5)
ax.tick_params(which='minor', size=5)

ax.xaxis.set_major_locator(MultipleLocator(100))
ax.xaxis.set_minor_locator(MultipleLocator(50))

ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.125))

ax.set_ylim(5.25,7)
ax.set_xlim(900,1900)

ax.plot(mapp_cut_pass['gas'],mapp_cut_pass['q-source'], 
        marker='o', markersize=18, color='#C52B38',linewidth=0, 
        label = 'Successful', alpha=0.2)
ax.plot(mapp_cut_okay['gas'],mapp_cut_okay['q-source'], 
        marker='^', markersize=18, color='#4CB91D',linewidth=0, 
        label = 'Flawed', alpha=0.2)
ax.plot(mapp_cut_fail['gas'],mapp_cut_fail['q-source'], 
        marker='s', markersize=18, color='#1b34a9',linewidth=0, 
        label = 'Failed', alpha=0.2)
ax.plot(mapp_cut_final['gas'],mapp_cut_final['q-source'], 
        marker='o', markersize=25, linewidth=0,
        mfc='none', mec='k',
        label = 'Final', alpha=1)

legend = ['Failed', 'Flawed', 'Successful']
ax.legend(legend, prop={'size': 18}, framealpha=1)

rectangle = patches.Rectangle((1094.9, 5.95), 398.144, 0.9, linestyle = '--', 
                         linewidth=2.5, edgecolor='grey', facecolor='none')

ax.vlines(1194, ymin=5.25,ymax=7, linestyle='-.',linewidth=2.5, color='k')
ax.vlines(1742, ymin=5.25,ymax=7, linestyle='-.',linewidth=2.5, color='k')

ax.add_patch(rectangle)

ax.figure.savefig(save + '2021-09-17-tube-mapp.pdf', bbox_inches='tight', dpi=300)

