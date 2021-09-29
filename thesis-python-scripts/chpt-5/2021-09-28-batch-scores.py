import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

scores = np.zeros((17,9))
save = 'C:/Users/php17wgy/Documents/thesis-chapter-5/images/'

def score(path, transpose=False):
    npI = 0
    for i in range(280,365,5):
        npJ = 0
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
    #                print(line)
                    if line == 'pass':
                        score = 4
                    elif line == 'fail	collapse':
                        score = 3
                    elif line == 'fail	burst':
                        score = 2
                    elif line == 'fail	penn':
                        score = 1
            scores[npI,npJ] = score
#            print(score,' ', end='')
            npJ += 1
        npI += 1
#        print()
    if (transpose):
        udFlip = np.flipud(scores)
        transpose = np.transpose(udFlip)
        return transpose
    return np.flipud(scores)

def plotOne(data, name='NO-NAME'):
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
#    hm.set_title('Tube weld simulation results at ' + name, fontsize=22)
    
    colorbar = ax.collections[0].colorbar 
    r = colorbar.vmax - colorbar.vmin 
    colorbar.set_ticks([colorbar.vmin + r / 4 * (0.5 + i) for i in range(4)])
    colorbar.set_ticklabels(['Lack of\nPenetration',
                             'Burst',
                             'Collapsed',
                             'Successful'])
    colorbar.ax.tick_params(labelsize=14)

def plotMany(data, axis, name='NO-NAME', transpose=False):
    pressures = [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20]
    qsources = [36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5,
                32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28]
    
    #cmap = sns.color_palette("Pastel1", 4)
    cmap = sns.color_palette("hls", 4)  

    if transpose:
        hm = sns.heatmap(data, cmap = cmap, ax = axis,
                         cbar_kws={"orientation": "horizontal"},
                         cbar_ax=cbar_ax,
                         linewidths=1, linecolor='k',
                         yticklabels=pressures, xticklabels=qsources)
    else:
        hm = sns.heatmap(data, cmap = cmap, ax = axis,
                         cbar_kws={"orientation": "horizontal"},
                         cbar_ax=cbar_ax,
                         linewidths=1, linecolor='k',
                         xticklabels=pressures, yticklabels=qsources)
        hm.set_yticklabels(hm.get_yticklabels(), rotation=0, fontsize=14)
    hm.set_xticklabels(hm.get_xticklabels(), rotation=0, fontsize=14)
    hm.set_xlabel("Internal pressure [a.u.]", fontsize=20)
    hm.set_ylabel("$Q_{Source}$ [a.u.]", fontsize=20)
    hm.set_title(name, fontsize=22)

    colorbar = axis.collections[0].colorbar 
    r = colorbar.vmax - colorbar.vmin 
    colorbar.set_ticks([colorbar.vmin + r / 4 * (0.5 + i) for i in range(4)])
    colorbar.set_ticklabels(['Lack of\nPenetration',
                             'Burst',
                             'Collapsed',
                             'Successful'])
    colorbar.ax.tick_params(labelsize=20)


path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults30/'    
fig, ax = plt.subplots(figsize=(9,10))
test_scores = score(path)  
plotOne(score(path), '0.3 [mm]')
ax.figure.savefig(save + '2021-09-26-tubeScore-300.pdf', bbox_inches='tight', dpi=300)


#gridkw = dict(height_ratios=[1, 1])
#fig2, ((ax_1, ax_2), (ax_3, ax_4), (ax_5, ax_6)) = plt.subplots(3,2, figsize=(15,22))#, gridspec_kw=gridkw)
fig2, ((ax_3, ax_4), (ax_5, ax_6)) = plt.subplots(2,2, figsize=(15,22))#, gridspec_kw=gridkw)
sns.set(font_scale=1.5)
#cbar_ax = fig2.add_axes([.92, .15, .03, .7])
cbar_ax = fig2.add_axes([.1, .05, .8, .02])

#path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults30/'    
#plotMany(score(path, True), ax_1, '0.2875 [mm]')

#path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults30/'    
#plotMany(score(path, True), ax_2, '0.3125 [mm]')

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults25/'    
plotMany(score(path), ax_3, '0.25 [mm]')

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults275/'    
plotMany(score(path), ax_4, '0.275 [mm]')

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults325/'    
plotMany(score(path), ax_5, '0.325 [mm]')

path = 'C:/Users/php17wgy/Documents/thesis-chapter-5/data/batchResults35/'    
plotMany(score(path), ax_6, '0.35 [mm]')   

fig.tight_layout()
#fig2.suptitle('Test')
fig2.savefig(save + '2021-09-27-tubeScore-combined.pdf', bbox_inches='tight', dpi=300)


