import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

address_ga = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/ga-melt-ray.txt'
address_tS6 = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/tin-solid-506-Ray.txt'
address_tS7 = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/tin-solid-507-Ray.txt'
address_tm6 = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/tin-melt-506-Ray.txt'
address_tm7 = 'C:/Users/php17wgy/Documents/thesis-chapter-4/data/tin-melt-507-Ray.txt'

def getData(address):
    time, ray, alpha = [], [], []

    with open(address) as file:
        for line in csv.reader(file, delimiter="\t"):
            time.append(int(line[0]))
            ray.append(float(line[1]))
            alpha.append(float(line[2]))      

    time_arr = Norm(np.array(time))
    ray_arr = Norm(np.array(ray))
    alpha_arr = np.array(alpha)
            
    return time_arr, alpha_arr, ray_arr 

def plotData(time, alpha, ray, color, label):
    ax1.plot(time, alpha, color=color, linestyle='--', label=label)  
    ax2.plot(time, ray, color=color)

def Norm(values):
    return values/np.max(values)
      
fig, ax1 = plt.subplots()
ax1.set_xlabel('Normalized time [s]')
ax1.set_ylabel('Alpha fraction')
ax1.tick_params(axis='y')    
ax2 = ax1.twinx()
#ax2.set_ylim(0,1)
ax2.set_ylabel('Normalized Rayleigh number')  
ax2.tick_params(axis='y')

#fig.tight_layout()  


legend_elements = [Line2D([0], [0], color='k', linestyle='--', label='Alpha Fraction'),
                   Line2D([0], [0], color='k', label='Normalized Rayleigh number')]

#melt = True
melt = False

if (melt):
    time_ga, alpha_ga, ray_ga = getData(address_ga)
    plotData(time_ga, alpha_ga, ray_ga, 'r', 'ga')
    time_tm6, alpha_tm6, ray_tm6 = getData(address_tm6)
    plotData(time_tm6, alpha_tm6, ray_tm6, 'b', 'tin1')
    time_tm7, alpha_tm7, ray_tm7 = getData(address_tm7)
    plotData(time_tm7, alpha_tm7, ray_tm7, 'g', 'tin2')
    ax1.legend(handles=legend_elements, loc='lower right')
    plt.savefig('../images/2021-06-23-alpha-ray-melt.pdf', bbox_inches='tight', dpi=300)
#    plt.savefig('../images/2021-06-23-alpha-ray-melt.png', bbox_inches='tight', dpi=300)    
else:
    time_tS6, alpha_tS6, ray_tS6 = getData(address_tS6)
    plotData(time_tS6, alpha_tS6, ray_tS6, 'b', 'tin1')
    time_tS7, alpha_tS7, ray_tS7 = getData(address_tS7)
    plotData(time_tS7, alpha_tS7, ray_tS7, 'g', 'tin2')
    ax1.legend(handles=legend_elements, loc='upper right')
    plt.savefig('../images/2021-06-23-alpha-ray-solid.pdf', bbox_inches='tight', dpi=300)
#    plt.savefig('../images/2021-06-23-alpha-ray-solid.pdf', bbox_inches='tight', dpi=300)


