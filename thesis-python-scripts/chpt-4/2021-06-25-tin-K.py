import numpy as np
import matplotlib.pyplot as plt

def peralTinK(a, b, c, T):
    return (a + b*(T/273.15) - c*(T/273.15)**2)

a = 10.204
b = 32.063
c = 5.686

T_range = np.arange(450, 550, 0.5)
peralT_values = peralTinK(a, b, c, T_range)

# wolfram alpha gives max at 3.98838 C
# = 277.13838 K
fontThesis = {'fontname' : 'Arial', 'fontsize' : 12}

plt.plot(T_range, peralT_values, 'r--', linewidth = 2.0,
         label = 'Function')
plt.tick_params('both', labelsize=10)
plt.xlabel("Temperature " + "[K]", #+ r'[$\degree$K]', 
           fontdict=(fontThesis)) 
plt.ylabel("Tin conductivity",fontdict=(fontThesis)) 
#plt.axvline(x=277.13838, c='green')
plt.axvline(x=505.1, c='k', label = 'Melting point')
plt.legend(loc="upper left", prop={'size': 12})

plt.savefig('../images/2021-06-25-TinK.pdf', bbox_inches='tight', dpi=300)
#plt.savefig('../images/2021-06-25-TinK.png', bbox_inches='tight', dpi=300)