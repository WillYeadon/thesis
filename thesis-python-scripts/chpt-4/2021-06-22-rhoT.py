# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

def rhoT(a, b, c, d, e, T):
    return (a + b*T - c*T**2 + d*T**3 - e*T**4)

def boussT(rho, beta, T):
    return rho*(1.0 - beta*T)

T_range = np.arange(-2.5, 27.5, 0.5)
K_range = T_range + 273.15
rhoT_values = rhoT(999.840281167108,
                   0.0673268037314653,
                   0.00894484552601798,
                   8.7846286650041**(-5), 
                   6.6213979262754**(-7),
                   T_range)

boussT_values = boussT(999.85, 0.000214, T_range)


# wolfram alpha gives max at 3.98838 C
# = 277.13838 K
fontThesis = {'fontname' : 'Arial', 'fontsize' : 12}

plt.axvline(x=273.15, c='k', linewidth=2.0, label = 'Melting point')
plt.axvline(x=277.13838, c='green', linewidth=2.0, label = 'Peak density')
plt.plot(K_range, boussT_values, 'b:', linewidth = 2.0,
         label = 'Boussinesq')
plt.plot(K_range, rhoT_values, 'r--', linewidth = 2.0,
         label = 'Function')
plt.tick_params('both', labelsize=10)
plt.xlabel("Temperature " + "[K]", #+ r'[$\degree$K]', 
           fontdict=(fontThesis)) 
plt.ylabel("Density [kg/m3]",fontdict=(fontThesis)) 
plt.legend(loc="upper right", prop={'size': 12})

plt.savefig('../images/2021-06-22-rhoT.pdf', bbox_inches='tight', dpi=300)
plt.savefig('../images/2021-06-22-rhoT.png', bbox_inches='tight', dpi=300)