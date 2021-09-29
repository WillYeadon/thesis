# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# Ti has a molar mass of 47.87 g/mol

def asm(a, b, c, T):
    # DOI: 10.1361/asmhba0005240
    # 4.184 J/cal
    cp_asm = (4.184/0.04787)*(a + b*T
                      + c/(T**2))
    return cp_asm

# asm alpha
T_range_asm_alpha = np.arange(298, 1155, 1)
a_asm_alpha = 5.28
b_asm_alpha  = 0.0024
c_asm_alpha  = 0
cp_asm_alpha = asm(a_asm_alpha, b_asm_alpha, c_asm_alpha, T_range_asm_alpha)

# asm beta
T_range_asm_beta = np.arange(1155, 1350, 1)
a_asm_beta = 4.74
b_asm_beta = 0.0019
c_asm_beta = 0

cp_asm1 = asm(a_asm_beta, b_asm_beta, c_asm_beta, T_range_asm_beta)

def NIST(a, b, c, d, e, T):
    T_local = T/1000
    cp_nist = (1/0.04787)*(a + b*T_local + c*T_local**2
            + d*T_local**3 + e/(T_local**2))
    return cp_nist

# NIST 298 - 700 alpha
T_range_nist_alpha_0 = np.arange(298, 700, 1)
a_nist_alpha_0 = 22.61942	
b_nist_alpha_0 = 18.98795
c_nist_alpha_0 = -18.18735
d_nist_alpha_0 = 7.080792	
e_nist_alpha_0 = -0.143457
cp_nist_alpha_0 = NIST(a_nist_alpha_0, b_nist_alpha_0, c_nist_alpha_0,
                       d_nist_alpha_0, e_nist_alpha_0, T_range_nist_alpha_0)

# NIST 700 - 1700 beta
T_range_nist_alpha_1 = np.arange(700, 1700, 1)
a_nist_alpha_1 = 44.37174	
b_nist_alpha_1 = -44.09225
c_nist_alpha_1 = 31.70602
d_nist_alpha_1 = 0.052209	
e_nist_alpha_1 = 0.036168
cp_nist_alpha_1 = NIST(a_nist_alpha_1, b_nist_alpha_1, c_nist_alpha_1, 
                       d_nist_alpha_1, e_nist_alpha_1, T_range_nist_alpha_1)

# NIST 298 - 1939 beta phase Ti
T_range_nist_beta = np.arange(298, 1939, 1)
a_nist_beta = 23.05660
b_nist_beta = 5.541331
c_nist_beta = -2.055881
d_nist_beta = 1.611745
e_nist_beta = -0.056075
cp_nist_beta = NIST(a_nist_beta, b_nist_beta, c_nist_beta, 
                    d_nist_beta, e_nist_beta, T_range_nist_beta)

fontThesis = {'fontname' : 'Arial', 'fontsize' : 12}

plt.plot(T_range_nist_alpha_0, cp_nist_alpha_0, 'r--', linewidth = 2.0, 
         label = "NIST " + r'$Ti_{(\alpha)}$')
plt.plot(T_range_nist_alpha_1, cp_nist_alpha_1, 'r--', linewidth = 2.0)
plt.plot(T_range_nist_beta, cp_nist_beta, 'g--', linewidth = 2.0, 
         label = "NIST " + r'$Ti_{(\beta)}$')
plt.plot(T_range_asm_alpha, cp_asm_alpha, 'b', linewidth = 2.0, 
         label = "ASM " + r'$Ti_{(\alpha)}$')
plt.plot(T_range_asm_beta, cp_asm1, 'm', linewidth = 2.0, 
         label = "ASM " + r'$Ti_{(\beta)}$')
plt.tick_params('both', labelsize=10)
plt.xlabel("Temperature " + "[K]", #+ r'[$\degree$K]', 
           fontdict=(fontThesis)) 
plt.ylabel("Specific Heat Capacity [J/kg $\cdot$ K]",
           #"Specific Heat Capacity J · kg⁻¹ K⁻¹", 
           #r'Specific Heat Capacity [ $J \cdot kg^{-1} \cdot K^{-1}$ ]', 
           fontdict=(fontThesis)) 
plt.legend(loc="upper left", prop={'size': 12})

plt.savefig('../images/cp-nist-asm.png', bbox_inches='tight', dpi=300)