# -*- coding: utf-8 -*-
"""
This code search proper parameter to observe the time delay which originated form Kerr effect.

The candidate paramters are n2 (nonlinearity index), and N_z (resolution of z direction).

@author: Gibaek
"""


import numpy as np
import pandas as pd
from main import THzPulse

# sweep parameters
n2_candi = [1e-12, 1e-13, 1e-14, 1e-15]
N_z_candi = [100,200,300,400,500]

# Define simulation parameters
t = np.linspace(-5e-12, 5e-12, 256)  # Time array from (0ps to 15ps)
x = np.linspace(-1e-2, 1e-2, 128)  # x spatial array (2cm total)
y = np.linspace(-1e-2, 1e-2, 128)  # y spatial array (2cm total)
z_max = 0.02  # Propagation distance (2cm)

x_idx, y_idx = 63, 63

columns = ['Linear E', 'Kerr E', 'Linear P', 'Kerr P', 'Linear t ave', 'Kerr t ave']

for ii in range(len(n2_candi)):
    for jj in range(len(N_z_candi)):
        
        n2 = n2_candi[ii]
        N_z = N_z_candi[jj]
        
        record = np.zeros([max(len(t), N_z), 6], dtype=np.complex128)

        linear = THzPulse(t, x, y, z_max, N_z, I = 1e13, w0 = 0.015, tau = 0.5e-12,
                     lambda_0 = 1.55e-6, n2 = 0, chi_3 = 1,  alpha = 1, omega_R = 1.56e13, delta_R = 5e11)
        E_linear = linear.solve()
        time_center_linear = linear.calculate_time_center(E_linear[:,64,64], t)
        
        kerr = THzPulse(t, x, y, z_max, N_z, I = 1e13, w0 = 0.015, tau = 0.5e-12,
                     lambda_0 = 1.55e-6, n2 = n2, chi_3 = 1, alpha = 1, omega_R = 1.56e13, delta_R = 5e11)
        E_kerr = kerr.solve()
        time_center_kerr = kerr.calculate_time_center(E_kerr[:,64,64], t)
        
        record[:len(t),0] = E_linear[:, x_idx, y_idx]
        record[:len(t),1] = E_kerr[:, x_idx, y_idx]
        record[:N_z,2] = linear.P_record[:N_z]
        record[:N_z,3] = kerr.P_record[:N_z]
        record[0,4] = time_center_linear
        record[0,5] = time_center_kerr
        
        df_tot = pd.DataFrame(record, columns = columns)
        df_tot.to_csv(f'./result/n2_{n2}_Nz_{N_z}.csv')
        
        print(f"Linear Time Center: {time_center_linear * 1e12:.8f} ps")
        print(f"Kerr Time Center: {time_center_kerr * 1e12:.8f} ps")

        
