# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:26:11 2024

@author: Gibaek
"""

import os
import sys
import ast

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


t = np.linspace(-5e-12, 5e-12, 256)  # Time array from (0ps to 15ps)
x = np.linspace(-1e-2, 1e-2, 128)  # x spatial array (2cm total)
y = np.linspace(-1e-2, 1e-2, 128)  # y spatial array (2cm total)
z_max = 0.02  # Propagation distance (2cm)

n2_candi = [1e-12, 1e-13, 1e-14, 1e-15]
# N_z_candi = [100,200,300,400,500]
N_z_candi = [100, 200, 400, 500]

n2 = n2_candi[2]

_, ax = plt.subplots(2, 4, figsize=(20, 9))  # Create a 4x2 grid for all subplots

for ii in range(len(N_z_candi)):
    fig_idx = ii  # Each N_z_candi gets one row for two subplots
    N_z = N_z_candi[ii]

    # Read CSV and parse complex numbers
    df_tot = pd.read_csv(f'./result/n2_{n2}_Nz_{N_z}.csv')
    col = df_tot.columns
    for col in col:
        if col != 'Unnamed: 0':
            df_tot[col] = df_tot[col].apply(ast.literal_eval)

    E_linear, E_kerr = df_tot['Linear E'].values[:int(len(t))], df_tot['Kerr E'].values[:int(len(t))]
    P_linear, P_kerr = df_tot['Linear P'].values[:N_z], df_tot['Kerr P'].values[:N_z]
    t_linear, t_kerr = abs(df_tot['Linear t ave'].values[0]), abs(df_tot['Kerr t ave'].values[0])

    # First subplot: Intensity comparison
    ax[0,fig_idx].plot(t * 1e12, np.abs(E_linear)**2, 'r', label="Linear")
    ax[0,fig_idx].plot(t * 1e12, np.abs(E_kerr)**2, 'b', label="Kerr")
    ax[0,fig_idx].set_xlabel("Time (ps)", fontsize = 15)
    ax[0,fig_idx].set_ylabel("Intensity", fontsize = 15)
    ax[0,fig_idx].set_title(f"Pulse Intensity (step: {N_z})", fontsize = 18)
    ax[0,fig_idx].legend()
    ax[0,fig_idx].grid(True)

    # Second subplot: Polarization magnitude
    z_step = np.linspace(0, z_max, N_z)
    ax[1,fig_idx].plot(z_step * 1e3, np.abs(P_linear), 'r', label="Linear")
    ax[1,fig_idx].plot(z_step * 1e3, np.abs(P_kerr), 'b', label="Kerr")
    ax[1,fig_idx].set_xlabel("Propagation Distance (mm)", fontsize = 15)
    ax[1,fig_idx].set_ylabel("Polarization Magnitude", fontsize = 15)
    ax[1,fig_idx].set_title(f"Kerr Pol Magnitude\n (N_z: {N_z}, Time Delay: {t_kerr*1e12:.6f} ps)"
                            , fontsize = 18)
    ax[1,fig_idx].legend()
    ax[1,fig_idx].grid(True)

plt.tight_layout()
plt.show()
