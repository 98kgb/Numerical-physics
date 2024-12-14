"""
This script visualizes THz pulse propagation, analyzing the effects of different 

GVD coefficients (beta2) and nonlinear indices (n2).

Results are loaded from previously saved data, showing time/frequency domain evolution

and peak field/phase shift trends.

Generates detailed plots to compare propagation dynamics for various parameter sweeps.
"""


import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft, ifft

# Add current directory to system path to enable module imports
dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)

#%% visualization for diffraction calculation
I0 = 1e8  # Peak intensity [W/m^2]
w0 = 0.5e-4  # Beam waist [m]
tau = 0.01e-12  # Temporal width (short pulse, in seconds)
lambda_0 = 1.55e-6  # Central wavelength [m]
Lx = 1e-3  # Spatial simulation window [m]
T = 0.1e-12  # Temporal simulation window [s]
N = 1024 # Grid size (for spatial and time resolution)
Lz = 0.02  # Propagation length [m]
dz = 0.0005  # Step size [m]
z_steps = int(Lz/dz)
t = np.linspace(-T/2, T/2, N)
x = np.linspace(-Lx/2, Lx/2, N)
n0 = 3.48 # Si refractive index
n2 = 0
beta2 = 0

cut_point = 1 # visualization parameter (cutting inphysical range)

# loading data
E0 = np.load(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy')    
spatial_map = np.load(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy')
temporal_map = np.load(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy')
phase_shift = np.load(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_phase.npy')

# reconstruct gaussian pulse
E0 = E0[:, N//2]
E_con = fft(E0)
E_con[int(N/2):] = 0+1j*0
E_con = ifft(E_con)

Ef = temporal_map[int(z_steps//cut_point)-1,:]
Ef = np.real(Ef)
Ef_con = fft(Ef)
Ef_con[int(N/2):] = 0 + 1j*0
Ef_con = ifft(Ef_con)

map_con = np.zeros(np.shape(temporal_map))

for kk in range(temporal_map.shape[0]):
    E_temp = temporal_map[kk, :]
    E_temp = np.real(E_temp)
    E_temp = fft(E_temp)
    E_temp[int(N/2):] = 0 + 1j*0
    E_temp = ifft(E_temp)
    map_con[kk,:] = 2*abs(E_temp)

# visualizing
_, ax = plt.subplots(2, 2, figsize = (12,10))

window_ratio = [0,1]
t_plot = t[int(len(t)*window_ratio[0]):int(len(t)*window_ratio[1])]*1e12
f = np.fft.fftfreq(N, T/N)*1e-12

fontsize = 14
ax[0,0].plot(t*1e12, E0, 'r')
ax[0,0].set_xlabel('Time [ps]', fontsize = fontsize)
ax[0,0].set_ylabel('E-field [V/m]', fontsize = fontsize)

E_plot = fft(E0)
ax[0,1].plot(f, E_plot, 'r')
ax[0,1].set_xlabel('Frequency [THz]', fontsize = fontsize)

E_plot[int(N/2):] = 0+1j*0

ax[1,0].plot(f, E_plot, 'r')
ax[1,0].set_xlabel('Frequency [THz]', fontsize = fontsize)
ax[1,0].set_ylabel('E-field [V/m]', fontsize = fontsize)

E_plot = ifft(E_plot)

E_plot = 2*np.abs(E_plot)
ax[1,1].plot(t*1e12, E_plot, 'r')
ax[1,1].set_xlabel('Time [ps]', fontsize = fontsize)

ax[0,0].tick_params(axis='both', direction='in', length=5, width=1, labelsize=12)
ax[0,1].tick_params(axis='both', direction='in', length=5, width=1, labelsize=12)
ax[1,0].tick_params(axis='both', direction='in', length=5, width=1, labelsize=12)
ax[1,1].tick_params(axis='both', direction='in', length=5, width=1, labelsize=12)

plt.tight_layout()

_, ax = plt.subplots(1,3, figsize = (24,8))

# Plot pulse
ax[0].plot(t_plot, E0[int(len(t)*window_ratio[0]):int(len(t)*window_ratio[1])], 'r', linestyle='--')
ax[0].plot(t_plot, 2*np.abs(E_con[int(len(t)*window_ratio[0]):int(len(t)*window_ratio[1])]), 'r', label="Initial Pulse")    
ax[0].plot(t_plot, Ef[int(len(t)*window_ratio[0]):int(len(t)*window_ratio[1])], 'b', linestyle='--')
ax[0].plot(t_plot, 2*np.abs(Ef_con[int(len(t)*window_ratio[0]):int(len(t)*window_ratio[1])]), 'b', label="Final Pulse")
ax[0].set_xlabel("Time [s]", fontsize = 18)
ax[0].set_ylabel(r"$E$", fontsize = 18)
ax[0].legend(fontsize = 16)

# Calculate global vmin and vmax across both datasets
vmin = min(abs(spatial_map).min(), abs(map_con).min())
vmax = max(abs(spatial_map).max(), abs(map_con).max())

# Spatial evolution plot
ax[1].imshow(abs(spatial_map[:int(z_steps//cut_point), :]), extent=[x.min()*1e3, x.max()*1e3, Lz/cut_point*1e3, 0],
                 aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
ax[1].set_xlabel(r"$X$ direction [mm]", fontsize=18)
ax[1].set_ylabel(r"Propagation Distance [mm]", fontsize=18)

# temporal evolution
ax[2].imshow(abs(map_con[:int(z_steps//cut_point), :]), extent=[t.min()*1e12, t.max()*1e12, Lz/cut_point*1e3, 0],
                  aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
ax[2].set_xlabel("Time [ps]", fontsize=18)
ax[2].set_ylabel("Propagation Distance [mm]", fontsize=18)

#%% Visualization of beta sweep
I0 = 1e8  # Peak intensity [W/m^2]
w0 = 0.5e-4  # Beam waist [m]
tau = 0.01e-12  # Temporal width (short pulse, in seconds)
lambda_0 = 0.6e-6  # Central wavelength [m]
Lx = 1e-3  # Spatial simulation window [m]
T = 0.1e-12  # Temporal simulation window [s]
N = 1024 # Grid size (for spatial and time resolution)
Lz = 0.02  # Propagation length [m]
dz = 0.0005  # Step size [m]
z_steps = int(Lz/dz)
t = np.linspace(-T/2, T/2, N)
x = np.linspace(-Lx/2, Lx/2, N)
n0 = 3.48 # Si refractive index
n2 = 0  # Nonlinear refractive index [m^2/W]

beta2_candi = np.array([1e-32, 1e-31, 1e-30, 1e-29])[::-1]

cut_point = 1

_, ax = plt.subplots(3,len(beta2_candi), figsize = (6*len(beta2_candi),18))

# initializing arrays
E_peak = []
Phase_shift_mag = []

for ii in range(len(beta2_candi)):
        
    beta2 = beta2_candi[ii]
    
    E0 = np.load(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy')    
    spatial_map = np.load(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy')
    temporal_map = np.load(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy')
    phase_shift = np.load(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_phase.npy')
    
    Phase_shift_mag.append(sum(phase_shift[0,:]))
    
    # Reconstruct pulse
    E0 = E0[:, N//2]
    E_con = fft(E0)
    E_con[int(N/2):] = 0+1j*0
    E_con = ifft(E_con)
    ax[0,ii].plot(t, E0, 'r', linestyle='--')
    ax[0,ii].plot(t, 2*np.abs(E_con), 'r', label="Initial Pulse")    
    
    Ef = temporal_map[int(z_steps//cut_point)-1,:]
    Ef = np.real(Ef)
    Ef_con = fft(Ef)
    Ef_con[int(N/2):] = 0 + 1j*0
    Ef_con = ifft(Ef_con)
    
    map_con = np.zeros(np.shape(temporal_map))

    for kk in range(temporal_map.shape[0]):
        E_temp = temporal_map[kk, :]
        E_temp = np.real(E_temp)
        E_temp = fft(E_temp)
        E_temp[int(N/2):] = 0 + 1j*0
        E_temp = ifft(E_temp)
        map_con[kk,:] = 2*abs(E_temp)
        
    ax[0,ii].plot(t, Ef, 'b', linestyle='--')
    ax[0,ii].plot(t, 2*np.abs(Ef_con), 'b', label="Final Pulse")
    ax[0,ii].set_title(rf'$\beta_2$: {beta2}', fontsize = 20)
    ax[0,ii].set_xlabel("Time [s]", fontsize=20)
    ax[0,ii].set_ylabel("E-field", fontsize=20)
    ax[0,ii].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax[0,ii].yaxis.get_offset_text().set_fontsize(16)
    ax[0,ii].legend()
    E_peak.append(max(2*np.abs(Ef_con)))
    
        
    # Calculate global vmin and vmax across both datasets
    vmin = min(abs(spatial_map).min(), abs(map_con).min())
    vmax = max(abs(spatial_map).max(), abs(map_con).max())
    
    # Spatial evolution plot
    ax[1, ii].imshow(abs(spatial_map[:int(z_steps//cut_point), :]), extent=[x.min()*1e3, x.max()*1e3, Lz/cut_point*1e3, 0],
                     aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
    ax[1, ii].set_xlabel("X direction [mm]", fontsize=20)
    ax[1, ii].set_ylabel("Propagation Distance [mm]", fontsize=20)
    
    # Temporal evolution plot
    ax[2, ii].imshow(abs(map_con[:int(z_steps//cut_point), :]), extent=[t.min()*1e12, t.max()*1e12, Lz/cut_point*1e3, 0],
                      aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
    ax[2, ii].set_xlabel("Time [ps]", fontsize=18)
    ax[2, ii].set_ylabel("Propagation Distance [mm]", fontsize=18)

for axis in ax.flatten():
    axis.tick_params(axis='both', which='major', labelsize=18)
    
plt.suptitle(rf"THz Pulse Propagation Analysis for Varying $\beta_2$ with $I_0: {I0:.1e}$  $n_2$: {n2}", fontsize=24, y=1.02)
plt.tight_layout()
plt.show()

plt.figure(1)
plt.semilogx(beta2_candi,E_peak, '--b', marker = 'o') if all(beta2_candi>0) else plt.semilogx(abs(beta2_candi),E_peak, '--b', marker = 'o')
plt.xlabel(r'$\beta_2$') if all(beta2_candi>0) else plt.xlabel(r'|$\beta_2$|')
plt.ylabel('Peak $E$ value')
plt.title(rf'$n_2$: {n2}')
plt.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

plt.figure(2)
plt.semilogx(beta2_candi,Phase_shift_mag, '--r', marker = 'o') if all(beta2_candi>0) else plt.semilogx(abs(beta2_candi),E_peak, '--r', marker = 'o')
plt.xlabel(r'$\beta_2$') if all(beta2_candi>0) else plt.xlabel(r'|$\beta_2$|')
plt.ylabel('Maximum Accumulated Phase shift')
plt.title(rf'$n_2$: {n2}')
plt.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

#%% visualization for n2 sweep
cut_point = 1
n2_candi = np.array([1e-13, 1e-14, 1e-15])[::-1]
beta2 = 1e-29

# initializing arrays
E_peak, phase = [], []

_, ax = plt.subplots(3,len(n2_candi), figsize = (24,18))
for ii in range(len(n2_candi)):
    
    n2 = n2_candi[ii]
    # Load data
    E0 = np.load(f'{dir_path}\\result\\sweep_n2\\I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy')
    spatial_map = np.load(f'{dir_path}\\result\\sweep_n2\\I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy')
    temporal_map = np.load(f'{dir_path}\\result\\sweep_n2\\I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy')
    phase_shift = np.load(f'{dir_path}\\result\\sweep_n2\\I0_{I0}_n2_{n2}_beta2_{beta2}_phase.npy')
    
    phase.append(np.max(phase_shift,axis=1))
    
    # Reconstruct pulse
    E0 = E0[:, N//2]
    E_con = fft(E0)
    E_con[int(N/2):] = 0+1j*0
    E_con = ifft(E_con)
    
    Ef = temporal_map[int(z_steps//cut_point)-1,:]
    Ef = np.real(Ef)
    Ef_con = fft(Ef)
    Ef_con[int(N/2):] = 0 + 1j*0
    Ef_con = ifft(Ef_con)
    
    map_con = np.zeros(np.shape(temporal_map))

    for kk in range(temporal_map.shape[0]):
        E_temp = temporal_map[kk, :]
        E_temp = np.real(E_temp)
        E_temp = fft(E_temp)
        E_temp[int(N/2):] = 0 + 1j*0
        E_temp = ifft(E_temp)
        map_con[kk,:] = 2*abs(E_temp)
    
    
    # Plot pulse
    ax[0,ii].plot(t, E0, 'r', linestyle='--')
    ax[0,ii].plot(t, 2*np.abs(E_con), 'r', label="Initial Pulse")    
    ax[0,ii].plot(t, Ef, 'b', linestyle='--')
    ax[0,ii].plot(t, 2*np.abs(Ef_con), 'b', label="Final Pulse")
    ax[0,ii].set_title(f'n2: {n2}', fontsize = 20)
    ax[0,ii].set_xlabel("Time [s]", fontsize=20)
    ax[0,ii].set_ylabel("E-field", fontsize=20)
    ax[0,ii].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax[0,ii].yaxis.get_offset_text().set_fontsize(16)
    ax[0,ii].legend()
    
    E_peak.append(max(2*np.abs(Ef_con)))
        
    # Calculate global vmin and vmax across both datasets
    vmin = min(abs(spatial_map).min(), abs(map_con).min())
    vmax = max(abs(spatial_map).max(), abs(map_con).max())
    
    # Spatial evolution plot
    ax[1, ii].imshow(abs(spatial_map[:int(z_steps//cut_point), :]), extent=[x.min()*1e3, x.max()*1e3, Lz/cut_point*1e3, 0],
                     aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
    ax[1, ii].set_xlabel("X direction [mm]", fontsize=20)
    ax[1, ii].set_ylabel("Propagation Distance [mm]", fontsize=20)
    
    # Temporal evolution plot
    ax[2, ii].imshow(abs(map_con[:int(z_steps//cut_point), :]), extent=[t.min()*1e12, t.max()*1e12, Lz/cut_point*1e3, 0],
                      aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
    ax[2, ii].set_xlabel("Time [ps]", fontsize=18)
    ax[2, ii].set_ylabel("Propagation Distance [mm]", fontsize=18)

for axis in ax.flatten():
    axis.tick_params(axis='both', which='major', labelsize=18)
    
plt.suptitle(rf"THz Pulse Propagation Analysis for Varying n2 with $I_0: {I0:.1e}$  $\beta_2$: {beta2}", fontsize=24, y=1.02)
plt.tight_layout()
plt.show()

plt.figure(1)
plt.semilogx(n2_candi,E_peak, '--b', marker = 'o')
plt.xlabel(r'$n_2$')
plt.ylabel('Peak $E$ value')
plt.title(rf'$\beta_2$: {beta2}')
plt.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

plt.figure(2)
phase = np.array(phase)
plt.plot(n2_candi,phase[:,0], '--b', marker = 'o', label = 'GVD')
plt.plot(n2_candi,phase[:,1], '--r', marker = 'o', label = 'Kerr')
plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$n_2$') 
plt.ylabel('Maximum Accumulated Phase shift')
plt.title(rf'$beta_2$: {beta2}')
plt.legend()
