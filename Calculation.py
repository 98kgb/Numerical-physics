"""
THz Propagation Simulation: The code simulates the propagation of a THz electromagnetic pulse,

including linear (diffraction, GVD) and nonlinear (Kerr effect) effects using specified parameters.

Parameter Sweeps: It evaluates the effect of varying GVD coefficient (beta2) and nonlinear index (n2)

on pulse dynamics and saves results for each parameter set.

Result Storage: The code saves the initial pulse, propagated fields, and phase shifts as .npy files for analysis.
"""

import os
import sys

# Add current directory to system path to enable module imports
dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)

import numpy as np
from main import THzProp

#%% Diffraction calculation
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
n2 = 0  # Nonlinear refractive index [m^2/W]
beta2 = 0 # GVD coefficient

# Calculating
model = THzProp(w0, tau, lambda_0, Lx, T, N, dz, Lz, n0, n2, beta2, I0)
spatial_map, temporal_map, E_list, phase_shift = model.propagate()

# Saving result
np.save(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy', model.E0)
np.save(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E_list.npy', E_list)
np.save(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy', spatial_map)
np.save(f'{dir_path}\\result\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy', temporal_map)

#%% Beta dependance calculation
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

# Define candidates
beta2_candi = np.array([1e-32, 1e-31, 1e-30, 1e-29])

# Ensuring directory existnace
if not os.path.exists(f'{dir_path}\\result\\sweep_beta\\'):
    os.mkdir(f'{dir_path}\\result\\sweep_beta\\')

# Sweep through beta2
for ii in range(len(beta2_candi)):
        
    beta2 = beta2_candi[ii]
    
    # Calculating
    model = THzProp(w0, tau, lambda_0, Lx, T, N, dz, Lz, n0, n2, beta2, I0)    
    spatial_map, temporal_map, E_list, phase_shift = model.propagate()
    
    # Saving result
    np.save(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy', model.E0)
    np.save(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E_list.npy', E_list)
    np.save(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy', spatial_map)
    np.save(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy', temporal_map)
    np.save(f'{dir_path}\\result\\sweep_beta\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_phase.npy', phase_shift)
    
    
#%% n2 sweep
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
beta2 = 1e-29 # GVD coefficient

# Define candidate range
n2_candi = np.array([1e-12, 1e-13, 1e-14, 1e-15])

# Ensuring directory existnace
if not os.path.exists(f'{dir_path}\\result\\sweep_n2\\'):
    os.mkdir(f'{dir_path}\\result\\sweep_n2\\')

# Sweep through n2
for ii in range(len(n2_candi)):
    
    n2 = n2_candi[ii]
    
    # Calculating
    model = THzProp(w0, tau, lambda_0, Lx, T, N, dz, Lz, n0, n2, beta2, I0)
    spatial_map, temporal_map, E_list, phase_shift = model.propagate()
    
    #Saving
    np.save(f'{dir_path}\\result\\sweep_n2\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E0.npy', model.E0)
    np.save(f'{dir_path}\\result\\sweep_n2\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_E_list.npy', E_list)
    np.save(f'{dir_path}\\result\\sweep_n2\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_spat.npy', spatial_map)
    np.save(f'{dir_path}\\result\\sweep_n2\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_temp.npy', temporal_map)
    np.save(f'{dir_path}\\result\\sweep_n2\\lamda_{lambda_0}_I0_{I0}_n2_{n2}_beta2_{beta2}_phase.npy', phase_shift)


