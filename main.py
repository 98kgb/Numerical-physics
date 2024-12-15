"""
This Python script simulates the propagation of a terahertz (THz)

electromagnetic pulse in space-time, including both linear effects

(diffraction and group velocity dispersion, GVD) and nonlinear effects (Kerr effect).

It uses the Unidirectional Pulse Propagation Equation (UPPE) to model these effects.
"""


import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
from tqdm import tqdm

class THzProp:
    def __init__(self, w0, tau, lambda_0, Lx, T, N, dz, Lz, n0, n2, beta2, I0, chirp = -1e28, verbose = False):
        
        """
        Initialize the parameters for space-time Kerr effect propagation.

        Parameters:
        w0 (float): Beam waist [m]
        tau (float): Temporal width [s]
        lambda_0 (float): Wavelength [m]
        Lx (float): Spatial grid length [m]
        T (float): Temporal grid length [s]
        N (int): Grid size (spatial and temporal resolution)
        dz (float): Propagation step size [m]
        Lz (float): Total propagation distance [m]
        n0 (float): Linear refractive index
        n2 (float): Nonlinear refractive index [m^2/W]
        beta2 (float): Group velocity dispersion coefficient [s^2/m]
        I0 (float): Peak intensity [W/m^2]
        """
        self.verbose = verbose # True to print the calculation procedure
        
        # Physical constants
        self.epsilon_0 = 8.854e-12 # Vacuum permittivity [F/m]
        self.mu0 = 4 * np.pi * 1e-7  # Vacuum permeability [H/m]
        self.c = 3e8 # Speed of light [m/s]
        
        # Material parameters
        self.n0 = n0  # Linear refractive index
        self.n2 = n2  # Nonlinear refractive index
        self.beta2 = beta2  # Group velocity dispersion coefficient
        # Third-order nonlinear susceptibility derived from n2
        self.chi3 = (4 * self.n0 * self.epsilon_0 * self.c * self.n2) / 3
        
        # Pulse parameters
        self.I0 = I0  # Peak intensity
        self.w0 = w0  # Beam waist
        self.tau = tau  # Temporal width
        self.lambda_0 = lambda_0  # Central wavelength
        
        # Dimensional parameters
        self.k0 = 2 * np.pi / lambda_0  # Central wavevector
        self.Lx = Lx  # Spatial grid length
        self.T = T  # Temporal grid length
        self.Nx = N  # Number of spatial grid points
        self.Nt = N  # Number of temporal grid points
        self.dx = Lx / N  # Spatial grid resolution
        self.dt = T / N  # Temporal grid resolution
        self.dz = dz  # Propagation step size
        
        self.Lz = Lz  # Total propagation distance
        
        # Initialize spatial and temporal grids
        self.x = np.linspace(-Lx/2, Lx/2, N)  # Spatial grid
        self.t = np.linspace(-T/2, T/2, N)   # Temporal grid
        self.X, self.T = np.meshgrid(self.x, self.t)  # Space-time grid
        self.r = self.X  # Radial distance (1D spatial dimension assumed)

        # Initialize spatial k-space (kx) and angular frequency (omega)
        self.kx = fftfreq(N, self.dx) * 2 * np.pi  # Spatial frequency
        self.kz = np.sqrt(np.maximum(self.k0**2 - self.kx**2, 0) + 0j)  # Longitudinal wavevector component
        self.omega_0 = 2 * np.pi * self.c / self.lambda_0  # Central angular frequency
        self.omega = fftfreq(N, self.dt) * 2 * np.pi - self.omega_0  # Angular frequency deviation

        # Meshgrid for k-space and angular frequency
        self.KZ, self.OMEGA = np.meshgrid(self.kz, self.omega)
        
        # Define chirp paramters
        self.chirp = chirp
        
        # Define initial Gaussian beam in space and time
        self.E0 = np.sqrt(2 * self.I0 * self.epsilon_0 * self.c * self.n0) \
                  * np.exp(-self.r**2 / self.w0**2) \
                  * np.exp(-self.T**2 / self.tau**2) \
                  * np.cos(self.omega_0 * self.T + self.chirp *self.T **2)
        
        self.z_steps = int(Lz / self.dz)  # Number of steps in propagation
        
        # self.z_steps = 1
    def propagate(self):
        """
        Simulate the propagation of the pulse over the specified distance.
        
        Returns:
        intensity_map (numpy.ndarray): 2D array of intensity over propagation.
        E (numpy.ndarray): Final electric field after propagation.
        """
        
        # Copy the initial electric field
        E = self.E0.copy()
        E_list = []
        spatial_map = []  # Store spatial intensity profiles for visualization
        temporal_map = []  # Store temporal intensity profiles for visualization
        phase_shift = np.zeros([2, self.z_steps]) # Store phase shift profile.
        
        # Loop through each propagation step
        for ii in tqdm(range(self.z_steps), desc = 'Solving UPPE...'):
            
            # Linear propagation
            E = fft2(E) # Fourier transform of the electric field
            E *= np.exp(-1j * self.KZ * self.dz) # Apply diffraction effect
            E *= np.exp(-1j * (self.OMEGA * self.n0 / self.c) * self.dz) # Apply refractive index phase shift
            E *= np.exp(-1j * (self.beta2 * self.OMEGA**2) * self.dz) # Apply group velocity dispersion (GVD)
            E = ifft2(E) # Inverse Fourier transform back to space-time domain
            
            # Kerr effect calculation
            d_phase, dE = self.Kerr(E)

            # Update electric field with Kerr effect polarization and phase
            E *= np.exp(1j * d_phase)
            E += dE * self.dz
            
            # Store the intensity profile at the central time and space slice
            spatial_map.append(E[self.Nt // 2, :])
            temporal_map.append(E[:, self.Nx//2])
            
            # Store phase shift information
            phase_shift[0,ii] = np.max((abs(self.beta2 * self.OMEGA**2)) * self.dz)
            phase_shift[1,ii] = np.max(abs(d_phase))
            
            # Calculate and store total energy at this step
            total_energy = np.sum(np.abs(E)**2) * self.dx * self.dt
            E_list.append(total_energy)
            
            # Print calculation procedure if verbose is True
            if self.verbose:
                print('\nGVD:', np.max((abs(self.beta2 * self.OMEGA**2)) * self.dz))
                print('Kerr:', np.max(abs(d_phase)),'\n')
                print(f'E max: {np.max(E):.1e}', )
                
        return np.array(spatial_map), np.array(temporal_map), E_list, phase_shift
        
    def Kerr(self, E):
        # Compute intensity and clip it to mitigate overflow.
        intensity = np.abs(E)**2
        clip_max = self.I0 * 2
        intensity = np.clip(intensity, 0, clip_max)
        
        # Kerr effect
        d_phase = self.k0 * self.n2 * intensity * self.dz
        
        P_kerr = self.epsilon_0 * self.chi3 * E * intensity
        dE = 1j * (self.mu0 * self.OMEGA**2 / np.sqrt(self.KZ**2 + 0j)) * P_kerr
        
        return d_phase, dE
        
# Main program for testing
if __name__ == '__main__':
    
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
    n0 = 3.48 # Si refractive index
    
    chirp = 1e28 # Define chirp paramters
    n2 = 0 # Nonlinear refractive index
    beta2 = -1e-28 # GVD coefficient
    
    # Define model
    for ii in range(2):
        chirp = chirp if ii == 0 else -chirp
        model = THzProp(w0, tau, lambda_0, Lx, T, N, dz, Lz, n0, n2, beta2, I0, chirp)
        
        # Propagate
        spatial_map, temporal_map, E_list, phase_shift = model.propagate()
    
        # visualization
        t = model.t
        x = model.x
        
        # Reconstruct Gaussian shape
        E0 = temporal_map[0, :]
        E_con = fft(E0)
        E_con[int(N/2):] = 0+1j*0
        E_con = ifft(E_con)
        
        Ef = temporal_map[-1,:]
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
        
        # plotting part
        _, ax = plt.subplots(1,3, figsize = (24,8))
        fontsize = 24
        ax[0].plot(t, E0, 'r', linestyle='--')
        ax[0].plot(t, 2*np.abs(E_con), 'r', label="Initial Pulse")    
        ax[0].plot(t, Ef, 'b', linestyle='--')
        ax[0].plot(t, 2*np.abs(Ef_con), 'b', label="Final Pulse")
        ax[0].set_xlabel("Time [s]", fontsize = fontsize)
        ax[0].set_ylabel("Amplitude", fontsize = fontsize)
        
        # Calculate global vmin and vmax across both datasets
        vmin = min(abs(spatial_map).min(), abs(map_con).min())
        vmax = max(abs(spatial_map).max(), abs(map_con).max())
        
        # Spatial evolution plot
        ax[1].imshow(abs(spatial_map), extent=[x.min(), x.max(), Lz, 0], aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
        ax[1].set_xlabel("X direction [m]", fontsize = fontsize)
        ax[1].set_ylabel("Propagation Distance [m]", fontsize = fontsize)
        
        # Temporal evolution plot
        ax[2].imshow(abs(map_con), extent=[t.min()*1e12, t.max()*1e12, Lz, 0], aspect='auto', cmap='jet', vmin=vmin, vmax=vmax)
        ax[2].set_xlabel("Time [ps]", fontsize=18)
        ax[2].set_ylabel("Propagation Distance [m]", fontsize = fontsize)
        
        # set tick parameters
        labelsize = 18
        ax[0].tick_params(axis='both', direction='in', length=5, width=1, labelsize=labelsize)
        ax[1].tick_params(axis='both', direction='in', length=5, width=1, labelsize=labelsize)
        ax[2].tick_params(axis='both', direction='in', length=5, width=1, labelsize=labelsize)
        
        plt.suptitle(rf"THz Pulse Propagation Analysis for $I_0: {I0:.1e}$  $\beta_2$: {beta2}  $n_2: {n2}$  Chirp: {chirp}", fontsize=24, y=1.02)
        plt.tight_layout()
        plt.show()
    
