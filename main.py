# -*- coding: utf-8 -*-
"""
THz Pulse Propagation using UPPE
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fft2, ifft2, fftfreq

class THzPulse:
    def __init__(self, t, x, y, z_max, N_z, lambda_0=1.55e-6, n2=1e-20, alpha=0.5, omega_R=1e13, delta_R=1e12):
        """
        Initialize the THzPulse simulation parameters.

        Parameters:
        t : array
            Time array (1D).
        x, y : array
            Spatial arrays for x and y directions (1D).
        z_max : float
            Maximum propagation distance (m).
        N_z : int
            Number of propagation steps.
        lambda_0 : float, optional
            Central wavelength (m). Default is 1.55e-6.
        n2 : float, optional
            Nonlinear refractive index (m^2/W). Default is 1e-20.
        alpha : float, optional
            Fraction of Kerr effect. Default is 0.5.
        omega_R : float, optional
            Raman resonance angular frequency. Default is 1e13.
        delta_R : float, optional
            Raman linewidth. Default is 1e12.
        """
        self.lambda_0 = lambda_0
        self.n2 = n2
        self.alpha = alpha
        self.omega_R = omega_R
        self.delta_R = delta_R
        self.z_max = z_max
        self.N_z = N_z
        self.dz = z_max / N_z
        self.t = t
        self.dt = self.t[1] - self.t[0]
        self.x = x
        self.y = y
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.c = 3e8
        self.k_0 = 2 * np.pi / self.lambda_0
        self.omega = 2 * np.pi * fftfreq(len(self.t), self.dt)
        
        # Transverse spatial frequency arrays
        self.k_x = 2 * np.pi * fftfreq(len(self.x), self.dx)
        self.k_y = 2 * np.pi * fftfreq(len(self.y), self.dy)
        self.k_x, self.k_y = np.meshgrid(self.k_x, self.k_y)
        self.k_perp2 = self.k_x**2 + self.k_y**2

        # Initial electric field (Gaussian pulse in time and space)
        T, X, Y = np.meshgrid(self.t, self.x, self.y, indexing="ij")
        self.E_0 = np.exp(-T**2 / (2 * (50e-15)**2)) * np.exp(-(X**2 + Y**2) / (2 * (100e-6)**2))
        
    def linear_operator(self):
        """
        Calculate the linear propagation operator in the frequency domain.
        """
        k2_omega = (self.omega / self.c)**2
        k2_omega = k2_omega[:, None, None]  # Expand dimensions to (N_t, 1, 1)
        k_perp2 = self.k_perp2[None, :, :]  # Expand dimensions to (1, N_x, N_y)

        diff = k2_omega - k_perp2
        diff[diff < 0] = 0  # Set negative values to zero
    
        return -1j * np.sqrt(diff)
        

    def kerr_polarization(self, E):
        """
        Calculate Kerr effect polarization in the time domain.
        """
        intensity = np.abs(E)**2
        return self.alpha * intensity * E

    def raman_polarization(self, E):
        """
        Calculate Raman effect polarization in the frequency domain and convert to time domain.
        """
        intensity = np.abs(E)**2
        intensity_w = fft(intensity, axis=0)  # Fourier Transform along the time axis
        
        # Raman response in frequency domain
        g_R_w = self.omega_R**2 / (self.omega_R**2 + 2j * self.omega * self.delta_R - self.omega**2)
        raman_polarization_w = g_R_w[:, None, None] * intensity_w

        # Convert back to time domain
        return ifft(raman_polarization_w, axis=0) * (1 - self.alpha)

    def solve(self):
        """
        Solve the THz pulse propagation using the Split-Step Fourier Method.
        """
        E = self.E_0.copy()
        lin_op = self.linear_operator()

        for _ in tqdm(range(self.N_z), desc = 'UPPE solving...'):
            # Step 1: Linear propagation in the frequency domain
            E_w = fft2(E, axes=(1, 2))  # Spatial Fourier Transform
            E_w *= np.exp(lin_op * self.dz)
            E = ifft2(E_w, axes=(1, 2))

            # Step 2: Calculate nonlinear polarization
            P_kerr = self.kerr_polarization(E)
            P_raman = self.raman_polarization(E)
            P_total = P_kerr + P_raman

            # Step 3: Apply nonlinear change in the time domain
            E += 1j * P_total * self.dz

        return E

# Define simulation parameters
t = np.linspace(-5e-13, 5e-13, 256)  # Time array
x = np.linspace(-1e-3, 1e-3, 128)  # x spatial array
y = np.linspace(-1e-3, 1e-3, 128)  # y spatial array
z_max = 0.01  # Propagation distance (m)
N_z = 100  # Number of steps

# Create the THzPulse model and solve it
model = THzPulse(t, x, y, z_max, N_z)
E = model.solve()
#%% plot
x_pos = 64
y_pos = 65

plt.plot(t*10e12, np.abs(E[:,x_pos,y_pos])**2)

plt.xlabel('Time (ps)')

#%% Plot the final electric field intensity (time-integrated)
intensity = np.sum(np.abs(E)**2, axis=0)
plt.figure(figsize=(8, 6))
plt.imshow(intensity, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
plt.title("Time-Integrated Intensity")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.colorbar(label="Intensity")
plt.show()
