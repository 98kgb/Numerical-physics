# -*- coding: utf-8 -*-
"""
THz Pulse Propagation using UPPE
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fft2, ifft2, fftfreq

class THzPulse:
    def __init__(self, t, x, y, z_max, N_z, I = 1e13, w0 = 0.015, tau = 0.5e-12,
                 lambda_0=1.55e-6, n2=4.5e-18, chi_3 = 1e-20, alpha=0.5, omega_R=1.56e13, delta_R=5e11):
        """
        Initialize the THzPulse simulation parameters.

        Parameters:
        ===========================DOMAINS===========================
        t : array
            Time array (1D).
        x, y : array
            Spatial arrays for x and y directions (1D).
        z_max : float
            Maximum propagation distance (m).
        N_z : int
            Number of propagation steps.
        ===========================LIGHT SOURCE===========================
        I: float, optional
            Input pulse intensity (W/m^2)
        w0: float, optional
            Beam waist (m)
        tau: float, optional
            input pulse duration(s)
        lambda_0 : float, optional
            Central wavelength (m). Default is 1.55e-6.
        ===========================MATERIAL===========================
        n2 : float, optional
            Nonlinear refractive index (m^2/W). Default is 4.5e-18.
        alpha : float, optional
            Fraction of Kerr effect. Default is 0.5.
        omega_R : float, optional
            Raman resonance angular frequency. Default is 1.56e13 (for Silicon).
        delta_R : float, optional
            Raman linewidth. Default is 5e11 (for Silicon).
        """
        self.lambda_0 = lambda_0
        self.n2 = n2
        self.chi_3 = chi_3
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
        self.c = 3e8 # (m/s)
        self.epsilon_0 = 8.854e-12 # (F/m)
        self.k_0 = 2 * np.pi / self.lambda_0
        self.omega = 2 * np.pi * fftfreq(len(self.t), self.dt)
        
        # Transverse spatial frequency arrays
        self.k_x = 2 * np.pi * fftfreq(len(self.x), self.dx)
        self.k_y = 2 * np.pi * fftfreq(len(self.y), self.dy)
        self.k_x, self.k_y = np.meshgrid(self.k_x, self.k_y)
        self.k_perp2 = self.k_x**2 + self.k_y**2

        # Initial electric field (Gaussian pulse in time and space)
        T, X, Y = np.meshgrid(self.t, self.x, self.y, indexing="ij")
        
        E_amp = np.sqrt(2 * I / self.epsilon_0 /self.c) # initial amplitude of pulse
        w0 = w0
        tau = tau
        
        self.E_0 = E_amp * np.exp(-T**2 / (2 * (tau)**2)) * np.exp(-(X**2 + Y**2) / (2 * (w0)**2))

    
    def linear_operator(self):
        """
        Calculate the linear propagation operator in the frequency domain.
        """
        k2_omega = (self.omega / self.c)**2
        k2_omega = k2_omega[:, None, None]  # Expand dimensions to (N_t, 1, 1)
        k_perp2 = self.k_perp2[None, :, :]  # Expand dimensions to (1, N_x, N_y)

        diff = k2_omega - k_perp2
        diff = np.where(diff > 0, diff, 0)  # eliminate negative values

        return -1j * np.sqrt(diff)
        

    def kerr_polarization(self, E):
        """
        Calculate Kerr effect polarization in the time domain.
        """
        intensity = np.abs(E)**2
        return self.n2  * intensity * E * self.chi_3

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
        return ifft(raman_polarization_w, axis=0) * self.chi_3

    def solve(self):
        """
        Solve the THz pulse propagation using the Split-Step Fourier Method.
        """
        E = self.E_0.copy()
        lin_op = self.linear_operator()
        
        P_record = np.zeros(self.N_z)
        
        for ii in tqdm(range(self.N_z), desc = 'UPPE solving...'):
            # Step 1: Linear propagation in the frequency domain
            E_w = fft2(E, axes=(1, 2))  # Spatial Fourier Transform
            E_w *= np.exp(lin_op * self.dz)
            E = ifft2(E_w, axes=(1, 2))

            # Step 2: Calculate nonlinear polarization
            P_kerr = self.kerr_polarization(E)
            P_raman = self.raman_polarization(E)
            
            P_total = P_kerr * self.alpha + P_raman * (1 - self.alpha)
            P_record[ii] = np.max(np.abs(P_total))

            # Step 3: Apply nonlinear change in the time domain
            E += 1j * np.clip(P_total * self.dz, -1e10, 1e10)  # 변화 제한
            
        self.P_record = P_record
        
        return E
        
    def calculate_time_center(self, E_time, t):
        intensity = np.abs(E_time)**2
        return np.sum(t * intensity) / np.sum(intensity)
    

