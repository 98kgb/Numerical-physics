# -*- coding: utf-8 -*-
"""
This function is prepared for SCG, rewrite to ODE45 interaction picture method.
T: time vector
A: input field in time domain
w0: pump center freq in radius
gamma: nonlinear coef
Beta_interp, beta0, beta1: dispersion parameters to include all higher order dispersion contributions.
alpha: /m, attenuation coefficient
fr, RT: for Raman response
betas: do simulations based on betas, real geo are ignored
DV: dispersion variation function, it should work together with betas
tau_shock: gamma1/gamma0, to include self-steepening effect
flength: the length of waveguide
nsaves: numbers of step when doing simulation
nplot: numbers of frames to plot to save graphic demands.

@author: Gibaek
"""
import numpy as np
import math
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d, PchipInterpolator
eps = 2.2204e-16
# === define function to print ODE integrator status
# N_flag = 0 # count the flag number

def report(z, y, flag, N_flag, flength):
    status = 0
    
    # Check if flag is empty or None
    if flag is None or not flag:
        print(f"{z / flength * 100:05.1f} % complete")
        N_flag += 1

    # Update status every 256 steps
    if N_flag % 256 == 0:
        print(f"{z / flength * 100:05.1f} % complete")
    
    return status, N_flag  # Return N_flag to allow it to update externally


# Only 3 useful parameters: Beta_interp, beta1, beta2!
# Beta_interp: interpolate Beta = f(W) with W = w_interp;
# Beta1_interp: interpolate Beta1 = f(W) with W = w_interp;
# Beta2_interp: interpolate Beta2 = f(W) with W = w_interp;

# w0: pump frequency
# beta0: the closest interpolated value @ pump wl
# beta1: the closest interpolated value @ pump wl
# beta2: the closest interpolated value @ pump wl

def interp_beta(W,w_interp,Beta,Beta1,Beta2,lambda0):
    # Define speed of light constant
    c = 299792458
    
    # Define the interpolation functions using PCHIP
    Beta_interp_func = PchipInterpolator(W, Beta, extrapolate=True)
    Beta1_interp_func = PchipInterpolator(W, Beta1, extrapolate=True)
    Beta2_interp_func = PchipInterpolator(W, Beta2, extrapolate=True)
    # print("w_interp",w_interp.shape)
    # Evaluate the interpolators at w_interp
    Beta_interp = Beta_interp_func(w_interp)
    Beta1_interp = Beta1_interp_func(w_interp)
    Beta2_interp = Beta2_interp_func(w_interp)
    
    # Find the minimum absolute value in w_interp and its index
    pp = np.argmin(np.abs(w_interp))
    
    # Set values before index `pp` to zero
    Beta_interp[:pp] = 0
    
    # Flip section from pp+1 to 2*pp-1
    Beta_interp[pp+1:2*pp] = Beta_interp[pp+1:2*pp][::-1]
    
    # Calculate w0
    w0 = 2 * np.pi * c / lambda0
    
    # Find index of the closest value in w_interp to w0
    uw = np.argmin(np.abs(w_interp - w0))
    
    # Get beta values at the point w0
    beta0 = Beta_interp[uw]
    beta1 = Beta1_interp[uw]
    beta2 = Beta2_interp[uw]

    
    return [Beta_interp, beta0, beta1, beta2]




def gnlse_YY_ODE45_TaperAdapt(t, A_0, w0, w_interp, W_dat, d_us, us, L_tot, Betas_Dummy, Beta1s_Dummy, 
                              Beta2s_Dummy,A_effs_Dummy,n2_effs_Dummy,n_effs_Dummy, alpha, fr, RT, betas,DV):
    
    
    Z_z = []
    AT_z = []
    AW_z = []
    
    isHybridMode = 0

    ## Method Runge-Kutta technique and Interaction Picture
    AT_0 = A_0 # the initial AT_0 for NL Eq. is pulse A_0
    c = 299792458
    C = c
    lambda0=  2 * np.pi * c / w0 
    
    for ii in range(len(d_us)):
        # This method takes into account beta(w) without approximation (no Taylor expansion), it needs to define beta(w)
        # exact length of each segement, depend on us(i)
        flength = L_tot*d_us[ii]  
        nplot = np.round(d_us[ii]/0.002)# try to record at each 0.002(normalised propa. dis.)
        
        # slip step method: need to set the steps for dispersion
        nsaves = np.round(flength/(10*1e-6))# by default 10µm, one dispersion operator.
    
        if nsaves < 3: # for ODE, a minimum number of 3 is needed?
            nsaves = 3
        
        # Read pseudo simulation data, so called 'Dummy'.
        B = np.flip(Betas_Dummy[:,ii])
        B1 = np.flip(Beta1s_Dummy[:,ii])
        B2 = np.flip(Beta2s_Dummy[:,ii])
        
        A_eff = A_effs_Dummy[:,ii]
        n2_eff = n2_effs_Dummy[:,ii]
        n_eff = n_effs_Dummy[:,ii]

        # for gamma and tau_shock 
        n_interp=2**12
        n2 = 2.4*1e-19  # Nonlinear index (m^2/W)
        #n2 = 1.43*1e-17 # AlGaAs
        pw = (max(W_dat)-min(W_dat))/n_interp
        W_interp = np.flip(np.arange(min(W_dat),max(W_dat),pw)) # Need to be improved?
            
        A_eff_interp_func = interp1d(W_dat, A_eff, kind='cubic')
        n2_eff_interp_func = interp1d(W_dat, n2_eff, kind='cubic')
        n_eff_interp_func = interp1d(W_dat, n_eff, kind='cubic')

        # Perform the interpolation at the desired points W_interp
        A_eff_fit = A_eff_interp_func(W_interp)
        n2_eff_fit = n2_eff_interp_func(W_interp)
        n_eff_fit = n_eff_interp_func(W_interp)
        
        # @ pump wl:
        u = np.min(np.abs(W_interp - w0))   # Minimum value
        v = np.argmin(np.abs(W_interp - w0))  # Index of the minimum value
        
        Aeff1_0 = A_eff_fit[v]# find Aeff @ pump wl
        n2_eff = n2_eff_fit[v]# find n2 @ pump wl
        neff_0 = n_eff_fit[v] # find neff @ pump wl
            
        gamma = w0 * n2 / (C * Aeff1_0) # calculate nl coeff. @ pump wl
      
            
        der_Aeff1_0 = (A_eff_fit[v+1] - A_eff_fit[v-1]) / (W_interp[v+1] - W_interp[v-1])
        der_neff1_0 = (n_eff_fit[v+1] - n_eff_fit[v-1]) / (W_interp[v+1] - W_interp[v-1])
        der_n2eff1_0 = (n2_eff_fit[v+1] - n2_eff_fit[v-1]) / (W_interp[v+1] - W_interp[v-1])

        if isHybridMode == 1:
            
            print('isHybridMode: #d \n',isHybridMode)
            tau_shock = 1/w0-der_Aeff1_0/Aeff1_0-der_neff1_0/neff_0 + der_n2eff1_0/n2_eff
            gamma = w0*n2_eff/(C*Aeff1_0) # calculate nl coeff. @ pump wl from effective n2 curve
        
        else:
            tau_shock = 1/w0-der_Aeff1_0/Aeff1_0-der_neff1_0/neff_0
    
    
        # for Raman effect, which is negligible for SiN
        
        hr = 0
        fr = 0
        ###################################################
            
        ### Prepare for linear operator ###
        # interpolate dispersion to each frequency point.
        W_dat_flip = np.flip(W_dat)      # W_dat and Beta are sorted with lambda, here we sort them with omega
        [Beta_interp, beta0, beta1, beta2] = interp_beta(W_dat_flip,w_interp,B,B1,B2,lambda0) 
        
        # The orign part of gnlse_YY_ODE45
        #  [Z, AT, AW, W] = gnlse_YY_ODE45(t, AT_0, w0, gamma, Beta_interp,beta0,beta1, ...
        #                                 alpha,fr,hr,tau_shock,flength, nsaves, nplot_seg)
        n, dT = len(t), t[1]-t[0] # grid parameters
        V = 2*np.pi*np.arange(-n/2,n/2)/(n*dT) # frequency grid
        
        ############ L operator #############
        
        L = 1j*(Beta_interp - V * beta1 - beta0) - (alpha/2) # Linear operator (dispersion): D=-a/2+i(beta-w*beta_1-beta_0)
        if betas:
            # betas = false, do simu based on real geo 
            # betas = [beta2, beta3,...], do simu based on betas
            B = 0
            betas_V = betas# copy betas to be betas_V so that it can be used with freq V can changed by DV. betas always keeps the same. 
            
            if 'DV' in locals() and 'var' in locals():  # Use `globals()` if they are global variables
                betas_V[0] = betas[0] * DV(us[ii] * L_tot)  # Apply dispersion variation to betas
                print(f"Theory study disp vari at {us[ii] * L_tot * 100:05.1f} cm")
                
            else:
                betas_V[0] = betas[0]
            
            for ibeta in range(len(betas_V)): #ibeta represent index of each element of betas_V
                B = B + betas_V(ibeta) / math.factorial(ibeta+1) * V ** (ibeta+1)
                
                tau_shock = 0 
                alpha = 0 # m-1
                # benchmark with paper 'Parabolic pulse gene by use of a disp-decre
                # fiber with ND' doi: 10.1364/OL.29.000498
    
            L = 1j * (B) - (alpha/2)
        
            L = np.fft.fftshift(L)# fftshift to freq domain
        
        ############ N operator #############
        if abs(w0) > eps:              # if w0>0 then include shock
            #gamma = gamma/w0    
            #W = V + w0                # for shock W is true freq
            shock_W = 1 + V*tau_shock   # FT{1+i*d(tau_shock)/dt}
        else:
            #W = 1                     # set W to 1 when no shock
            shock_W = 1               # set shock term to 1 when no shock
        RW = n * np.fft.ifft(np.fft.fftshift(np.array([RT]).T))   # frequency domain Raman
        shock_W = np.fft.fftshift(shock_W)    # shift to fft space
        ############ parameter set for fractional Raman ############
        #RT = 0 # no raman effect, fr = 0, define function to return the RHS of Eq. (3.13)
    
        N_flag = 0 # count the flag number
        
        # === setup and run the ODE integrator
        
        Z = np.linspace(0, flength, int(nsaves))  # select output z points
        
        print("Z 0",Z[0])
        print("Z end",Z[-1])
        print("Z shape",Z.shape)
        y0 = np.fft.ifft(AT_0).squeeze()
        print("y0", y0.shape)
        print("RW", RW.shape)
        print("shock_W", shock_W.shape)
        
        # === set error control options
        def rhs(z, AW):
           
            AW = AW.squeeze()
            AT = np.fft.fft(AW*np.exp(L*z))         # time domain field
            IT = abs(AT)**2                # time domain intensity
            
            if (RT == int) and (abs(fr) < eps): # no Raman case
                M = np.fft.ifft(AT*IT)             # response function
            else:
                RS = dT*fr*np.fft.fft(np.fft.ifft(IT)*RW) # Raman convolution
                M = np.fft.ifft(AT*((1-fr)*IT + RS))# response function

            R = 1j*gamma*shock_W*M*np.exp(-L*z)
            
            return R
        
        sol = solve_ivp(
        rhs, # RHS function
        t_span = (Z[0], Z[-1]),                      # Integration range
        y0 = np.fft.ifft(AT_0).squeeze(),            # Initial conditions
        method='RK45',                # Similar to MATLAB's ode45
        vectorized=True               # Set to True if rhs is vectorized
        # args=(L, dT, gamma, RW, shock_W, 0, [0])
        )

        Z = sol.t       # The points where the solution was evaluated
        AW = sol.y.T
        
        
        # === process output of integrator
        if nplot <1:
            nplot = 1 # at least 1 frame for each segement
        
        plot_select = np.round(np.linspace(2, int(nsaves), int(nplot))).astype(int) # select only limited simulation results to be ploted to avoid huge graphic demands.
        
        AW = AW[plot_select,:]
        Z = Z[plot_select]
        AT = np.zeros(np.shape(AW[1,:]))
        for j in range(AW.shape[0]):  # Loop over each row
            # Element-wise multiplication and exponentiation (equivalent to MATLAB's `.*`)
            AW[j, :] = AW[j, :] * np.exp(L * Z[j])  # Change variables for AW
            AT[j, :] = np.fft.fft(AW[j, :])         # Perform FFT and store in AT
            AW[j, :] = np.fft.fftshift(AW[j, :]) * dT * n  # Apply fftshift and scale AW
            
        
        W = V + w0 # the absolute frequency grid
        
        """Treatment needed to be done. Simulations of each segment needed to be append to each other 
        in order to have the whole temporal-spatial simulation. Z_z, AT_z,AW_z are results concatenated by diff. WGs with diff. widths"""
        # Inside your loop:
        Z_new = Z + (us[ii] - d_us[ii]) * L_tot  # Calculate new Z values
        
        # Concatenate along rows (axis=0)
        Z_z = np.concatenate((Z_z, Z_new), axis=0)
        AT_z = np.concatenate((AT_z, AT), axis=0)
        AW_z = np.concatenate((AW_z, AW), axis=0)
        
        # Get the last row of AT
        AT_0 = AT[-1, :]  
        
        print('Propa. Dis.:#d \n',us[ii])

    return [Z_z, AT_z, AW_z, W]


    
# def rhs(z, AW, L, dT, gamma, RW, shock_W, fr = 0, RT = 0):
   
#     AW = AW.squeeze()
#     AT = np.fft.fft(AW*np.exp(L*z))         # time domain field
#     IT = abs(AT)**2                # time domain intensity
    
#     if (len(RT) == 1) and (abs(fr) < eps): # no Raman case
#         M = np.fft.ifft(AT*IT)             # response function
#     else:
#         RS = dT*fr*np.fft.fft(np.fft.ifft(IT)*RW) # Raman convolution
#         M = np.fft.ifft(AT*((1-fr)*IT + RS))# response function

#     R = 1j*gamma*shock_W*M*np.exp(-L*z)
    
#     return R
