# -*- coding: utf-8 -*-
"""
Import files containing dispersion data from mode solver

@author: Gibaek
"""
import os
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

dir_path = os.path.dirname(os.path.abspath(__file__))

def Dummy_simu(Width_Dummy, W_dat,Lambda,Fs_Beta,Fs_Beta1,Fs_Beta2,Fs_n_eff,Fs_A_eff,Fs_n2_eff, isPlot):
    simu_num = len(W_dat)
    Beta_Dummy = np.zeros(simu_num)
    Beta1_Dummy = np.zeros(simu_num)
    Beta2_Dummy = np.zeros(simu_num)
    n_eff_Dummy = np.zeros(simu_num)
    A_eff_Dummy = np.zeros(simu_num)
    n2_eff_Dummy =np. zeros(simu_num)

    for i in range(len(Fs_Beta2)):
        Beta2_Dummy[i] = Fs_Beta2[i](Width_Dummy)
        Beta_Dummy[i] = Fs_Beta[i](Width_Dummy)
        Beta1_Dummy[i] = Fs_Beta1[i](Width_Dummy)
        n_eff_Dummy[i] = Fs_n_eff[i](Width_Dummy)
        A_eff_Dummy[i] = Fs_A_eff[i](Width_Dummy)
        n2_eff_Dummy[i] = Fs_n2_eff[i](Width_Dummy)

    if isPlot == 1:
        lgd = ('width {} nm'.format(Width_Dummy))
        plt.plot(Lambda*1e9,Beta2_Dummy*1e24,'o', label = lgd);
        plt.ylabel(['\beta_{2} ','ps^{2} m^{-1}'])
        plt.xlabel('Wavelength (nm)')
        plt.title('TE polarisation')
        plt.legend()
        plt.plot(Lambda*1e9, 0*Lambda*1e9,'--')
    
    return [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy]


def DispersionVariation(z):

    # g = 46 # m-1
    # g = 125
    g=0
    
    D_z = 1/(1 + g*z)
    
    return D_z


isPlot = 1
isPlotGVD = 1

# alpha = math.log(10)/10*400
alpha = 0
C = 299792458    # Speed of light [m/s]

isAND_DW = 0

if isAND_DW == 1:
    print('WARNING!! all GVD curve are modified for all normal dispersion DW generation!\n')

#%% Load all data from simulations.

isEff = 0

if isEff == 1:
    print('WARNING!! effective geo is on, width is not real value!\n')
else:
    FolderName = 'SCG simulation JH'
    SampleName = np.array([600,1200,2000]) # Sample variable


isTM = 0

file_name_list = []

for ii in range(len(SampleName)):
    if isTM == 0:
        file_name_list.append('SiN_Strip_{}nm_h_700nm_TE_240731.txt'.format(SampleName[ii]))
      
    else:
        file_name_list.append('SiN_Strip_{}nm_h_700nm_TM_240731.txt'.format(SampleName[ii]))
        
        print('WARNING!! for TM mode!\n')
    
    if isEff ==1:
        print(file_name_list[ii])

Lambda_list = []
W_dat_list = []
Beta_list = []
Beta1_list = []
Beta2_list = []
n_eff_list = []
A_eff_list = []
n2_eff_list = []

for ii in range(len(file_name_list)):
   
    file_name = 'width 640nm D 3.xlsx'
    print('loading.. ',file_name,'\n')
    
    Mb = pd.read_excel(dir_path + '\\' + file_name)
    Mb = Mb.to_numpy()
    
    Lambda = Mb[:, 0] * 1e-9
    W_dat = Mb[:, 1]
    Beta = Mb[:, 2]
    Beta1 = Mb[:, 3]
    Beta2 = Mb[:, 4]
    n_eff = Mb[:, 5]
    A_eff = Mb[:, 6]
    n2_eff = Mb[:, 7]
    
    Lambda_list.append(Mb[:,0] * 1e-9)    
    W_dat_list.append(Mb[:, 1])
    Beta_list.append(Mb[:, 2])
    Beta1_list.append(Mb[:, 3])
    Beta2_list.append(Mb[:, 4])
    n_eff_list.append(Mb[:, 5])
    A_eff_list.append(Mb[:, 6])
    n2_eff_list.append(Mb[:, 7])

# convert list to numpy array
Lambda_list, W_dat_list = np.array(Lambda_list), np.array(W_dat_list)
Beta_list, Beta1_list, Beta2_list = np.array(Beta_list), np.array(Beta1_list), np.array(Beta2_list)
n_eff_list, A_eff_list, n2_eff_list  = np.array(n_eff_list), np.array(A_eff_list), np.array(n2_eff_list)

print('All data from simulation has been loaded.\n')
#%% Creation all regression functions to generate dummy simulation.

Widths_simu = SampleName # numer of variation in width 

Fs_Beta = []
Fs_Beta2 = []
Fs_Beta1 = []

Fs_n_eff = []
Fs_A_eff = []
Fs_n2_eff = []

for ii in range(len(Mb)):
    # Ensure the number of sample data for interpolation
    if len(Widths_simu) >= 2 and len(Beta2_list[:,0]) >= 2:
        
        F = make_interp_spline(Widths_simu, Beta_list[:,ii], k = 2)
        Fs_Beta.append(F)
        
        F = make_interp_spline(Widths_simu, Beta2_list[:,ii], k = 2)
        Fs_Beta2.append(F)
        
        F = make_interp_spline(Widths_simu, Beta1_list[:,ii], k = 2)
        Fs_Beta1.append(F)
        
        F = make_interp_spline(Widths_simu, n_eff_list[:,ii], k = 2)
        Fs_n_eff.append(F)
        
        F = make_interp_spline(Widths_simu, A_eff_list[:,ii], k = 2)
        Fs_A_eff.append(F)
        
        F = make_interp_spline(Widths_simu, n2_eff_list[:,ii], k = 2)
        Fs_n2_eff.append(F)
        
    else:
        print('Not enough data points for interpolation.')
    
    
    if isPlot == 1:
        plt.plot(Widths_simu, Beta2_list[:,ii] * 1e24, 'o')
        plt.ylabel(['\beta_{2} ', 'ps^{2} m^{-1}'])
        plt.xlabel('Width (nm)')
        plt.title('TE polarisation')
        plt.legend()
        
        Widths_interp = np.linspace(700, 1500, 81)
        plt.plot(Widths_interp, Fs_Beta2[ii](Widths_interp) * 1e24, '-')
        
print('All regression functions to generate dummy simulation data created.\n')

#%% Define waveguide geometry, Dispersion managed WG geo creation
# with unified parametry u from 0 to 1

L = 0.006 # m

# ParaFun = lambda u: 3200+u*(-2000)
# ParaFun = lambda u: (np.where(u <= 0.3, 3200, 3200 + u * -2000))
ParaFun = lambda u: 640 + u * 0

isTailor = 0
isOHY8P13 = 0
isDM_design = 0
betas = -9.661465159831380e-25
betas= False

"""
Define DV as a function handle.
It should only work when betas has values which means we are in theoritical study mode.
later it will be used as betas * DV as a function of propagation distance us*L
"""
if betas:
    DV = DispersionVariation
    
else:
    DV = lambda u: 1

Width_Tol = 25                                                               # tolereance in nm
num_discretization = 2**15
Lp = np.linspace(0,1,num_discretization)# propa. dis. discretization.
ws_discretization = ParaFun(Lp)

### For periodical structure ###
num_period = 1
Lp = np.linspace(0,1,num_discretization*num_period)
ws_discretization = np.tile(ws_discretization,(1,num_period))

# New version of discretization of WG geo
d_Lp = Lp[1]-Lp[0]

j = 0
d_us = [0]
us = [0]

while us[j] < 1.0:
    D_Lp = d_Lp #initialize by basic step defined by d_Lp
    while abs(ParaFun(us[j]) - ParaFun(us[j] + D_Lp)) <= Width_Tol:
        if us[j] + D_Lp > 1:
            break
        
        D_Lp = D_Lp+d_Lp    
    D_Lp = D_Lp-d_Lp # go back to tolerance
#     print('Width difference right now is:%d nm.\n',abs(ParaFun(us(j))-ParaFun(us(j)+D_Lp)))
    if D_Lp <= 0:
        print('Warning!! Cannot fulfill tolerance demand, please change to more indense Lp.\n')
        print('Problem happens at us(j) = {}.\n'.format(us[j]))
        break
    else:
        us.append(us[j]+D_Lp)
        d_us.append(D_Lp)
       
    j = j+1
    # make sure propagation distance do not exceed 1.
    if us[j] >= 1:
        us[j] = 1
        d_us[j] = 1 - us[j-1]
    
    print('Propa Dis u: {}.\n'.format(us[j]))


del us[0]# to remove the first element 0
del d_us[0] # to remove the first element 0

us = np.array(us).conj().T
d_us = np.array(d_us).conj().T
widths = ParaFun(us)

### consider periodical stucture###
d_us = np.tile(d_us,(num_period,1))
d_us = d_us/num_period #renormalised to 1
widths = np.tile(widths,(1,num_period))

if num_period>1:
    us_p = us
    for k in range(1,num_period): # 1이 아니라 2일 수도
        us = np.hstack(1,us,(us_p+k-1)) # vstack 일수도

us = us/num_period # renormalised to 1

plt.title('Waveguide length: {} cm'.format(L*100))#,', tolerance:',num2str(Width_Tol),'nm'))
plt.plot(us*L,widths,'-o')
plt.xlim([0, L*2])
print('{} segements at tolerence of {} nm.\n'.format(len(d_us),Width_Tol))
print('Dispersion managed WG geo created.\n')
plt.legend()


#%% find dummy simulation data at different given width.
Widths_Dummy = widths
Width_num = len(Widths_Dummy)
W_num = len(Mb)
Betas_Dummy = np.zeros([W_num,Width_num])
Beta1s_Dummy = np.zeros([W_num,Width_num])
Beta2s_Dummy = np.zeros([W_num,Width_num])
n_effs_Dummy = np.zeros([W_num,Width_num])
A_effs_Dummy = np.zeros([W_num,Width_num])
n2_effs_Dummy = np.zeros([W_num,Width_num])

for i in range(len(Widths_Dummy)):
    Width_Dummy = Widths_Dummy[i]
    [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy] = Dummy_simu(Width_Dummy, W_dat,Lambda
                                                                                           ,Fs_Beta, Fs_Beta1,Fs_Beta2,
                                                                                           Fs_n_eff,Fs_A_eff,Fs_n2_eff, isPlot)
    Betas_Dummy[:,i] = Beta_Dummy
    Beta1s_Dummy[:,i] = Beta1_Dummy
    Beta2s_Dummy[:,i] = Beta2_Dummy
    n_effs_Dummy[:,i] = n_eff_Dummy
    A_effs_Dummy[:,i] = A_eff_Dummy
    n2_effs_Dummy[:,i] = n2_eff_Dummy

print('All dummy simu data created.\n')
