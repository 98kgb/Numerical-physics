# -*- coding: utf-8 -*-
"""
Test for the conversion (Matlab to python)

@author: Gibaek
"""

import numpy as np

C = 299792458    #Speed of light

save_enable = 1 # Saving options, 1 stands for save, 0 doesn't save. always put it 1
isMat = 0
plot_evo_enable = 1
plot_spectrum_enable = 1
plot_beta2_enable = 1

##
isTM = 0
Nsim = 3  # set simulation iteration times
lambda0 = 1300*1e-9   #m
Ch = 0# Initial chirp
betas = [-1.172554609910727e-25, -4.291112339397303e-39]
betas = False

TaperDispersion # TaperDispersion script check, for dispersive

P0_list=[600]

for k in range(len(P0_list)):
    P0 = P0_list(k)
    FolderName = 'SCG_simu_results\SiN_240801'

    isPulseCompression = 0
    isHyperbolic = 1
    isFilter = 0
    isGaussian = 0

    if isTM == 1:
        FolderName = FolderName+'\\'+'TM'
    else :
        FolderName = FolderName, '\\', 'TE'
    
    
    DeviceName = '190fs_SiN_W640_Degree 3 L {}',foramt(L)


    if isMat == 1:
        SaveName = FolderName + '//' + DeviceName + '_AVE{}'.format(Nsim) + '_lbd0nm_{}'.format(lambda0*1e9)+'_P0W_{}'.format(P0)+'.mat'
    else:
        SaveName = FolderName+'//'+DeviceName+'_AVE{}'.format(Nsim)+'_lbd0nm_{}'.format(lambda0*1e9)+'_P0W_{}'.format(P0)
    

    Spectra = FolderName+'/'+DeviceName+'_AVE{}'.format(Nsim)+'_lbd0nm_{}'.format(lambda0*1e9)+'_P0W_{}'.format(P0)+'.png'
    SpectraEvo = FolderName+'/'+DeviceName+'_AVE{}'.format(Nsim)+'_lbd0nm_{}'.format(lambda0*1e9)+'_P0W_{}'.format(P0)+'_Evo.png'
    
    isSave = mkdir(SaveName)    # creat a folder to restore raw data
    if isSave == 0:
        fprintf('Warning!!! Simulation results are not able to be saved!')
        
    
        
#         
    # the effective parameter of a str WG equivlent to real geo
    Omega_eff = 2*pi*C/Lambda
    Beta_eff = Betas_Dummy*d_us
    Beta1_eff = Beta1s_Dummy*d_us
    Beta2_eff=Beta2s_Dummy*d_us#effective beta2, the average GVD of all waveguide
    n_eff_eff = n_effs_Dummy*d_us
    A_eff_eff = A_effs_Dummy*d_us
    ng_eff = Beta_eff*0
    n2_eff_eff = Beta_eff*0
    D_eff = Beta_eff*0

    # Table_geo_eff = table(Lambda*1e9,Omega_eff,Beta_eff,Beta1_eff,Beta2_eff,n_eff_eff,A_eff_eff,ng_eff,n2_eff_eff,D_eff, ...
    #                     'VariableNames',{'Wavelength(nm)','Omega(rad/s)','Beta0(1/m)',...
    #                     'Beta1(s/m)','Beta2(s^2/m)','neff','Aeff(m^2)','ng','n2_eff','D(ps/(nm*km)'})
    # writetable(Table_geo_eff,strcat(SaveName,'/geo_eff','.txt'),'Delimiter','\t')    
    
    SCG_average

n2 = 2.4e-19
Tfwhm = 190e-15
Dlength = Tfwhm**2/abs(Beta2)
gamma = w0*n2/(C*A_eff_eff)
NLlength = 1/(P0_list*gamma)
NofSoliton = (Dlength/NLlength)**(1/2)
Ratiooflengths = Dlength/NLlength
properties = [gamma, NLlength, Dlength, NofSoliton, Ratiooflengths]


