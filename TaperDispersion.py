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

isPlot = 0
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
        print(strcat(file_name_list(i),'\n'))

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
        figure()
        lgd = strcat('@', num2str(Lambda_list[ii] * 1e9), 'nm')
        plt.plot(Widths_simu, Beta2_list[:,ii] * 1e24, 'o', 'DisplayName', label = lgd)
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
betas= false

"""
Define DV as a function handle.
It should only work when betas has values which means we are in theoritical study mode.
later it will be used as betas * DV as a function of propagation distance us*L
"""
if betas != False:
    # DV = @DispersionVariation
    print("need to be checked....")
else:
    DV = lambda u: u = 1

"""여기부터 고쳐야함!"""

Width_Tol = 25                                                               # tolereance in nm
plt.figure(0)
num_discretization = 2**15
Lp = np.linspace(0,1,num_discretization)# propa. dis. discretization.
ws_discretization = ParaFun(Lp)

### For periodical structure ###
num_period = 1
Lp = np.linspace(0,1,num_discretization*num_period)
ws_discretization = repmat(ws_discretization,1,num_period)
###


# plt.plot(Lp,ws_discretization,'-','LineWidth',0.5,'DisplayName','Geo. origin')
plt.xlabel('Distance m')
plt.ylabel('Width nm')
# plt.title(strcat('Waveguide length:  ',num2str(L*100),'cm',', tolerance:',num2str(Width_Tol),'nm'))

# h = 0.02       # step size
# X = 0:h:1     # domain
# f = ParaFun(X)# range
# Y = diff(f)/h   # first derivative
# # Z = diff(Y)/h   # second derivative
# plot(X(:,1:length(Y)),Y,'r',X,f,'b','LineWidth',1.5)#, X(:,1:length(Z)),Z,'k')
# ### linear regression for first derivative ###
# F_dY = griddedInterpolant(X(:,1:length(Y)),Y,'linear')


### calculate steps from u=0 ###

# New version of discretization of WG geo
d_Lp = Lp[1]-Lp[0]
j = 1
d_us=[]
d_us(1) = 0
us = []
us(1) = 0

while us(j)<1.0
    D_Lp = d_Lp #initialize by basic step defined by d_Lp
    while abs(ParaFun(us(j))-ParaFun(us(j)+D_Lp))<=Width_Tol
        if us(j)+D_Lp >1
            break
        end

        D_Lp = D_Lp+d_Lp
    end
    D_Lp = D_Lp-d_Lp # go back to tolerance
#     print('Width difference right now is:%d nm.\n',abs(ParaFun(us(j))-ParaFun(us(j)+D_Lp)))
    if D_Lp<=0
        print('Warning!! Cannot fulfill tolerance demand, please change to more indense Lp.\n')
        print('Problem happens at us(j) = %f.\n',us(j))
        break
    else
        us(j+1) = us(j)+D_Lp
        d_us(j+1) = D_Lp
    end
    j = j+1
    # make sure propagation distance do not exceed 1.
    if us(j)>=1
        us(j) = 1
        d_us(j) = 1-us(j-1)
    end
    print('Propa Dis u:%d.\n',us(j))
end





us(1) =[] # to remove the first element 0
d_us(1) =[] # to remove the first element 0
us = us'
d_us = d_us'
widths = ParaFun(us)
### consider periodical stucture###
    d_us = repmat(d_us,num_period,1)
    d_us = d_us./num_period #renormalised to 1
    widths = repmat(widths,1,num_period)
    if num_period>1
        us_p = us
        for k = 2:num_period
            us = cat(1,us,(us_p+k-1))
        end
    end
    us = us./num_period # renormalised to 1
###
# for w= 1400 str wg 0.008m or 0.01m around 700W
# for w= 1100 str wg 0.006m for power around 700W 
# for w= 1100 str wg 0.015m for power 200W @1310nm
# 
if isTailor == 1
#     L = 1152%m
%     us = linspace(0,L,257)'
    L = 0.02#m
    us = linspace(0,L,129)'
    d_us = diff(us)
    us(1) = []
    widths = 999.*ones(size(us))% no real waveguide widths
    us = us/L
    d_us = d_us/L
end

# 
if isTailor == 1
# #     simuation 2 segment: 1 str and 1 tapered WG for soliton shifting(dispersion decreasing design).
    L1 = 40*5*1e-6#m
    L2 = 10*5*1e-6
    W1 = 1000#nm
    W2 = 800
    
    
    L1 = 50*1e-6#m
    L2 = 50*1e-6
    W1 = 1700#nm
    W2 = 2300
    
    L1 = 50*1e-6#m
    L2 = 50*1e-6
    W1 = 1700#nm
    W2 = 1000
    W1 = 1400#nm
    W2 = 1000
 
    
    W1 = 1000#nm
    W2 = 1400
    d_us = []
    widths = []
    for i = 1:200
        d_us = cat(1,d_us,L1,L2)
        widths = cat(1,widths,W1,W2)
    end
#     d_us = d_us'
#     widths = widths'
    us = cumsum(d_us)
    L = us(end)
    us = us/L
    d_us = d_us/L

end

if isTailor == 1


    Period = 50 * 1e-6# period of QPM waveguide
    Ampli = 250# nm
    Bias = 1450# nm
    
    Period = 50*5 * 1e-6# period of QPM waveguide
    Ampli = 200# nm
    Bias = 1200# nm
    
    d_us = Period/8 * ones(160*4,1)# each sampling is 1/8 of the period
    us = cumsum(d_us)
    widths = Ampli* sin(2*pi/Period * us) + Bias
    L = us(end)
    us = us/L
    d_us = d_us/L


end


# if isTailor == 1
# #     L = 0.06
#     us = cat(1,linspace(0.01,0.05,64)')
#     d_us = cat(1,0.008,diff(us))
#     widths_interm = us
#     widths_interm(1) = []
# #     widths = cat(1,1400,widths_interm.*(-700)./(L-0.008)+1550)
#     widths = cat(1,linspace(1400,900,length(us))')
#     us = cat(1,us,L)
#     d_us = cat(1,d_us,L-0.05)
#     widths = cat(1,widths,2500)
#     us = us/L
#     d_us = d_us/L
# end
if isOHY8P13 == 1
    # try to simulation the quais periodical width modulation of OHY8P13
    [propa_dis,widths,L] = ParaFunCreator_us()
#     
#     if L>0.015*1e6
#         L = 0.015*1e6
#         widths(propa_dis>0.015*1e6)=[]
#         propa_dis(propa_dis>0.015*1e6)=[]
#     end
#     
    us = propa_dis./L
    d_us = diff(us)
    widths(1) = []
    us(1) = []
    L = L*1e-6# convert to m
end
if isDM_design == 1
    [propa_dis,widths,L] = DM_design()
    us = propa_dis
    d_us = diff(us)
    widths(1) = []
    us(1) = []

end
title(strcat('Waveguide length:  ',num2str(L*100),'cm'))#,', tolerance:',num2str(Width_Tol),'nm'))
plot(us.*L,widths,'-o','LineWidth',1,'DisplayName','Geo. simu.')
xlim([0 L])
print('%d segements at tolerence of %d nm.\n',length(d_us),Width_Tol)
print('Dispersion managed WG geo created.\n')
legend()
p_Geo = gcf
hold off

# find dummy simulation data at different given width.
Widths_Dummy = widths
Width_num = length(Widths_Dummy)
W_num = length(W_dat)
Betas_Dummy = zeros(W_num,Width_num)
Beta1s_Dummy = zeros(W_num,Width_num)
Beta2s_Dummy = zeros(W_num,Width_num)
n_effs_Dummy = zeros(W_num,Width_num)
A_effs_Dummy = zeros(W_num,Width_num)
n2_effs_Dummy = zeros(W_num,Width_num)

for i = 1:length(Widths_Dummy)
    Width_Dummy = Widths_Dummy(i)
    [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy] = Dummy_simu(Width_Dummy, ...
        W_dat,Lambda,Fs_Beta,Fs_Beta1,Fs_Beta2,Fs_n_eff,Fs_A_eff,Fs_n2_eff, ...
        isPlot)
    Betas_Dummy(:,i) = Beta_Dummy
    Beta1s_Dummy(:,i) = Beta1_Dummy
    Beta2s_Dummy(:,i) = Beta2_Dummy
    n_effs_Dummy(:,i) = n_eff_Dummy
    A_effs_Dummy(:,i) = A_eff_Dummy
    n2_effs_Dummy(:,i) = n2_eff_Dummy

end

hold off
print('All dummy simu data created.\n')


if isPlot == 1
    figure(5)
    lgd = strcat('width ',num2str(Width_Dummy),'nm')
    plot(Lambda*1e9,Beta2_Dummy*1e24,'o','DisplayName',lgd)
    ylabel(['\beta_{2} ','ps^{2} m^{-1}'])
    xlabel('Wavelength (nm)')
    title('TE polarisation')
    legend()
    hold on
    plot(Lambda*1e9,0.*Lambda*1e9,'Color','black','LineWidth',1.5,'LineStyle','--')
    % compare to real simulaion results FDTD
    Mb = importdata('Dispersion_file_Lumerical_ModeSolver/OHY8P13_1000nm_height_800nm_TE_221011.txt')        # Waveguide properties from simulation
    Mb_dat = Mb.data
    Lambda = Mb_dat(:,1)*1e-9
    u = find(Lambda>=500e-9)
    Lambda = Mb_dat(u,1)*1e-9
    W_dat = Mb_dat(u,2)
    Beta = Mb_dat(u,3)
    Beta1 = Mb_dat(u,4)
    Beta2 = Mb_dat(u,5)
    plot(Lambda*1e9,Beta2*1e24,'r','DisplayName','simulation 1000nm')
    hold on
end

function [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy] = Dummy_simu(Width_Dummy, ...
    W_dat,Lambda,Fs_Beta,Fs_Beta1,Fs_Beta2,Fs_n_eff,Fs_A_eff,Fs_n2_eff, ...
    isPlot)
    simu_num = length(W_dat)
    Beta_Dummy = zeros(simu_num,1)
    Beta1_Dummy = zeros(simu_num,1)
    Beta2_Dummy = zeros(simu_num,1)
    n_eff_Dummy = zeros(simu_num,1)
    A_eff_Dummy = zeros(simu_num,1)
    n2_eff_Dummy = zeros(simu_num,1)
    for i = 1:length(Fs_Beta2)
        Beta2_Dummy(i) = Fs_Beta2{i}(Width_Dummy)
        Beta_Dummy(i) = Fs_Beta{i}(Width_Dummy)
        Beta1_Dummy(i) = Fs_Beta1{i}(Width_Dummy)
        n_eff_Dummy(i) = Fs_n_eff{i}(Width_Dummy)
        A_eff_Dummy(i) = Fs_A_eff{i}(Width_Dummy)
        n2_eff_Dummy(i) = Fs_n2_eff{i}(Width_Dummy)

    end
    if isPlot == 1
        figure(5)
        lgd = strcat('width ',num2str(Width_Dummy),'nm')
        plot(Lambda*1e9,Beta2_Dummy*1e24,'o','DisplayName',lgd)
        ylabel(['\beta_{2} ','ps^{2} m^{-1}'])
        xlabel('Wavelength (nm)')
        title('TE polarisation')
        legend()
        hold on
        plot(Lambda*1e9,0.*Lambda*1e9,'Color','black','LineWidth',1.5,'LineStyle','--')
    end
end

# functions

function D_z = DispersionVariation(z)

    g = 46 %m-1
    g = 125
    g=0
    D_z = 1./(1 + g.*z)
end
