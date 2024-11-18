"""
Supercontinuum generation code


@author: Gibaek
"""

import numpy as np
import matplotlib.pyplot as plt
from functions import gnlse_YY_ODE45_TaperAdapt

P0 = 600
Ch = 0
L = 0.006 # m


# FREQUENCY AND TIME PARAMETERS                     #

C = 299792458    #Speed of light
h_bar = 6.626e-34 / (2 * np.pi)

N = 2**15                   # number of grid points, enable to plot until 500nm for lbd0 = 1060nm
twidth = 12.5e-12          # width of time window [s]

dt = twidth/N
df = 1/(twidth)
t = np.arange(-N/2, N/2) * dt  # time grid
f = np.arange(-N/2,N/2) *df
w = 2 * np.pi * f
w_s = np.fft.fftshift(w)       # Pulsation vector, shifted
Te = dt
fe = 1 / dt
Nsim = 1

# WAVEGUIDE AND PULSE PARAMETERS (Define them manually here or load them from a database)

# Change parameters here:
lambda0 = 1300*1e-9 
f0 = C/lambda0
w0 = 2*np.pi*f0
Tfwhm = 190e-15 

print('Tfwhm {} fs'.format(Tfwhm*1e15))


isPulseCompression = 0
if isPulseCompression == 1:
    Tfwhm = 10000e-15 
    # Tfwhm = 35e-15 

w_interp = (w+w0)                                                          # Pulsation vector to interpolate beta(w) on the absolute simulation frequency grid

# Temporal FWHM (s)
T0 = Tfwhm/1.76


# Main part
# AT: complex field amplitude in scale of time
# AW: complex field amplitude in scale of frequency not in scale of
# wavelength!!!
N_current = 0
AT_z_acc=[]
AW_z_acc = []
IT_acc = [] # for sum time intensity
AT_acc =[] # for complex field
IW_acc = [] # for sum frequency intensity
A_0_acc = []
A_0_freq = []
Z_z_acc = []
A_T_Nsim = np.zeros([N,Nsim]) #save complex field after each iteration
M = 0 # specifies maximum number of workers


for k_sim in range(Nsim):

    # N_current = N_current+1

    # Add noise: quantum noise and amplitude noise
    # rng('shuffle')
    qsn_f = np.sqrt(h_bar*(w+w0)/df)*np.exp(2 * 1j * np.pi * np.random.rand(N,1).T)               # Quantum Noise, one photon per mode
    qsn_f = qsn_f.squeeze()
    qsn_f[w+w0 <= 0] = 0
    # qsn_f[(w + w0) <= 0] = 0
    rin = 0.97e-2                                                       # Amplitude noise from laser source itself, Relative intensity noise
    # rin = 0.8
    amp_noise = (1+0.5*rin*np.random.rand(1,1))     
    
    # trigger the simulation
    # origin from iSimulation_supercontinuum
    ###########################################################################
    #                                                                         #
    #                       INITIAL PULSE ENVELOPE                            #
    #                                                                         #
    ###########################################################################
    # remove the noise for the moment 
    # amp_noise = 1
    # qsn_f = 0

    isHyperbolic = 0
    isGaussian = 0
    isFilter = 0
    
    if isHyperbolic == 1:
        A_0_acc.append(np.sqrt(P0)*amp_noise*1/np.cosh(t/T0)*np.exp(-1j*Ch*t**2/(2*T0**2)))  # Hyperbolic secant type envelope, add amplitude noise
    elif isGaussian == 1:
        A_0_acc.append(np.sqrt(P0)*amp_noise*np.sqrt(2**(-(2**t/Tfwhm)**2)))           # ref: https://www.rp-photonics.com/gaussian_pulses.html
    else:
        # by defalut use hyperbolic
        A_0_acc.append(np.sqrt(P0) * amp_noise * 1/np.cosh(t/T0) * np.exp(-1j * Ch * t ** 2 / (2 * T0 ** 2)))  # Hyperbolic secant type envelope, add amplitude noise

    A_0_freq.append(np.fft.ifft(A_0_acc[k_sim])*N*Te + np.fft.ifftshift(qsn_f))               # add quantum noise, qsn = quantum shot noise
    # isFilter = 1
    if isFilter == 1:
        UpEdge_f = C/((1550-3.75) * 1e-9) - C/((1550) * 1e-9)
        DownEdge_f = C/((1550+3.75) * 1e-9) - C/((1550) * 1e-9)
        Filter = rectangularPulse(DownEdge_f,UpEdge_f,f)# f not yet ifftshift
        A_0_freq[:,k_sim] = A_0_freq[:,k_sim] * np.fft.ifftshift(Filter)
    
    
    A_0_acc[k_sim] = np.fft.fft(A_0_freq[k_sim])/(N*Te) # return to time domaine, alway taking into account the 'scale' problem
    
    plt.plot(t*1e12,abs(A_0_acc[k_sim].squeeze())**2)
    plt.xlabel('time (ps)')
    plt.ylabel('power (W)')#ylabel('10log(|A_{freq}|^2) (dB)')
    
    # INITIALIZING THE ENVELOPE FOR SIMULATION
    # SIMULATION PARAMETER 
    
    n_h = 20000  # Number of steps
    h = L/n_h  # Step size

    #### ALGORITHM CHOICE ####
    Steep = 1 # Taking into account self-steepening and modal area wavelength dependance (0 or 1)

    #### NUMBER OF SAVED FRAMES ####
    r_lim = 100 # The algorithm saves the results every r_lim steps (important for global simulation time). If r_lim<600, it plots the maps.#}

    ### Report simulation condition ###
    print('Peak power: {} W.'.format(P0))
    print('Pump Wavelength: {} nm.'.format(lambda0*1e9))
    print('{} nm is the short wavelength edge.'.format(round(1e9*1/(1/(2*C*dt)+1/(lambda0))))) # the minimum WL that can be accuratly calculated due to Nyquist sampling thoerem
    # numbers of steps in z to calculate
    nsaves = 2**8


    # numbers of steps to plot evolution through the waveguide, too many steps to plot many cause graphic problem
    #     nplot = 2^8
    # each segment will have to plot nplot_seg.
    #     nplot_seg = round(nplot/length(d_us))
    # prepare for saving simulation data with respect to propagation.

    ## Method Runge-Kutta technique and Interaction Picture
    isImproved = 1
    
    if isImproved == 1:
        RT = 0
        fr = 0

        [Z_z, AT_z, AW_z, W] = gnlse_YY_ODE45_TaperAdapt(t,A_0_acc[k_sim].T, w0, w_interp, W_dat,
                               d_us,us,L,
                               Betas_Dummy,Beta1s_Dummy,Beta2s_Dummy,A_effs_Dummy,n2_effs_Dummy,n_effs_Dummy,
                               alpha, fr, RT,
                               betas,DV)
    
    # report how many simulations have been implemented
    print('{} in {} complete\n'.format(k_sim,Nsim))
    
  
    Z_z_acc.append(Z_z)
    IT_acc.append(abs(AT_z)**2)
    AT_acc.append(AT_z)# for complex field 
    IW_acc.apend(abs(AW_z)**2)
    
    A_T_Nsim.append(AT_z[-1,:]) # the complex field at the end of propagation

# Calculate the average intensity, the captor only mesure average intensity
# because of the integration time.
IT_mean = mean(IT_acc,3)
IW_mean = mean(IW_acc,3)
AT_mean = mean(AT_acc,3)
A_0_mean = mean(A_0_acc,2)
##
### plot ###
# === plot output
T = t# time grid   
dT = T(2)-T(1)
W = w_interp # W and w_interp are the same, except for litter differnet caused by eps
dW = W(2)-W(1)
WL_min = 400*1e-9   # min WL for plot
WL_max = 2000*1e-9  # max WL for plot

### plot the spectrum and pulse evolution through propagation ###
WL, iis = 2*pi*C/W, (WL>WL_min & WL<WL_max) # wavelength grid

IW = IW_mean* 2*pi*C/WL.T**2             # linear scale spectral intensity.
lIW = 10*log10(IW)                        # log scale spectral intensity

# save the start and end spectrum for comparasion.
# log scale
lIW_start = lIW[0,:]
lIW_end = lIW[-1,:]
# linear scale
IW_end = IW[-1,:]

## Save simulation
mlIW = max(max(lIW))       # max value, for scaling plot
# observe with propagation, defined by is_observation. the observation
# point is L/3 2L/3 and L repectively, they can also be changed.
Z = Z_z_acc[:,0]

is_observation = [
    np.argmin(np.abs(Z - L / 6)),
    np.argmin(np.abs(Z - L * 2 / 6)),
    np.argmin(np.abs(Z - L * 3 / 6)),
    np.argmin(np.abs(Z - L * 4 / 6)),
    np.argmin(np.abs(Z - L * 5 / 6)),
    np.argmin(np.abs(Z - L * 6 / 6))
]


#%%
# creat a folder to restore raw data

save_enable = 1

if save_enable == 1:
    try:
        # Create the folder
        os.mkdir(SaveName)
        isSave = True  # Successfully created
    except FileExistsError:
        print("Warning: Folder already exists. Proceeding anyway.")
        isSave = True  # Folder exists, so we proceed
    except Exception as e:
        print(f"Warning!!! Simulation results are not able to be saved! Error: {e}")
        isSave = False  # Failed to create the folder

    # For input pulse: both spectra and complex field is recorded
    IT_0 = abs(A_0_mean)**2
    AW_0 = fftshift(ifft(A_0_mean))*dT*N
    IW_0 = abs(AW_0)**2*2*pi*C/WL.T**2                        # linear scale spectral intensity.
    lIW_0 = 10*np.log10(IW_0)                                     # log scale spectral intensity
    
    PulseIncident.WL = WL(iis)*1e9
    PulseIncident.W = W(iis)
    PulseIncident.I_linear = IW_0(iis)
    PulseIncident.I_dB = lIW_0(iis)
    Table = table(PulseIncident.WL,PulseIncident.W,PulseIncident.I_dB,PulseIncident.I_linear, ...
                'VariableNames',{'wavelength(nm)','omega(Hz*2pi)','transmission(dB)','transmission'})
    writetable(Table,strcat(SaveName,'/spectrum',num2str(0),'.txt'))
    PulseIncident.A_T_real = real(A_0_mean)
    PulseIncident.A_T_imag = imag(A_0_mean)
    PulseIncident.T = T

    Table_pulse = table(PulseIncident.T,PulseIncident.A_T_real,PulseIncident.A_T_imag, ...
                'VariableNames',{'Delay(s)','A_T real','A_T imag'})
    writetable(Table_pulse,strcat(SaveName,'/pulse',num2str(0),'.txt'))

    print('Simu. at #d saved in the path: #s \n',0,SaveName)
    
 # for different propagation distance: only spectra is recorded
    for i = 1:length(is_observation)
            index = is_observation(i)
            # through propagation...
            lIW_us = lIW(index,:)
            AT_us = AT_mean(index,:)# complex field
            # linear scale
            IW_us = IW(index,:)
            # create the struct SCG_mean to save simulation results
            SCG_mean.WL = WL(iis)*1e9 # in nm
            SCG_mean.W = W(iis)
            SCG_mean.I_linear = IW_us(iis)' # linear scale spectral intensity
            SCG_mean.I_dB = lIW_us(iis)'    # log scale spectral intensity
#                 SCG_mean.I_linear_norm = IW_us(iis)'./sum(IW_us*dW) # normalised to 1
#                 SCG_mean.I_dB_norm = (lIW_us(iis)-10*log10(sum(IW_us*dW)))'# normalised to 0dB
        Table = table(SCG_mean.WL,SCG_mean.W,SCG_mean.I_dB,SCG_mean.I_linear, ...
                'VariableNames',{'wavelength(nm)','omega(Hz*2pi)','transmission(dB)','transmission(arbi. unit)'})
        writetable(Table,strcat(SaveName,'/spectrum',num2str(Z(index)),'.txt'))
        Table_FieldComplex_Z = table(T,AT_us',...
                               'VariableNames',{'T(s)','Complex'})
        writetable(Table_FieldComplex_Z,strcat(SaveName,'/ComplexField',num2str(Z(index)),'.txt'))
        print('Simu. at #d saved in the path: #s \n',Z(index),SaveName)

#           plot the spectrum of supercontinuum generation process
         if plot_spectrum_enable == 1
                plot(WL*1e9,lIW_us','LineWidth',2,'DisplayName',num2str(Z(index)))
                xlabel('\lambda (nm)')
                ylabel('P (dB) 10log(|A_{freq}|^2)')#ylabel('10log(|A_{freq}|^2) (dB)')
#                     legend('Input spectrum','Output spectrum','Location','northeast')
                axis([WL_min*1e9 WL_max*1e9 mlIW-50 mlIW])   
                title(SaveName,'Interpreter','none')
                grid minor
                p_Spectra = gcf
                hold on

if plot_spectrum_enable == 1
    plot(PulseIncident.WL,PulseIncident.I_dB,'LineWidth',2,'DisplayName','incident pulse')

legend()
hold off
# Save different complex field at the end of propagation, there will be
# Nsim simulations
Table_ComplexField = table(T,A_T_Nsim)
writetable(Table_ComplexField,strcat(SaveName,'/complex_field','.txt'))

# Save Geo. d_us, us and widths
Table_geo = table(us*L,d_us*L,widths, ...
        'VariableNames',{'Propa. dis.(m)','Seg. dis.(m)','width(nm)'})
writetable(Table_geo,strcat(SaveName,'/Geo','.txt'))


 # plot spectrum evolution
mlIW = max(max(lIW))       # max value, for scaling plot
 
if plot_evo_enable == 1
    f_evo = figure()
    f_evo.Position(3:4)=[280*4 210*2]
    subplot(1,2,1)             
    pcolor(WL(iis)*1e9, Z, lIW(:,iis)) # plot as pseudocolor map
    colormap turbo
    caxis([mlIW-60.0, mlIW])  xlim([WL_min*1e9,WL_max*1e9]) shading interp 
    xlabel('Wavelength / nm') ylabel('Distance / m')
    colorbar


    # pulse evolution
    lIT = 10*log10(IT_mean) # log scale temporal intensity
    mlIT = max(max(lIT))       # max value, for scaling plot

    subplot(1,2,2)
    pcolor(T*1e12, Z, lIT)    # plot as pseudocolor map
    caxis([mlIT-40.0, mlIT])  xlim([-5,5]) shading interp
    xlabel('Delay / ps') ylabel('Distance / m')
    colorbar
    p_SpectraEvo = gcf
    hold off
    
    
# strcat(SaveName,'/','Spectra','.png')
if save_enable == 1
exportgraphics(p_Spectra,strcat(SaveName,'/','Spectra','.png'),'Resolution',600)
exportgraphics(p_SpectraEvo,strcat(SaveName,'/','SpectraEvo','.png'),"Resolution",600)
exportgraphics(p_Geo,strcat(SaveName,'/','Geo','.png'),'Resolution',600)

if Nsim >=21 # only save raw data in SERIOUS simulation
    isSaveEvoData = 1
else
    isSaveEvoData = 0
end

if isSaveEvoData == 1
    # save the raw evolution data
     filename_Evo =strcat(SaveName,'/','Evo','.mat')
     WL_Evo = WL(iis)*1e9
     lIW_Evo = lIW(:,iis)
     T_Evo = T*1e12
     lIT_Evo = lIT
     save(filename_Evo,'WL_Evo','Z','lIW_Evo','T_Evo','lIT_Evo')
     
#          iis_T = (T_Evo>-5 & T_Evo<5) # time grid, only save time within +-5ps
     iis_T = (T_Evo>-10 & T_Evo<10) # time grid, only save time within +-5ps

     T_Evo = T_Evo(iis_T)
     lIT_Evo = lIT_Evo(:,iis_T)
     n_T = round(linspace(1,length(T_Evo),2048))# only choose 2048 different z to save space.
     # special matrix prepared for origin
     if length(Z)< 512
         Z_steps = length(Z)
     else
         Z_steps = 512
     end
     n_z = round(linspace(1,length(Z),Z_steps))# only choose 512 max different z to save space.
     n_wl = round(linspace(1,length(WL_Evo),2048))# only choose 2048 different z to save space.
#          nImages = length(n)
     
     filename_Evo_XYZ =strcat(SaveName,'/','Evo_W_XYZ','.txt')
     [X,Y] = meshgrid(WL_Evo(n_wl),Z(n_z))
     lIW_Evo = lIW_Evo(n_z,n_wl)
     lIW_XYZ = [X(:) Y(:) lIW_Evo(:)]
     writematrix(lIW_XYZ,filename_Evo_XYZ)
     
     filename_Evo_XYZ =strcat(SaveName,'/','Evo_T_XYZ','.txt')
     [X,Y] = meshgrid(T_Evo(n_T),Z(n_z))
     lIT_Evo = lIT_Evo(n_z,n_T)
     lIT_XYZ = [X(:) Y(:) lIT_Evo(:)]
     writematrix(lIT_XYZ,filename_Evo_XYZ)

end
# 

#
end
## create spectrogram 
isCreateSpectrogram = 0
if isCreateSpectrogram == 1
signal_sg = A_T_Nsim(:,1)
spectrogram_Yijun(signal_sg)
end

## create gif
isCreatGif = 1
if isCreatGif == 1
# spectrum
n = round(linspace(1,length(Z),256))
nImages = length(n)
fig = figure('Renderer', 'painters', 'Position', [10 10 1350 600])
ParaFunInterp = griddedInterpolant(cat(1,0,us),cat(1,widths(1),widths),'linear')# interp geo vs unitied propa dis. in order to show in gif
for idx = 1:nImages
    plot(WL(iis)*1e9,lIW(n(idx),iis),"LineWidth",2,'DisplayName',strcat('width: ',num2str(round(ParaFunInterp(Z(n(idx))/L))),' nm'))
    legend()
    title("Propa. dis. = " + num2str(Z(n(idx))) + 'm')
    xlabel('\lambda (nm)')
    ylabel('P (dB) 10log(|A_{freq}|^2)')#ylabel('10log(|A_{freq}|^2) (dB)')
    axis([WL_min*1e9 WL_max*1e9 mlIW-50 mlIW])
    drawnow
    frame = getframe(fig)
    im{idx} = frame2im(frame)

filename = strcat(SaveName,'/','lIW','.gif') # Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256)
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.05)
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.05)


# temporal pulse
nImages = length(n)
fig = figure('Renderer', 'painters', 'Position', [10 10 1350 600])
for idx = 1:nImages
    plot(T*1e12,lIT(n(idx),:),"LineWidth",2,'DisplayName',strcat('width: ',num2str(round(ParaFunInterp(Z(n(idx))/L))),' nm'))
    legend()
#     caxis([mlIT-40.0, mlIT])  xlim([-5,5])
#         axis([-5,5 mlIT-50 mlIT])
    axis([-10,10 mlIT-50 mlIT])

    title("Propa. dis. = " + num2str(Z(n(idx))) + 'm')
    xlabel('Delay / ps')
    ylabel('I (dB)')#ylabel('10log(|A_{freq}|^2) (dB)')
    drawnow
    frame = getframe(fig)
    im{idx} = frame2im(frame)

end

close
filename = strcat(SaveName,'/','lIT','.gif') # Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256)
    if idx == 1:
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.05)
    else:
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.05)

### plot the input and output spectrum comparison chart in frequency domain
### 

### plot the input and ouput spectrum comparison chart in wavelength
### domaine ###
if plot_spectrum_enable == 1
    figure()
    plot(WL*1e9,lIW_start,WL*1e9,lIW_end,'LineWidth',2)
    #title('Wavelength spectrum (input and output)')
    xlabel('\lambda (nm)')
    ylabel('P (dB) 10log(|A_{freq}|^2)')#ylabel('10log(|A_{freq}|^2) (dB)')
    legend('Input spectrum','Output spectrum','Location','northeast')
    axis([WL_min*1e9 WL_max*1e9 mlIW-40 mlIW])   
    title(SaveName,'Interpreter','none')
    grid minor
    p_Spectra = gcf

if save_enable == 1
        exportgraphics(p_Spectra,Spectra,'Resolution',600)
        exportgraphics(p_SpectraEvo,SpectraEvo,"Resolution",600)
