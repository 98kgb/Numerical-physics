function [Z_z, AT_z, AW_z, W] = gnlse_YY_ODE45_TaperAdapt(t, A_0, w0,w_interp,W_dat, ...
                            d_us,us,L_tot,...
                            Betas_Dummy,Beta1s_Dummy,Beta2s_Dummy,A_effs_Dummy,n2_effs_Dummy,n_effs_Dummy, ...
                            alpha, fr, RT, ...
                            betas,DV)
% This function is prepared for SCG, rewrite to ODE45 interaction
% picture method.
% T: time vector
% A: input field in time domain
% w0: pump center freq in radius
% gamma: nonlinear coef
% Beta_interp, beta0, beta1: dispersion parameters to include all higher
% order dispersion contributions.
% alpha: /m, attenuation coefficient
% fr, RT: for Raman response
% betas: do simulations based on betas, real geo are ignored
% DV: dispersion variation function, it should work together with betas
% tau_shock: gamma1/gamma0, to include self-steepening effect
% flength: the length of waveguide
% nsaves: numbers of step when doing simulation
% nplot: numbers of frames to plot to save graphic demands.
%%
Z_z = [];
AT_z = [];
AW_z = [];
isHybridMode = 0;
%% Method Runge-Kutta technique and Interaction Picture
tic;
AT_0 = A_0;% the initial AT_0 for NL Eq. is pulse A_0
c = physconst('LightSpeed');
C=c;
lambda0= 2*pi*c/w0 ;
for i = 1:length(d_us)
    % This method takes into account beta(w) without approximation (no Taylor expansion), it needs to define beta(w)
    % exact length of each segement, depend on us(i)
    flength = L_tot*d_us(i);  
    nplot = round(d_us(i)/0.002);% try to record at each 0.002(normalised propa. dis.)
    
    % slip step method: need to set the steps for dispersion
%     nsaves = round(flength/(20*1e-6));% by default 20µm, one dispersion operator.
      nsaves = round(flength/(10*1e-6));% by default 10µm, one dispersion operator.

%     nsaves = round(flength/(200*1e-6));% by default 0.2mm, one dispersion operator.
%     nsaves = round(flength/(8*1e-2));% by default 8cm, one dispersion operator. for benchmark with DOI10.1364/OL.29.000498


    if nsaves<=3 % for ODE, a minimum number of 3 is needed?
        nsaves = 3;
    end
    % Read pseudo simulation data, so called 'Dummy'.
    B = flip(Betas_Dummy(:,i));
    B1 = flip(Beta1s_Dummy(:,i));
    B2 = flip(Beta2s_Dummy(:,i));
    A_eff = A_effs_Dummy(:,i);
    n2_eff = n2_effs_Dummy(:,i);
    n_eff = n_effs_Dummy(:,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for gamma and tau_shock 
    n_interp=2^12;
    n2 = 2.4*1e-19;                                                             % Nonlinear index (m^2/W)
    %n2 = 1.43*1e-17;%AlGaAs
    pw = (max(W_dat)-min(W_dat))/n_interp;
    W_interp = flip([min(W_dat):pw:max(W_dat)]);                                % Need to be improved?
        
    A_eff_fit = interp1(W_dat,A_eff,W_interp,'spline');
    n2_eff_fit = interp1(W_dat,n2_eff,W_interp,'spline');
    n_eff_fit = interp1(W_dat,n_eff,W_interp,'spline');
    % @ pump wl:
    [u, v] = min(abs(W_interp-w0));% find pump wl
    Aeff1_0 = A_eff_fit(v);% find Aeff @ pump wl
    n2_eff = n2_eff_fit(v);% find n2 @ pump wl
    neff_0 = n_eff_fit(v); % find neff @ pump wl
        
    gamma = w0*n2/(C*Aeff1_0); % calculate nl coeff. @ pump wl
  
        
        der_Aeff1_0 = (A_eff_fit(v+1)-A_eff_fit(v-1))/(W_interp(v+1)-W_interp(v-1));
        der_neff1_0 = (n_eff_fit(v+1)-n_eff_fit(v-1))/(W_interp(v+1)-W_interp(v-1));
        der_n2eff1_0 = (n2_eff_fit(v+1)-n2_eff_fit(v-1))/(W_interp(v+1)-W_interp(v-1));

%     tau_shock = 1/w0-der_Aeff1_0/Aeff1_0-der_neff1_0/neff_0;
    if (isHybridMode == 1)
        
        fprintf('isHybridMode: %d \n',isHybridMode);

        tau_shock = 1/w0-der_Aeff1_0/Aeff1_0-der_neff1_0/neff_0 + der_n2eff1_0/n2_eff;
        gamma = w0*n2_eff/(C*Aeff1_0); % calculate nl coeff. @ pump wl from effective n2 curve
    else
        tau_shock = 1/w0-der_Aeff1_0/Aeff1_0-der_neff1_0/neff_0;
    end


    % for Raman effect, which is negligible for SiN
    
    hr = 0;
    fr = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%% Prepare for linear operator %%%
    % interpolate dispersion to each frequency point.
    W_dat_flip = flip(W_dat);      % W_dat and Beta are sorted with lambda, here we sort them with omega
    [Beta_interp, beta0, beta1, ~] = interp_beta(W_dat_flip,w_interp,B,B1,B2,lambda0); 
    
    
    
    % The orign part of gnlse_YY_ODE45
    %  [Z, AT, AW, W] = gnlse_YY_ODE45(t, AT_0, w0, gamma, Beta_interp,beta0,beta1, ...
    %                                 alpha,fr,hr,tau_shock,flength, nsaves, nplot_seg);
    n = length(t); dT = t(2)-t(1); % grid parameters
    V = 2*pi*(-n/2:n/2-1)'/(n*dT); % frequency grid
    
    %%%%%%%%%%%% L operator %%%%%%%%%%%%%
    L = 1i*(Beta_interp - V*beta1-beta0) - (alpha/2); % Linear operator (dispersion): D=-a/2+i(beta-w*beta_1-beta_0)
    if betas ~= false 
        % betas = false, do simu based on real geo; 
        % betas = [beta2, beta3,...], do simu based on betas
        B = 0;
        betas_V = betas;% copy betas to be betas_V so that it can be used with freq V can changed by DV. betas always keeps the same. 
        if exist('DV','var') ~= false
            betas_V(1) = betas(1).* DV(us(i)*L_tot);% DV: dispersion variation
            fprintf('Theroy study disp vari at %05.1f % cm\n', us(i)*L_tot*100);
            
        else
            betas_V(1) = betas(1);
        end
        
        for ibeta = 1:length(betas_V) %ibeta represent index of each element of betas_V
            B = B + betas_V(ibeta)/factorial(ibeta+1).*V.^(ibeta+1);
        end


        tau_shock = 0 ;
        alpha = 0;% m-1
        % benchmark with paper 'Parabolic pulse gene by use of a disp-decre
        % fiber with ND' doi: 10.1364/OL.29.000498
   

      
       

        L = 1i*(B) - (alpha/2);
    end
    
    L = fftshift(L);% fftshift to freq domain
    
    %%%%%%%%%%%% N operator %%%%%%%%%%%%%
     if abs(w0) > eps               % if w0>0 then include shock
         %gamma = gamma/w0;    
         %W = V + w0;                % for shock W is true freq
         shock_W = 1 + V*tau_shock;   % FT{1+i*d(tau_shock)/dt}
     else
         %W = 1;                     % set W to 1 when no shock
         shock_W = 1;               % set shock term to 1 when no shock
     end
    
    RW = n*ifft(fftshift(RT.'));   % frequency domain Raman
    %W = fftshift(W); % shift to fft space
    shock_W = fftshift(shock_W);    % shift to fft space
    %%%%%%%%%%%% parameter set for fractional Raman %%%%%%%%%%%%
    %RT = 0; % no raman effect
    %fr = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % === define function to return the RHS of Eq. (3.13)

%     function R = rhs(z, AW)
%       AT = fft(AW.*exp(L*z));         % time domain field
%       IT = abs(AT).^2;                % time domain intensity
%       if (length(RT) == 1) || (abs(fr) < eps) % no Raman case
%         M = ifft(AT.*IT);             % response function
%       else
%         RS = dT*fr*fft(ifft(IT).*RW); % Raman convolution
%         M = ifft(AT.*((1-fr).*IT + RS));% response function
%       end
%       %R = 1i*gamma*W.*M.*exp(-L*z);   % full RHS of Eq. (3.13)
%       R = 1i*gamma*shock_W.*M.*exp(-L*z);
%     end
%     
%     
%     % === define function to print ODE integrator status
%     N_flag = 0;% count the flag number
%     function status = report(z, y, flag) % 
%       status = 0;
%       if isempty(flag)
%     %   clc;% clear command windows
%         %fprintf('%05.1f %% complete\n', z/flength*100);
%         N_flag = N_flag + 1;
%         if (mod(N_flag,256) == 0)       % update the status every 256 steps.
%                 fprintf('%05.1f %% complete\n', z/flength*100);
%         end
%       end
%     end

    N_flag = 0;% count the flag number
    % === setup and run the ODE integrator
    Z = linspace(0, flength, nsaves);  % select output z points
    Z = Z';
    % === set error control options
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-12, ...
                     'NormControl', 'on', ...
                     'OutputFcn', @report);
    [Z, AW] = ode45(@rhs, Z, ifft(AT_0), options); % run integrator
    
    
    % === process output of integrator
    if nplot <1
        nplot = 1; % at least 1 frame for each segement
    end
    plot_select = round(linspace(2,nsaves,nplot)); % select only limited simulation results to be ploted to avoid huge graphic demands.
    AW = AW(plot_select,:);
    Z = Z(plot_select);
    AT = zeros(size(AW(1,:)));
    for j = 1:length(AW(:,1))
      AW(j,:) = AW(j,:).*exp(L.'*Z(j)); % change variables
      AT(j,:) = fft(AW(j,:));           % time domain output
      AW(j,:) = fftshift(AW(j,:)).*dT*n;  % scale
      % Explanation of scaling problem:
      % AW is Fourier Transform of AT but realized by ifft in Matlab.
      % ifft in Matlab involved a 1/n * iFFT(AT); iFFT is reverse Fourier transform for
      % an array.
      % iFFT(AT) = ifft(AT) * n
      % A(w) = dT * iFFT(AT) = dT * n * ifft(AT)
    end
    
    W = V + w0; % the absolute frequency grid
    
    % Treatment needed to be done.
    % simulations of each segment needed to be append to each other in order to
    % have the whole temporal-spatial simulation.
    % Z_z, AT_z,AW_z are results concatenated by diff. WGs with diff. widths
    Z_z = cat(1,Z_z,Z+(us(i)-d_us(i))*L_tot);%us(i)-d_us(i): propa. dis. - last segment
    AT_z = cat(1,AT_z,AT);
    AT_0 = AT(end,:);% for loop
    AW_z = cat(1,AW_z,AW);
    fprintf('Propa. Dis.:%d \n',us(i));
    toc;
end
    
 function R = rhs(z, AW)
      AT = fft(AW.*exp(L*z));         % time domain field
      IT = abs(AT).^2;                % time domain intensity
      if (length(RT) == 1) || (abs(fr) < eps) % no Raman case
        M = ifft(AT.*IT);             % response function
      else
        RS = dT*fr*fft(ifft(IT).*RW); % Raman convolution
        M = ifft(AT.*((1-fr).*IT + RS));% response function
      end
      %R = 1i*gamma*W.*M.*exp(-L*z);   % full RHS of Eq. (3.13)
      R = 1i*gamma*shock_W.*M.*exp(-L*z);
    end
    
    
    % === define function to print ODE integrator status
%     N_flag = 0;% count the flag number
    function status = report(z, y, flag) % 
      status = 0;
      if isempty(flag)
    %   clc;% clear command windows
        %fprintf('%05.1f %% complete\n', z/flength*100);
        N_flag = N_flag + 1;
        if (mod(N_flag,256) == 0)       % update the status every 256 steps.
                fprintf('%05.1f %% complete\n', z/flength*100);
        end
      end
    end
end




