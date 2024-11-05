%%% Import files containing dispersion data from mode solver %%%

tic;
isPlot=0;
isPlotGVD = 1;
alpha = log(10)/10*400;
alpha = 0;
C = physconst('LightSpeed');    %Speed of light

isAND_DW = 0;
if isAND_DW == 1
    fprintf('WARNING!! all GVD curve are modified for all normal dispersion DW generation!\n')
end


%% Load all data from simulations.

isEff = 0;
if isEff == 1
    fprintf('WARNING!! effective geo is on, width is not real value!\n')
else 
    FolderName = 'SCG simulation JH';
    SampleName = [600;1200;2000]; % sample variation
end

isTM = 0;

file_name_list=strings();

for i = 1:length(SampleName)
    if isTM == 0
        file_name_list(i) = strcat('SiN_Strip_',num2str(SampleName(i)),'nm_h_700nm_TE_240730.txt');
      
    else
        file_name_list(i) = strcat('SiN_Strip_',num2str(SampleName(i)),'nm_h_700nm_TM_240726.txt');
        
        fprintf('WARNING!! for TM mode!\n')
    
    end
    
    if isEff ==1
        %file_name_list(i) = 'geo_eff_OHY8P13_1100S_9pts.txt';
        fprintf(strcat(file_name_list(i),'\n'))
    end
end

Beta_list = zeros(0);
W_dat_list = zeros(0);
Lambda_list = zeros(0);
Beta1_list = zeros(0);
Beta2_list = zeros(0);
n_eff_list = zeros(0);
A_eff_list = zeros(0);
n2_eff = zeros(0);

for i = 1:length(file_name_list)
   
    file_name = 'width 640nm D 11.txt'; %같은 폴더의 txt를 열어야댐.
    fprintf(strcat('loading.. ',file_name,'\n'));
    Mb = importdata(file_name);   % 시뮬레이션에서 가져온 웨이브가이드 속성
    Mb_dat = Mb.data;
    Lambda = Mb_dat(:, 1) * 1e-9;
    u = find(Lambda >= 500e-9);
    Lambda = Mb_dat(u, 1) * 1e-9;
    W_dat = Mb_dat(u, 2);
    Beta = Mb_dat(u, 3);
    Beta1 = Mb_dat(u, 4);
    Beta2 = Mb_dat(u, 5);
    n_eff = Mb_dat(u, 6);
    A_eff = Mb_dat(u, 7);
    n2_eff = Mb_dat(u, 8);

    if isempty(Beta_list)
        % 배열 크기 초기화
        Beta_list = zeros(length(Beta), length(file_name_list));
        W_dat_list = zeros(size(Beta_list));
        Lambda_list = zeros(size(Beta_list));
        Beta1_list = zeros(size(Beta_list));
        Beta2_list = zeros(size(Beta_list));
        n_eff_list = zeros(size(Beta_list));
        A_eff_list = zeros(size(Beta_list));
        n2_eff_list = zeros(size(Beta_list));
    end
    
    % 데이터 할당
    W_dat_list(:, i) = W_dat;
    Lambda_list(:, i) = Lambda;
    Beta_list(:, i) = Beta;
    Beta1_list(:, i) = Beta1;
    Beta2_list(:, i) = Beta2;
    n_eff_list(:, i) = n_eff;
    A_eff_list(:, i) = A_eff;
    n2_eff_list(:, i) = n2_eff;
end

fprintf('All data from simulation has been loaded.\n')
    
hold off

%% Creation all regression functions to generate dummy simulation.
Widths_simu = SampleName; %numer of variattion of width 

Fs_Beta2 = cell(length(W_dat), 1);
for i = 1:length(W_dat)
    % 보간을 위한 최소 두 개의 포인트가 있는지 확인
    if length(Widths_simu) >= 2 && length(Beta2_list(i, :)) >= 2
        F = griddedInterpolant(Widths_simu, Beta2_list(i, :), 'spline');
        Fs_Beta2{i} = F;
    else
        error('Not enough data points for interpolation.');
    end
    
    if isPlot == 1
        figure()
        lgd = strcat('@', num2str(Lambda_list(i) * 1e9), 'nm');
        plot(Widths_simu, Beta2_list(i, :) * 1e24, 'o', 'DisplayName', lgd);
        ylabel(['\beta_{2} ', 'ps^{2} m^{-1}']);
        xlabel('Width (nm)');
        title('TE polarisation');
        legend()
        hold on
        Widths_interp = 700:10:1500;
        plot(Widths_interp, Fs_Beta2{i}(Widths_interp) * 1e24, '-');
    end
end
hold off

Fs_Beta = cell(length(W_dat),1);
for i = 1:length(W_dat)
    F = griddedInterpolant(Widths_simu,Beta_list(i,:),'spline');
    Fs_Beta{i} = F;
end

Fs_Beta1 = cell(length(W_dat),1);
for i = 1:length(W_dat)
    F = griddedInterpolant(Widths_simu,Beta1_list(i,:),'spline');
    Fs_Beta1{i} = F;
end

Fs_n_eff = cell(length(W_dat),1);
for i = 1:length(W_dat)
    F = griddedInterpolant(Widths_simu,n_eff_list(i,:),'spline');
    Fs_n_eff{i} = F;
end

Fs_A_eff = cell(length(W_dat),1);
for i = 1:length(W_dat)
    F = griddedInterpolant(Widths_simu,A_eff_list(i,:),'spline');
    Fs_A_eff{i} = F;
end

Fs_n2_eff = cell(length(W_dat),1);
for i = 1:length(W_dat)
    F = griddedInterpolant(Widths_simu,n2_eff_list(i,:),'spline');
    Fs_n2_eff{i} = F;
end
fprintf('All regression functions to generate dummy simulation data created.\n')


%% Define waveguide geometry, Dispersion managed WG geo creation
% with unified parametry u from 0 to 1
clear ParaFun


L = 0.006;%m

% ParaFun = @(u) 3200+u.*-2000;
ParaFun = @(u) 640+u.*0;

isTailor = 0
isOHY8P13 = 0
isDM_design = 0
betas = -9.661465159831380e-25;
betas= false;

if betas ~= false
% define DV as a function handle
% it should only work when betas has values which means we are in
% theoritical study mode.
% later it will be used as betas .* DV as a function of propagation
% distance us*L
    DV = @DispersionVariation;
else
    DV = @ (x) 1;
end

%  ParaFun = @(u) ((u<=0.3).*3200 + (u>0.3).*(3200+u.*-2000))                                      % in nm
% Create a new parametric function


Width_Tol = 25;                                                               % tolereance in nm
figure()
num_discretization = 2^15;
Lp = linspace(0,1,num_discretization);% propa. dis. discretization.
ws_discretization = ParaFun(Lp);

%%% For periodical structure %%%
    num_period = 1;
    Lp = linspace(0,1,num_discretization*num_period);
    ws_discretization = repmat(ws_discretization,1,num_period);
%%%
% plot(Lp,ws_discretization,'-','LineWidth',0.5,'DisplayName','Geo. origin');
xlabel('Distance m')
ylabel('Width nm')
% title(strcat('Waveguide length:  ',num2str(L*100),'cm',', tolerance:',num2str(Width_Tol),'nm'))
% hold on
% h = 0.02;       % step size
% X = 0:h:1;     % domain
% f = ParaFun(X);% range
% Y = diff(f)/h;   % first derivative
% % Z = diff(Y)/h;   % second derivative
% plot(X(:,1:length(Y)),Y,'r',X,f,'b','LineWidth',1.5)%, X(:,1:length(Z)),Z,'k')
% %%% linear regression for first derivative %%%
% F_dY = griddedInterpolant(X(:,1:length(Y)),Y,'linear');
hold on

%%% calculate steps from u=0 %%%

% New version of discretization of WG geo
d_Lp = Lp(2)-Lp(1);
j = 1;
d_us=[]; d_us(1) = 0;
us = []; us(1) = 0;
while us(j)<1.0
    D_Lp = d_Lp; %initialize by basic step defined by d_Lp
    while abs(ParaFun(us(j))-ParaFun(us(j)+D_Lp))<=Width_Tol
        if us(j)+D_Lp >1
            break;
        end

        D_Lp = D_Lp+d_Lp;
    end
    D_Lp = D_Lp-d_Lp; % go back to tolerance
%     fprintf('Width difference right now is:%d nm.\n',abs(ParaFun(us(j))-ParaFun(us(j)+D_Lp)))
    if D_Lp<=0
        fprintf('Warning!! Cannot fulfill tolerance demand, please change to more indense Lp.\n')
        fprintf('Problem happens at us(j) = %f.\n',us(j))
        break;
    else
        us(j+1) = us(j)+D_Lp;
        d_us(j+1) = D_Lp;
    end
    j = j+1;
    % make sure propagation distance do not exceed 1.
    if us(j)>=1
        us(j) = 1;
        d_us(j) = 1-us(j-1);
    end
    fprintf('Propa Dis u:%d.\n',us(j))
end



%
us(1) =[]; % to remove the first element 0
d_us(1) =[]; % to remove the first element 0
us = us';
d_us = d_us';
widths = ParaFun(us);
%%% consider periodical stucture%%%
    d_us = repmat(d_us,num_period,1);
    d_us = d_us./num_period; %renormalised to 1
    widths = repmat(widths,1,num_period);
    if num_period>1
        us_p = us;
        for k = 2:num_period
            us = cat(1,us,(us_p+k-1));
        end
    end
    us = us./num_period; % renormalised to 1
%%%
% for w= 1400 str wg 0.008m or 0.01m around 700W
% for w= 1100 str wg 0.006m for power around 700W 
% for w= 1100 str wg 0.015m for power 200W @1310nm
% 
if isTailor == 1
%     L = 1152;%m
%     us = linspace(0,L,257)';
    L = 0.02;%m
    us = linspace(0,L,129)';
    d_us = diff(us);
    us(1) = [];
    widths = 999.*ones(size(us));% no real waveguide widths
    us = us/L;
    d_us = d_us/L;
end

% 
if isTailor == 1
% %     simuation 2 segment: 1 str and 1 tapered WG for soliton shifting(dispersion decreasing design).
    L1 = 40*5*1e-6;%m
    L2 = 10*5*1e-6;
    W1 = 1000;%nm
    W2 = 800;
    
    
    L1 = 50*1e-6;%m
    L2 = 50*1e-6;
    W1 = 1700;%nm
    W2 = 2300;
    
    L1 = 50*1e-6;%m
    L2 = 50*1e-6;
    W1 = 1700;%nm
    W2 = 1000;
    W1 = 1400;%nm
    W2 = 1000;
 
    
    W1 = 1000;%nm
    W2 = 1400;
    d_us = [];
    widths = [];
    for i = 1:200
        d_us = cat(1,d_us,L1,L2);
        widths = cat(1,widths,W1,W2);
    end
%     d_us = d_us';
%     widths = widths';
    us = cumsum(d_us);
    L = us(end);
    us = us/L;
    d_us = d_us/L;

end

if isTailor == 1


    Period = 50 * 1e-6;% period of QPM waveguide
    Ampli = 250;% nm
    Bias = 1450;% nm
    
    Period = 50*5 * 1e-6;% period of QPM waveguide
    Ampli = 200;% nm
    Bias = 1200;% nm
    
    d_us = Period/8 * ones(160*4,1);% each sampling is 1/8 of the period
    us = cumsum(d_us);
    widths = Ampli* sin(2*pi/Period * us) + Bias;
    L = us(end);
    us = us/L;
    d_us = d_us/L;


end


% if isTailor == 1
% %     L = 0.06;
%     us = cat(1,linspace(0.01,0.05,64)');
%     d_us = cat(1,0.008,diff(us));
%     widths_interm = us;
%     widths_interm(1) = [];
% %     widths = cat(1,1400,widths_interm.*(-700)./(L-0.008)+1550);
%     widths = cat(1,linspace(1400,900,length(us))');
%     us = cat(1,us,L);
%     d_us = cat(1,d_us,L-0.05);
%     widths = cat(1,widths,2500);
%     us = us/L;
%     d_us = d_us/L;
% end
if isOHY8P13 == 1
    % try to simulation the quais periodical width modulation of OHY8P13
    [propa_dis,widths,L] = ParaFunCreator_us();
%     
%     if L>0.015*1e6
%         L = 0.015*1e6;
%         widths(propa_dis>0.015*1e6)=[];
%         propa_dis(propa_dis>0.015*1e6)=[];
%     end
%     
    us = propa_dis./L;
    d_us = diff(us);
    widths(1) = [];
    us(1) = [];
    L = L*1e-6;% convert to m
end
if isDM_design == 1
    [propa_dis,widths,L] = DM_design();
    us = propa_dis;
    d_us = diff(us);
    widths(1) = [];
    us(1) = [];

end
title(strcat('Waveguide length:  ',num2str(L*100),'cm'))%,', tolerance:',num2str(Width_Tol),'nm'))
plot(us.*L,widths,'-o','LineWidth',1,'DisplayName','Geo. simu.')
xlim([0 L])
fprintf('%d segements at tolerence of %d nm.\n',length(d_us),Width_Tol);
fprintf('Dispersion managed WG geo created.\n')
legend()
p_Geo = gcf;
hold off

%% find dummy simulation data at different given width.
Widths_Dummy = widths;
Width_num = length(Widths_Dummy);
W_num = length(W_dat);
Betas_Dummy = zeros(W_num,Width_num);
Beta1s_Dummy = zeros(W_num,Width_num);
Beta2s_Dummy = zeros(W_num,Width_num);
n_effs_Dummy = zeros(W_num,Width_num);
A_effs_Dummy = zeros(W_num,Width_num);
n2_effs_Dummy = zeros(W_num,Width_num);

for i = 1:length(Widths_Dummy)
    Width_Dummy = Widths_Dummy(i);
    [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy] = Dummy_simu(Width_Dummy, ...
        W_dat,Lambda,Fs_Beta,Fs_Beta1,Fs_Beta2,Fs_n_eff,Fs_A_eff,Fs_n2_eff, ...
        isPlot);
    Betas_Dummy(:,i) = Beta_Dummy;
    Beta1s_Dummy(:,i) = Beta1_Dummy;
    Beta2s_Dummy(:,i) = Beta2_Dummy;
    n_effs_Dummy(:,i) = n_eff_Dummy;
    A_effs_Dummy(:,i) = A_eff_Dummy;
    n2_effs_Dummy(:,i) = n2_eff_Dummy;

end
% % the effective parameter of a str WG equivlent to real geo
% Omega_eff = 2*pi.*C./Lambda;
% Beta_eff = Betas_Dummy*d_us;
% Beta1_eff = Beta1s_Dummy*d_us;
% Beta2_eff=Beta2s_Dummy*d_us;%effective beta2, the average GVD of all waveguide
% n_eff_eff = n_effs_Dummy*d_us;
% A_eff_eff = A_effs_Dummy*d_us;
% ng_eff = Beta_eff*0;
% n2_eff_eff = Beta_eff*0;
% D_eff = Beta_eff*0;
% 
% Table_geo_eff = table(Lambda*1e9,Omega_eff,Beta_eff,Beta1_eff,Beta2_eff,n_eff_eff,A_eff_eff,ng_eff,n2_eff_eff,D_eff, ...
%                     'VariableNames',{'Wavelength(nm)','Omega(rad/s)','Beta0(1/m)',...
%                     'Beta1(s/m)','Beta2(s^2/m)','neff','Aeff(m^2)','ng','n2_eff','D(ps/(nm*km)'});
% writetable(Table_geo_eff,strcat(SaveName,'/geo_eff','.txt'),'Delimiter','\t');
% fprintf('Simu. at %d saved in the path: %s \n',0,SaveName);
hold off
fprintf('All dummy simu data created.\n')

% 
% Beta_Dummy = zeros(length(W_dat),1);
% Beta1_Dummy = zeros(length(W_dat),1);
% Beta2_Dummy = zeros(length(W_dat),1);
% n_eff_Dummy = zeros(length(W_dat),1);
% A_eff_Dummy = zeros(length(W_dat),1);
% 
% 
% for i = 1:length(Fs_Beta2)
%     Beta2_Dummy(i) = Fs_Beta2{i}(Width_Dummy);
%     Beta_Dummy(i) = Fs_Beta{i}(Width_Dummy);
%     Beta1_Dummy(i) = Fs_Beta1{i}(Width_Dummy);
%     n_eff_Dummy(i) = Fs_n_eff{i}(Width_Dummy);
%     A_eff_Dummy(i) = Fs_A_eff{i}(Width_Dummy);
% end
if isPlot == 1
    figure(5)
    lgd = strcat('width ',num2str(Width_Dummy),'nm');
    plot(Lambda*1e9,Beta2_Dummy*1e24,'o','DisplayName',lgd);
    ylabel(['\beta_{2} ','ps^{2} m^{-1}']);
    xlabel('Wavelength (nm)')
    title('TE polarisation');
    legend()
    hold on
    plot(Lambda*1e9,0.*Lambda*1e9,'Color','black','LineWidth',1.5,'LineStyle','--');
    % compare to real simulaion results FDTD
    Mb = importdata('Dispersion_file_Lumerical_ModeSolver/OHY8P13_1000nm_height_800nm_TE_221011.txt');        % Waveguide properties from simulation
    Mb_dat = Mb.data;
    Lambda = Mb_dat(:,1)*1e-9;
    u = find(Lambda>=500e-9);
    Lambda = Mb_dat(u,1)*1e-9;
    W_dat = Mb_dat(u,2);
    Beta = Mb_dat(u,3);
    Beta1 = Mb_dat(u,4);
    Beta2 = Mb_dat(u,5);
    plot(Lambda*1e9,Beta2*1e24,'r','DisplayName','simulation 1000nm')
    hold on
end
toc;
%% functions
function [Beta_Dummy,Beta1_Dummy,Beta2_Dummy,n_eff_Dummy,A_eff_Dummy,n2_eff_Dummy] = Dummy_simu(Width_Dummy, ...
    W_dat,Lambda,Fs_Beta,Fs_Beta1,Fs_Beta2,Fs_n_eff,Fs_A_eff,Fs_n2_eff, ...
    isPlot)
    simu_num = length(W_dat);
    Beta_Dummy = zeros(simu_num,1);
    Beta1_Dummy = zeros(simu_num,1);
    Beta2_Dummy = zeros(simu_num,1);
    n_eff_Dummy = zeros(simu_num,1);
    A_eff_Dummy = zeros(simu_num,1);
    n2_eff_Dummy = zeros(simu_num,1);
    for i = 1:length(Fs_Beta2)
        Beta2_Dummy(i) = Fs_Beta2{i}(Width_Dummy);
        Beta_Dummy(i) = Fs_Beta{i}(Width_Dummy);
        Beta1_Dummy(i) = Fs_Beta1{i}(Width_Dummy);
        n_eff_Dummy(i) = Fs_n_eff{i}(Width_Dummy);
        A_eff_Dummy(i) = Fs_A_eff{i}(Width_Dummy);
        n2_eff_Dummy(i) = Fs_n2_eff{i}(Width_Dummy);

    end
    if isPlot == 1
        figure(5)
        lgd = strcat('width ',num2str(Width_Dummy),'nm');
        plot(Lambda*1e9,Beta2_Dummy*1e24,'o','DisplayName',lgd);
        ylabel(['\beta_{2} ','ps^{2} m^{-1}']);
        xlabel('Wavelength (nm)')
        title('TE polarisation');
        legend()
        hold on
        plot(Lambda*1e9,0.*Lambda*1e9,'Color','black','LineWidth',1.5,'LineStyle','--');
%         lgd_ref = legend(ref);
%         lgd_ref = legend('AutoUpdate','off'); 
    end
end

%% functions
%%

function D_z = DispersionVariation(z)
%     g = 0.028; %m-1 from 10.1364/OL.29.000498 benchmark
    g = 46; %m-1
    g = 125;
    g=0;
    D_z = 1./(1 + g.*z);
end
