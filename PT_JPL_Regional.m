%==============%
% PT_JPL model %
%==============%
% This program is developed for calculating daily and annual T/ET at regional scale
% Code: Zhongen Niu, IGSNRR, CAS
% Questions to: niuze.16b@igsnrr.ac.cn
% 20200429
%-------------------------------------------------------
clc;clear;close all

%% parameters setting based on land cover
% read land cover dara
[LC,R] = geotiffread('.\parameters\LC.tif');
[m,n] = size(LC);

% constant paramters
k_Rn = 0.6;     %Fisher,2008; Imens&Lemur,1969
b1   = 0.95;    %Ruimy.1999,
b2   = 0.9355;  %Hipps. 1983

%Changed parameters
beta1 = zeros(m,n);
k1 = zeros(m,n);
k2 = zeros(m,n);
beta1(LC==1) = 1.28; k1(LC==1) = 0.57; k2(LC==1) = 0.81;  %Forest
beta1(LC==2) = 1.17; k1(LC==2) = 0.56; k2(LC==2) = 0.91;  %Shrub
beta1(LC==3) = 1.43; k1(LC==3) = 0.59; k2(LC==3) = 0.84;  %Crop
beta1(LC==4) = 0.80; k1(LC==4) = 0.59; k2(LC==4) = 0.80;  %Grassland

Topt = geotiffread('./parameters/Topt.tif');

% parameters to calculate Potential evapotranspiration
alfa =     1.26;
gamma1 =   0.066;

%% Annual loop
for year = 1981:2015
    %initial yearly output results
    year_ETc = zeros(m,n);
    year_ET  = zeros(m,n);
    
    % Calcualte f_APARmax
    f_APARmax = zeros(m,n);
    for mnth = 1:46
        if mnth < 10
            file_LAI = ['.\Drives\LAI_' num2str(year) '_0' num2str(mnth) '.tif'];
        else
            file_LAI = ['.\Drives\LAI_' num2str(year) '_' num2str(mnth) '.tif'];
        end
        LAI = geotiffread(file_LAI);
        f_APAR = b1.*(1.0 - exp(-k1.*LAI));
        hang = f_APARmax<f_APAR;
        f_APARmax(hang) = f_APAR(hang);
    end
    
    %daily loop
    for mnth = 1:46
        %Reading drives
        if mnth < 10
            file_Rn =  ['.\Drives\ras_'  num2str(year) '_0' num2str(mnth) '.tif']; %Net radiation
            file_Ta =  ['.\Drives\TAVG_' num2str(year) '_0' num2str(mnth) '.tif']; %Mean air temperature
            file_RH =  ['.\Drives\RHU_'  num2str(year) '_0' num2str(mnth) '.tif']; %Relative humidity
            file_LAI = ['.\Drives\LAI_'  num2str(year) '_0' num2str(mnth) '.tif']; %Leaf area index
        else
            file_Rn =  ['.\Drives\ras_'  num2str(year) '_' num2str(mnth) '.tif'];  %Net radiation
            file_Ta =  ['.\Drives\TAVG_' num2str(year) '_' num2str(mnth) '.tif'];  %Mean air temperature
            file_RH =  ['.\Drives\RHU_'  num2str(year) '_' num2str(mnth) '.tif'];  %Relative humidity
            file_LAI = ['.\Drives\LAI_'  num2str(year) '_' num2str(mnth) '.tif'];  %Leaf area index
        end
        
        RN =  geotiffread(file_Rn);
        Ta =  geotiffread(file_Ta);
        RH =  geotiffread(file_RH);
        LAI = geotiffread(file_LAI);
        
        %  Unit conversion and quality control of drives
        RN = RN*1000000.0; RN = RN/(24.0*3600.0); RN(RN<0) = 0;
        RH = RH*0.01;  RH(RH<0) = 0;
        Ta = Ta*0.1;
        
        %Intermediate variable
        es = 0.6108.*exp(17.27.*Ta./(Ta+237.3)); % satureation vapur pressure, kPa
        delta = 4098.*es./(Ta+237.3).^2; 
        ea = RH.*es;  % actual vapor pressure, kPa
        VPD = es-ea;  % vapor pressure deficit. kPa
        
        %Net radiation  partitioning
        Rns = RN.*exp(-k_Rn.*LAI);           % Net radiation to the soil  
        Rnc = RN - Rns;                      % Net radiation to the canopy
           
        %parameters restrict potential evaporation to the actual values
        f_APAR = b1.*(1.0 - exp(-k1.*LAI));  %Fraction of PAR absorbed by green vegetation cover
        f_IPAR = b2.*(1.0 - exp(-k2.*LAI));  %Fraction of PAR intercepted by total vegetation cover
        fc = f_IPAR;                         % Fraction total vegetation cover            
        fg = (1.0+f_APAR)./(1.0+f_IPAR);% Green canopy fraction
        ft = exp(-((Ta-Topt)./Topt).^2);      % Plant temperature constraint 
        fm = f_APAR./f_APARmax;              % Plant moisture constraint
        fsm = RH.^(VPD./beta1);              % Soil moisture constraint
        fwet = RH.^4;                        % Relative surface wetness

        %Calculate evapotranspiration
        PJ  = alfa.*delta./(delta+gamma1);
        ETc = (1-fwet).*fg.*ft.*fm.*PJ.*Rnc;
        ETs = (fwet+fsm.*(1-fwet)).*PJ.*Rns;
        ETi = fwet.*PJ.*Rnc;
        
        ETc(ETc<0) = 0;
        ETs(ETs<0) = 0;
        ETi(ETi<0) = 0;
       
        %mm s-1 to mm day
        ETc = ETc*24.0*3600.0/1000000.0/2.44;  %8d之和
        ETs = ETs*24.0*3600.0/1000000.0/2.44;
        ETi = ETi*24.0*3600.0/1000000.0/2.44;
        ET = ETc + ETs + ETi;
        T_ET = ETc./ET;
        T_ET(T_ET<0) = 0; T_ET(T_ET>1)=1;
        
        %output daily results
        doy = (mnth-1)*8+1;
        if doy < 10
            file_out_daily_T = ['.\T\Daily\T_' num2str(year) '_00' num2str(doy) '.tif'];
            file_out_daily_ET = ['.\ET\Daily\ET_' num2str(year) '_00' num2str(doy) '.tif'];
            file_out_daily_T_ET = ['.\T_ET\Daily\T_ET_' num2str(year) '_00' num2str(doy) '.tif'];
        elseif doy < 100
            file_out_daily_T = ['.\T\Daily\T_' num2str(year) '_0' num2str(doy) '.tif'];
            file_out_daily_ET = ['.\ET\Daily\ET_' num2str(year) '_0' num2str(doy) '.tif'];
            file_out_daily_T_ET = ['.\T_ET\Daily\T_ET_' num2str(year) '_0' num2str(doy) '.tif'];
        else
            file_out_daily_T = ['.\T\Daily\T_' num2str(year) '_' num2str(doy) '.tif'];
            file_out_daily_ET = ['.\ET\Daily\ET_' num2str(year) '_' num2str(doy) '.tif'];
            file_out_daily_T_ET = ['.\T_ET\Daily\T_ET_' num2str(year) '_' num2str(doy) '.tif'];
        end
        
        geotiffwrite(file_out_daily_T, ETc,R)
        geotiffwrite(file_out_daily_ET, ET,R)
        geotiffwrite(file_out_daily_T_ET, T_ET,R)
        
        %Calculate year results
        year_ETc = year_ETc + ETc*8;
        year_ET  = year_ET  + ET*8;
    end
    %Output year results
    year_T_ET = year_ETc./year_ET;
    year_T_ET(year_T_ET<0)=0;year_T_ET(year_T_ET>1)=1;
    
    file_out_year_T = ['.\T\Year\year_T_' num2str(year) '.tif'];
    file_out_year_ET = ['.\ET\Year\year_ET_' num2str(year) '.tif'];
    file_out_year_T_ET = ['.\T_ET\Year\year_T_ET_' num2str(year) '.tif'];

    geotiffwrite(file_out_year_T,year_ETc,R)
    geotiffwrite(file_out_year_ET,year_ET,R)
    geotiffwrite(file_out_year_T_ET,year_T_ET,R)
end
