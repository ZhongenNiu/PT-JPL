%==============%
% PT_JPL model %
%==============%
% This program is developed for calculating daily and annual T/ET at site scale
% Code: Zhongen Niu, IGSNRR, CAS
% Questions to: niuze.16b@igsnrr.ac.cn
% 20200429
%-------------------------------------------------------

function [LE_year,LEc_year,T_ET_year,LE,LEc,T_ET] = PT_JPL_Sites_MDF(in, para)
%Output
% LE_year: Annual evapotranspiration 
% LEc_year: Annual transpiration
% T_ET_year: Annual T/ET
% LE: Daily evapotranspiration
% LEc: Daily transpiration
% T_ET: Daily T/ET

%-----------
% Input:
% in(:,1) :: Rn -- Net Radiation (W m-2)
% in(:,2) :: Ta -- Air temperature (Celsius)
% in(:,3) :: RH -- Relative humidity (0.01)
% in(:,4) :: LAI--leaf area index
%-----------
% Forcing Data:
    Rn   = in(:,1); 
    Ta   = in(:,2); 
    RH   = in(:,3);
    LAI  = in(:,4);

% Changed Parameters:
    
    %b1   = para(1,1);      
    k1   = para(1,1);     
    %b2   = para(1,3);      
    k2   = para(1,2);        
    Topt = para(1,3);
    beta = para(1,4); 
    %k_Rn = para(1,7);    


% Constant Parameters:
    k_Rn = 0.6;       % Fisher,2008; Imens&Lemur,1969
    b1 = 0.95;         %Ruimy.1999, GCB
    b2 = 0.9355;       %Hipps. 1983
    %k1 = 0.5;          %Ruimy.1999, GCB
    %k2 = 0.91;         %Hipps. 1983
    
    %Topt = 25;
    %beta = 1;         % Fisher,2008; (kPa)
    
%   parameters to calculate Potential evapotranspiration
    alfa = 1.26;      % Fisher,2008; Priestley&Taylor,1972
    gamma = 0.066;    % psychrometric constant,(kPa C-1)

    % Intermediate variable
    es(:,1) = 0.6108.*exp(17.27.*Ta./(Ta+237.3)); % satureation vapur pressure, kPa
    delta = 4098.*es(:,1)./(Ta+237.3).^2; 
    ea(:,1) = RH.*es(:,1); % actual vapor pressure, kPa
    VPD(:,1) = es(:,1)-ea(:,1); % vapor pressure deficit. kPa

    % Cal Topt
    % TTT = Rn.*Ta.*EVI./VPD;
    % TTT(VPD==0,:) = [];
    % Topt = Ta(TTT == max(TTT),1); % Optmum growth temperature, Celsius
    % Topt=25;         % Garcia,2013; Yao,2013

    % Main program
    
    %changed by Niu Zhongen 20180319
    f_APAR = b1.*(1.0 - exp(-k1.*LAI));  %Fraction of PAR absorbed by green vegetation cover 
    f_IPAR = b2.*(1.0 - exp(-k2.*LAI));  % Fraction of PAR intercepted by total vegetation cover
    %end of changed
    
    % Cal f_APAR during the study time (Yearly or..)
    f_APARmax = max(f_APAR);
    fc = f_IPAR;               % Fraction total vegetation cover
    
    
    Rns = Rn.*exp(-k_Rn.*LAI); % Net radiation to the soil
    Rnc = Rn-Rns;              % Net radiation to the canopy
    fg = f_APAR./f_IPAR;       % Green canopy fraction
    ft = exp(-((Ta-Topt)./Topt).^2); % Plant temperature constraint
    fm = f_APAR./f_APARmax;    % Plant moisture constraint
    fsm = RH.^((VPD(:,1))./beta); % Soil moisture constraint
    fwet = RH.^4;              % Relative surface wetness

    % function cal
    [LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,fg,ft,fm,fsm,Rnc,Rns,G);
    LEc = LEc*24.0*3600.0/1000000.0/2.44; 
    LEs = LEs*24.0*3600.0/1000000.0/2.44;
    LEi = LEi*24.0*3600.0/1000000.0/2.44; 
    LE = LEc+LEs+LEi;
    T_ET = LEc./LE;
    

    LE_year=[];
    for year = 1:length(LE)/46
        temp1 = (year-1)*46+1;
        temp2 = year*46;
        LE_year = [LE_year;sum(LE(temp1:temp2))*8];
        LEc_year = [LEc_year;sum(LE(temp1:temp2))*8];
    end
    T_ET_year = LEc_year./LE_year;
    
end

function [LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,fg,ft,fm,fsm,Rnc,Rns,G)
    PJ  = alfa.*delta./(delta+gamma);
    LEc = (1-fwet).*fg.*ft.*fm.*PJ.*Rnc;
    LEs = (fwet+fsm.*(1-fwet)).*PJ.*(Rns-G);
    LEi = fwet.*PJ.*Rnc;
end
