%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------%
%            LAI  Preparation                %
%--------------------------------------------%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

% Interpolation of GLOBMAP LAI data during 1981 and 2000 (half-month to 8 day)
clc;clear;close all
%Calculate time series leaf area index (1~365)
HalfMonth_1 = {...
    '001','016','032','047','060','075','091','106','121','136',...
    '152','167','182','197','213','228','244','259','274','289',...
    '305','320','335','350'};
HalfMonth_2 = {...
    '001','016','032','047','061','076','092','107','122','137',...
    '153','168','183','198','214','229','245','260','275','290',...
    '306','321','336','351'};
for year = 1981:2000
    for num = 1:length(HalfMonth_1)-1
        %read LAI data
        if mod(year,4) == 0
            file = ['.\' num2str(year) '\GlobMapLAIV1.A'  num2str(year) HalfMonth_2{num} '.Global.tif']
           [data,R] = geotiffread(file);
            for doy = str2double(HalfMonth_2{num}):str2double(HalfMonth_2{num+1})-1
                if doy < 10
                    file_out =  ['.\TimeSeries\' num2str(year) '00' num2str(doy) '.tif']
                elseif doy < 100
                    file_out =  ['.\TimeSeries\' num2str(year) '0' num2str(doy) '.tif']
                else
                    file_out =  ['.\TimeSeries\' num2str(year) '' num2str(doy) '.tif']
                end
                geotiffwrite(file_out,data,R)
                if doy == str2double(HalfMonth_2{end})-1
                    for doy = str2double(HalfMonth_2{end}):365
                        file_out =  ['.\TimeSeries\' num2str(year) '' num2str(doy) '.tif']
                        geotiffwrite(file_out,data,R)
                    end
                end
                
            end
            doy;
        else
            file = ['.\' num2str(year) '\GlobMapLAIV1.A'  num2str(year) HalfMonth_1{num} '.Global.tif']
            [data,R] = geotiffread(file);
            for doy = str2double(HalfMonth_1{num}):str2double(HalfMonth_1{num+1})-1
                if doy < 10
                    file_out =  ['.\TimeSeries\' num2str(year) '00' num2str(doy) '.tif']
                elseif doy < 100
                    file_out =  ['.\TimeSeries\' num2str(year) '0' num2str(doy) '.tif']
                else
                    file_out =  ['.\TimeSeries\' num2str(year) '' num2str(doy) '.tif']
                end
                geotiffwrite(file_out,data,R)
            end
             if doy == str2double(HalfMonth_1{end})-1
                 for doy = str2double(HalfMonth_1{end}):365
                     file_out =  ['.\TimeSeries\' num2str(year) '' num2str(doy) '.tif']
                     geotiffwrite(file_out,data,R)
                 end
             end
        end

    end
end

%Calulate 8 day average leaf area index
clc;clear;close all
day8 = {...
    '001','009','017','025','033','041','049','057','065','073',...
    '081','089','097','105','113','121','129','137','145','153',...
    '161','169','177','185','193','201','209','217','225','233',...
    '241','249','257','265','273','281','289','297','305','313',...
    '321','329','337','345','353','361'};
for year = 1981:2000
    for num = 1:46
        data_out = zeros(2091,4950);
        ii = 1;
        for doy = str2double(day8{num}):str2double(day8{num})+7
            if doy == 366
                break;
            end
            if doy < 10
                file_in = ['.\TimeSeries\' num2str(year) '00' num2str(doy) '.tif']
            elseif doy < 100
                file_in = ['.\TimeSeries\' num2str(year) '0' num2str(doy) '.tif']
            else
                file_in = ['.\TimeSeries\' num2str(year) '' num2str(doy) '.tif']
            end
            [data,R] = geotiffread(file_in);
            data_out = data_out + double(data);
            ii = ii+1;
        end
        data_out = data_out/ii;
        file_out = ['.\GLOBMAP\GLOBMAP_' num2str(year) '_' day8{num} '.tif']
        geotiffwrite(file_out,data_out,R);
        data_out = []; ii = 1;
    end
end

%% Arithmetic mean value of GLOBMAO LAI and GLASS LAI during 1981 and 2015
clc;clear;close all
for year = 1981:2015
   for mnth = 1:46 
       doy = (mnth-1)*8+1
       % files name
       if doy < 10
           file_GLOBMAP = ['.\GLOBMAP\GLOBMAP_' num2str(year) '_00' num2str(mnth) '.tif'];
           file_GLASS   = ['.\GLASS\GLASS_' num2str(year) '_00' num2str(mnth) '.tif'];
       elseif doy < 100
           file_GLOBMAP = ['.\GLOBMAP\GLOBMAP_' num2str(year) '_0' num2str(mnth) '.tif'];
           file_GLASS   = ['.\GLASS\GLASS_' num2str(year) '_0' num2str(mnth) '.tif'];
       else
           file_GLOBMAP = ['.\GLOBMAP\GLOBMAP_' num2str(year) '_' num2str(mnth) '.tif'];
           file_GLASS   = ['.\GLASS\GLASS_' num2str(year) '_' num2str(mnth) '.tif'];
       end
       
       % read lai data
       [GLOBMAO,R] = geotiffread(file_GLOBMAO);
       [GLASS,R]   = geotiffread(file_GLASS);
       
       % Arithmetic mean value
       data_out = (GLOBMAO + GLASS)./2;
       
       % output data
       file_out = ['./meanLAI/LAI_' num2str(year) '_' num2str(mnth) '.tif'];
       geotiffwrite(file_out,data_out,R)
       data_out = [];
   end
end
