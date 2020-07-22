%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------%
%  Sobol' Global sensitivity analysis method %
%--------------------------------------------%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

function [Si,St] = Sobol_SA(interval, in, refout, Nss)
% ---------------
% function input:
%interval:    prior distribution for the selected parameters
%in:          model forcing data
%refout:      model reference output data
% Nss:        the length of loop times

% ---------------
% function output:
%Si:          First-order sensitivity index
%St:          Total-order sensitivity index
% ---------------

interval = interval';
obs = refout;
dem=size(interval,1); % dimension of parameter vector
% Generate matrices A and B using the LHS technique
[Am,Bm]=LHSsample(Nss,dem,interval);                
% calculate the sensitivity index
sa1=0;
sb1=0;
ss1=zeros(1,dem);
st1=zeros(1,dem);
% AB=zeros(1,dem);
h1=waitbar(0,'Please Wait...Have a cup of Tea...');
% preallocating
bfo = zeros(Nss,1);
VY  = bfo;
Si  = zeros(Nss, dem);
St  = Si;
for j=1:Nss
    

    [y_s1] = PT_JPL_Sites_MDF(in,Am(j,:));
    
    out1=abs(real(RMSE(obs,y_s1)));
    sa1=sa1+out1;
    sb1=sb1+out1^2;
    bfo(j,1)=sa1/j;
    VY(j,1)=sb1/j-bfo(j,1)^2;
    
    [y_s2] = PT_JPL_Sites_MDF(in,Bm(j,:));  
    
    
    out2=abs(real(RMSE(obs,y_s2)));
    for i=1:dem
        AB=Am(j,:);
        AB(i)=Bm(j,i);
        
        [y_s3] = PT_JPL_Sites_MDF(in,AB(1,:)); %
        
        
        out3=abs(real(RMSE(obs,y_s3)));
        ss1(i)=out2*(out3-out1)+ss1(i);
        Si(j,i)=ss1(i)/(j*VY(j,1));
        st1(i)=st1(i)+(out3-out1)^2;
        St(j,i)=st1(i)/(2*j*VY(j,1));
    end
    waitbar(j/Nss,h1,'Please Wait...Have a cup of Tea...');
end
waitbar(1,h1,'Done!!');
delete(h1);
end
%===============================================================
function[A,B]=LHSsample(N,d,interval)

interval=[interval;interval];
%Generates a LHS M1 containing N samples and 2*d dimension
M1=lhsdesign(N,2*d);
M = zeros(size(M1,1),size(M1,2));
for j=1:size(M1,2)
    int=interval(j,:);
    for i=1:size(M1,1)
        %transform to parameter space
        M(i,j)=unifinv(M1(i,j),int(1),int(2));
    end
end
A=M(:,1:d);             % the first d columns were designed to matrix A
B=M(:,d+1:end);        % the last d columns were designed to matrix B
end
%===============================================================
function [RMSE] = RMSE(y_obs,y_s)
RMSE=sqrt(sum((y_obs-y_s).^2)./length(y_obs));
end






