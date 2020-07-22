%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------%
% Differential Evolution Markov Chain (DEMC) %
%--------------------------------------------%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

% Code: Kun Zhang, Gaofeng Zhu, Lanzhou University, China
% Date:15/12/2016, Modify: 31/05/2017
% Questions to: zhangk12@lzu.edu.cn | kun.zhang322@gmail.com

% Reference:
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [nepa, optpara] = DEMCmain(Nss, crange, modelin, refout, N)
% ---------------
% function input:
% Nss:        the length of loop times
% crange:     prior distribution for the selected parameters, e.g. [1,1,1; 5,5,5]
% modelin:    model forcing data
% refout:     model reference output data
% N:          the numbers of Chains
% ---------------
% function output:
% nepa:       the optimized parameters for each loop
% optpara:    the final optimized parameters
% ---------------
% if nargin < 6
%     N = 30; % Set numbers of Chains
% end


cmin = crange(1,:);
cmax = crange(2,:);
in = modelin;
obs = refout;

% main program
% cont = 0;
d = size(cmin,2);
b = 1e-6;
xo = zeros(N,d);
for i = 1 : d
    for j = 1 : N
        xo(j,i) = cmin(i) + rand .* (cmax(i) - cmin(i)); % Initial paramerters
    end
end
pa = zeros(N,d,Nss);  
pa(:,:,1) = xo;
J_old = Nss .* ones(N,1);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~===============
h = waitbar(0,'Please wait, preparing calculation...');
for s = 1 : Nss
    if mod(s, 10) == 0
        c = 1;
    else
        c = 2.38 ./ sqrt(2 .* d);
    end
    
    % Choose two chains
    for i = 1 : N
        while true
            R1 = floor(rand .* N);
            if R1 ~= i && R1 ~= 0
                break;
            end
        end
        while true
            R2 = floor(rand.*N);
            if R2 ~= R1 && R2 ~= i && R2 ~= 0
                break;
            end
        end
        x_p = zeros(1,d);
        for j = 1 : d
            x_p(j) = xo(i,j) + c .* (xo(R1,j)-xo(R2,j)) + sqrt(b) .* randn; 
            if x_p(j) > cmax(j) || x_p(j) < cmin(j)         
                x_p(j) = cmin(j) + rand .* (cmax(j) - cmin(j));
            end
        end
        
        [result] = PT_JPL_Sites_MDF(in,x_p);

        
        %似然函数--------------------------
        e1 = (norm(result - obs)) .^ 2;  %
        DJ1 = 2 .* std(obs) .^ 2;
        J_new = e1 ./ DJ1;
        delta_J = J_new - J_old(i,1);
        if min(1,exp(-delta_J)) > rand  %MH法则
            xo(i,:) = x_p;
            J_old(i,1) = J_new;
        end
    end
    
    % save the parameters for each loop
    pa(:,:,s) = xo;
    
    cc = round (100 .* (s / Nss));
    str = ['Opt-ing...Running...',num2str(cc),'%'];
    waitbar(s / Nss, h, str);
end
waitbar(1, h, 'Done!');
close(h);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Get result
nepa = permute(pa,[3,1,2]);
optpara = zeros(1,d);
for ii = 1 : d
    xx = reshape(nepa(:,:,ii),Nss*N,1);
    muhat = median(xx);
    optpara(ii) = muhat;
end
end