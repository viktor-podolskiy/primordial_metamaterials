function [ epsOut] = epsDrudePl( lamArr,lamPl,eps0)
%EPSDRUDE Summary of this function goes here
%   Detailed explanation goes here
c0=3e8; 
omgArr=2e6*pi*c0./lamArr; 

tau=10e13; %[s]
% tau=1e13; %[s]
% tau=1e11; %[s]
% lamPl=8.5; %um
% eps0=12.15; 
omgPl=2e6*pi*c0./lamPl; 
tau=0.1*omgPl; %[s]

epsOut=eps0*(1-omgPl.^2./omgArr./(omgArr+1i*tau)); 


end

