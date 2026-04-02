function [ epsOut] = epsDrude( lamArr,lamPl,eps0)
%Calculates drude response
c0=3e8; 
omgArr=2e6*pi*c0./lamArr; 

tau=0.8e13; %[Hz]
omgPl=2e6*pi*c0./lamPl; 

epsOut=eps0*(1-omgPl.^2./omgArr./(omgArr+1i*tau)); 


end

