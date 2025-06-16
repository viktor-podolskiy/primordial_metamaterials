function [ epsOut] = epsDrude( lamArr,lamPl,eps0)
%EPSDRUDE Drude permittivity
%   expects array of wavelength (in microns), plasma wavelength, and
%   background permittivity; 
%   outputs Drude mode predictions


c0=3e8; %speed of light in vacuum
omgArr=2e6*pi*c0./lamArr; % frequency

tau=0.8e13; %[s]
omgPl=2e6*pi*c0./lamPl; % plasma frequency

epsOut=eps0*(1-omgPl.^2./omgArr./(omgArr+1i*tau)); 

end

