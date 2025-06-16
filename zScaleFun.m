function [scaleData] = zScaleFun(L0,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scaleData=data/L0; 
scaleData(abs(scaleData)>1)=sign(scaleData(abs(scaleData)>1)).*...
    (log10(abs(scaleData(abs(scaleData)>1)))+1); 
scaleData=scaleData*L0; 
end

