function [scaleData] = zScaleFun(L0,data)
%Nonlinear scale for 3D dispersion plots
scaleData=data/L0; 
scaleData(abs(scaleData)>1)=sign(scaleData(abs(scaleData)>1)).*...
    (log10(abs(scaleData(abs(scaleData)>1)))+1); 
scaleData=scaleData*L0; 
end

