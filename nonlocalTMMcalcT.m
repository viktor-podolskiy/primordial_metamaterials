function [Ttot,kzi,RMati,TMati] = nonlocalTMMcalcT(omg0,kx,di,epsPerp,epsZZ,alpZZ)
%implementation of transmission/reflection calculations in nonlocal
%multilayered composites; see "Primordial Media: the shrouded realm of
%composite materials" for details


% setting up material-dependent propagation constants, as well as T and R
% matricex
kzi=cell(length(epsPerp),1); 
RMati=cell(length(epsPerp)-1,1); 
TMati=RMati; 

% initialize coefficients/propagation constants
for ilr=1:length(epsPerp)
    if alpZZ(ilr)==0 
        kzi{ilr}=sqrt(epsPerp(ilr)*(omg0^2-kx^2/epsZZ(ilr))); 
        if imag(kzi{ilr})<0
            kzi{ilr}=-kzi{ilr}; 
        end 
        RMati{ilr}=0; 
    else 
        DD=sqrt((alpZZ(ilr)*epsPerp(ilr)+epsZZ(ilr))^2-4*epsPerp(ilr)*alpZZ(ilr)*kx^2/omg0^2); 
        pm=sqrt(epsZZ(ilr)^2)/epsZZ(ilr); 
        kzz=omg0*sqrt((alpZZ(ilr)*epsPerp(ilr)-epsZZ(ilr)+pm*[DD;-DD])/2/alpZZ(ilr));
        kzz(imag(kzz)<0)=-kzz(imag(kzz)<0); 
        kzi{ilr}=kzz; 
        RMati{ilr}=zeros(2); 
    end 
end 

% calculate R/T matrices
Ttot=1; 
for ilr=length(epsPerp)-1:-1:1
    mP=length(kzi{ilr+1}); mM=length(kzi{ilr}); 
    dC=di(ilr+1); 
    if mP==1 && mM==1
        FM=exp(1i*kzi{ilr+1}*dC);
        RTild=FM*RMati{ilr+1}*FM; 
        LMat=[1+RTild, -1; ...
            epsPerp(ilr+1)*omg0/kzi{ilr+1}*(1-RTild), epsPerp(ilr)*omg0/kzi{ilr}]; 
        RMat=[1; epsPerp(ilr)*omg0/kzi{ilr}]; 
        RT=LMat\RMat; 
        TMati{ilr}=RT(1);
        RMati{ilr}=RT(2); 
    elseif mP==2 && mM==1
        FM=diag(exp(1i*kzi{ilr+1}*dC));
        RTild=FM*RMati{ilr+1}*FM;

        %Modified ABC, based on Poynting flux conservation
        NP=[1, 1, 1, 1;...
            epsPerp(ilr+1)*omg0/kzi{ilr+1}(1), epsPerp(ilr+1)*omg0/kzi{ilr+1}(2), ...
            -epsPerp(ilr+1)*omg0/kzi{ilr+1}(1), -epsPerp(ilr+1)*omg0/kzi{ilr+1}(2); ...
            (kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1)), (kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1)), ...
            (kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1)), (kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1))]; 

        NP=NP*[eye(2);RTild]; %NPTilde from the notes
        LMat=[NP,[-1;epsPerp(ilr)*omg0/kzi{ilr};0]]; 
        RMat=[1; epsPerp(ilr)*omg0/kzi{ilr};0]; 

        RT=LMat\RMat; 
        RT=mat2cell(RT,[2 1]); 
        TMati{ilr}=RT{1}; 
        RMati{ilr}=RT{2}; 

        RT=LMat\RMat; 
        RT=mat2cell(RT,[2 1]); 
        TMati{ilr}=RT{1}; 
        RMati{ilr}=RT{2}; 
            
    elseif mP==1 && mM==2 
        FM=exp(1i*kzi{ilr+1}*dC);
        RTild=FM*RMati{ilr+1}*FM; 

        %modified ABC 
        LMat=[1+RTild, -1, -1; ...
            epsPerp(ilr+1)*omg0/kzi{ilr+1}*(1-RTild), ...
            epsPerp(ilr)*omg0/kzi{ilr}(1), epsPerp(ilr)*omg0/kzi{ilr}(2); ...
            0, -(kzi{ilr}(1)^2/omg0^2-epsPerp(ilr)),-(kzi{ilr}(2)^2/omg0^2-epsPerp(ilr))]; 
        RMat=[1 1; ...
            epsPerp(ilr)*omg0/kzi{ilr}(1), epsPerp(ilr)*omg0/kzi{ilr}(2);...
            (kzi{ilr}(1)^2/omg0^2-epsPerp(ilr)), (kzi{ilr}(2)^2/omg0^2-epsPerp(ilr))]; 

        RT=LMat\RMat; 
        RT=mat2cell(RT,[1 2]); 
        TMati{ilr}=RT{1}; 
        RMati{ilr}=RT{2}; 
    elseif mP==2 && mM==2 
        FM=diag(exp(1i*kzi{ilr+1}*dC));
        RTild=FM*RMati{ilr+1}*FM;

        %Modified ABC, based on Poynting flux conservation
        NP=[1, 1, 1, 1;...
            epsPerp(ilr+1)*omg0/kzi{ilr+1}(1), epsPerp(ilr+1)*omg0/kzi{ilr+1}(2), ...
            -epsPerp(ilr+1)*omg0/kzi{ilr+1}(1), -epsPerp(ilr+1)*omg0/kzi{ilr+1}(2); ...
            alpZZ(ilr+1)*(kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1)), alpZZ(ilr+1)*(kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1)), ...
            alpZZ(ilr+1)*(kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1)), alpZZ(ilr+1)*(kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1)); ...
            (kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1))*omg0/kzi{ilr+1}(1), ...
            (kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1))*omg0/kzi{ilr+1}(2), ...
            -(kzi{ilr+1}(1)^2/omg0^2-epsPerp(ilr+1))*omg0/kzi{ilr+1}(1), ...
            -(kzi{ilr+1}(2)^2/omg0^2-epsPerp(ilr+1))*omg0/kzi{ilr+1}(2)]; 

        NP=NP*[eye(2);RTild]; %NPTilde from the notes

        %modified ABC 
        LMat=[NP, ...
            [-1, -1;...
            epsPerp(ilr)*omg0/kzi{ilr}(1), epsPerp(ilr)*omg0/kzi{ilr}(2); ...
            -alpZZ(ilr)*(kzi{ilr}(1)^2/omg0^2-epsPerp(ilr)), -alpZZ(ilr)*(kzi{ilr}(2)^2/omg0^2-epsPerp(ilr)); ...
            (kzi{ilr}(1)^2/omg0^2-epsPerp(ilr))*omg0/kzi{ilr}(1), ...
            (kzi{ilr}(2)^2/omg0^2-epsPerp(ilr))*omg0/kzi{ilr}(2)]];             

        RMat=[1, 1; ... 
            epsPerp(ilr)*omg0/kzi{ilr}(1), epsPerp(ilr)*omg0/kzi{ilr}(2); ...
            alpZZ(ilr)*(kzi{ilr}(1)^2/omg0^2-epsPerp(ilr)), alpZZ(ilr)*(kzi{ilr}(2)^2/omg0^2-epsPerp(ilr)); ...
            (kzi{ilr}(1)^2/omg0^2-epsPerp(ilr))*omg0/kzi{ilr}(1), ...
            (kzi{ilr}(2)^2/omg0^2-epsPerp(ilr))*omg0/kzi{ilr}(2)]; 

        RT=LMat\RMat; 
        RT=mat2cell(RT,[2 2], 2); 
        TMati{ilr}=RT{1}; 
        RMati{ilr}=RT{2}; 
    else 
        disp('unexpected interface')
    end 
    if ilr==length(epsPerp)-1
        Ttot=FM*TMati{ilr}; 
    else
        Ttot=Ttot*FM*TMati{ilr}; 
    end 
end


end

