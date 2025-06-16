function [Ex,Hy] = nonlocalTMMfldZ(di,epsPerp,epsZZ,alp,kzi,aPl,aMin,zt,nlOnly)
%Calculates the field distributions based on pre-calculated mode amplitudes

Ex=0*zt; 
Hy=0*zt; 

% calculate positions of the interfaces
zi=0*di; zi(1)=di(1); 
for il=1:length(kzi)-1
    zi(il+1)=zi(il)+di(il+1); 
end 
zi(end)=max(zi(end),zt(end)); 

% zero-out local waves if needed 
if nlOnly
    for il=1:length(kzi)
        aPl{il}(1)=0; 
        aMin{il}(1)=0; 
    end 
end 

for il=length(kzi):-1:1
    if isscalar(kzi{il})
        Fi=exp(1i*kzi{il}*zt); 
        Fi1=exp(-1i*kzi{il}*zt); 
        Exi=1*(Fi*aPl{il}+Fi1*aMin{il}); 
        Hyi=epsPerp(il)/kzi{il}*(Fi*aPl{il}-Fi1*aMin{il}); 
    else 
        Exi=0*zt; Hyi=0*zt; 
        for im=[1,2]
            Fi=exp(1i*kzi{il}(im)*zt); 
            Fi1=exp(-1i*kzi{il}(im)*zt); 
            Exi=Exi+1*(Fi*aPl{il}(im)+Fi1*aMin{il}(im)); 
            Hyi=Hyi+epsPerp(il)/kzi{il}(im)*(Fi*aPl{il}(im)-Fi1*aMin{il}(im));
        end 
    end 
    Ex(zt<=zi(il))=Exi(zt<=zi(il)); 
    Hy(zt<=zi(il))=Hyi(zt<=zi(il)); 
end 

end

