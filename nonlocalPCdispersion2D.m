%-----------------------------
% Intrinsically nonlocal metamaterials manuscript 
% Supplementary Information
% The script calculates the figure describing 2D dispersion plot of 
% nonlocal composites (Fig.S5e). 
%
% copyright (2025) U Mass Lowell 
% Jacob LaMountain and Viktor Podolskiy 
%-----------------------------


clear 

nxArr=[1.]; % x-component of dimensionless wavevector

% 2D/wavelength plots
lamArr=(4:0.0125:10); 

[lamArr,nxArr]=meshgrid(lamArr,nxArr); 

omgArr=2*pi./lamArr; 

% permittivity and nonlocality parameters
eps1Arr=(epsDrude(lamArr,6.3,12.3));
eps2Arr=10+0*eps1Arr; 
alp1Arr=-12.3*4.1e-5*3/5/6.3^2./(1./lamArr.*(1./lamArr+0.015i)); alp1=alp1Arr(1); 
alp2=1e-4*(1+0.3i); 



dd=0.01*[1,1]; % layer thicknesses

ne1Arr=0*nxArr; 
ne2Arr=0*nxArr; 
nzlArr=0*nxArr; 

%% running

for il=1: numel(lamArr)

    nx=nxArr(il); 
    omg0=omgArr(il);

    epsZZ=[eps1Arr(il), eps2Arr(il)];
    epsPerp=epsZZ; 
    alp1=alp1Arr(il);
    alpZZ=[alp1,alp2];



    NN=cell(2,1); 
    Nnorm=cell(2,1); 
    Nninv=cell(2,1); 
    FN=cell(2,1); 
    
    NL=cell(2,1); 
    FL=cell(2,1); 
    kzi=cell(2,1); 
    kzl=cell(1,1); 
    
    for ilr=[1,2]
        DD=sqrt((alpZZ(ilr)*epsPerp(ilr)+epsZZ(ilr))^2-4*epsPerp(ilr)*alpZZ(ilr)*nx^2); 
        nzz=sqrt((alpZZ(ilr)*epsPerp(ilr)-epsZZ(ilr)+[DD;-DD])/2/alpZZ(ilr));
        nzz(imag(nzz)<0)=-nzz(imag(nzz)<0); 
    
        numr=nzz.^2-epsPerp(ilr); 

        kzi{ilr}=omg0*nzz; 
    
        NN{ilr}=[1, 1, 1, 1; ...
            epsPerp(ilr)/nzz(1), epsPerp(ilr)/nzz(2), -epsPerp(ilr)/nzz(1), -epsPerp(ilr)/nzz(2); ...
            numr(1)*alpZZ(ilr), numr(2)*alpZZ(ilr), numr(1)*alpZZ(ilr), numr(2)*alpZZ(ilr); ...
            numr(1)/nzz(1), numr(2)/nzz(2), -numr(1)/nzz(1), -numr(2)/nzz(2)]; 

        FN{ilr}=diag(exp(1i*omg0*[nzz(1), nzz(2), -nzz(1), -nzz(2)]*dd(ilr))); 
    
        nzl=sqrt(epsPerp(ilr)*(1-nx^2/epsZZ(ilr)));
        if(imag(nzl)<0)
            nzl=-nzl; 
        end 
        NL{ilr}=[1, 1; ...
            epsPerp(ilr)/nzl, -epsPerp(ilr)/nzl]; 
        FL{ilr}=diag(exp(1i*omg0*[nzl, -nzl]*dd(ilr))); 
        kzl{ilr}=omg0*nzl; 
    end 

    % nonlocal and local transfer matrices of one period, see Eq.(S6)
    MM=FN{1}*(NN{1}\NN{2})*FN{2}*(NN{2}\NN{1}); 
    ML=FL{1}*(NL{1}\NL{2})*FL{2}*(NL{2}\NL{1}); 

    ne=(log(eig(MM)))/1i/omg0/sum(dd); 
    nl=(log(eig(ML)))/1i/omg0/sum(dd); 
    [~,ind]=sort(abs(ne.^2)); 

    % calculating z-components of wavenumbers
    ne1Arr(il)=ne(ind(2)); 
    ne2Arr(il)=ne(ind(3)); 

    nzlArr(il)=nl(1); 
    
end 

%% plotting

ne1Arr(imag(ne1Arr)<0)=-ne1Arr(imag(ne1Arr)<0); 
ne2Arr(imag(ne2Arr)<0)=-ne2Arr(imag(ne2Arr)<0); 


epsPerpEMT=(eps1Arr+eps2Arr)/2; 
epsZZEMT=2./(1./eps1Arr+1./eps2Arr); 
neEMT=sqrt(epsPerpEMT.*(1-nxArr.^2./epsZZEMT)); 

RGB = orderedcolors("gem");

figure(23)
clf 

t=tiledlayout(1,1); 
ax1 = axes(t);
ax1.XColor = RGB(1,:);


hold on 

plot(ax1,real(ne1Arr(1,:)),omgArr(1,:), '-','LineWidth',2,'Color',RGB(1,:))
plot(ax1, imag(ne1Arr(1,:)),omgArr(1,:), '--','LineWidth',2,'Color',RGB(1,:))

% uncomment lines below to show local EMT calculations
% plot(ax1, real(neEMT(1:25:end)), omgArr(1,(1:25:end)), '^','LineWidth',2,'Color',RGB(1,:))
% plot(ax1, imag(neEMT(1:25:end)), omgArr(1,(1:25:end)), '<','LineWidth',2,'Color',RGB(1,:))

xlim(ax1,1*[-1.5 3])
set(ax1, 'FontSize', 18)

omgTicks=(0.6:0.2:1.6); 
yticks(ax1,omgTicks);
ylabel('\omega/c, \mum^-1')
xlabel('$\tilde{k}_z$','Interpreter','latex')
ax2 = axes(t);

plot(ax2,real(ne2Arr(1,:)),omgArr(1,:),'-','LineWidth',2,'Color',RGB(2,:))
hold on 
plot(ax2, imag(ne2Arr(1,:)),omgArr(1,:),'--','LineWidth',2,'Color',RGB(2,:))
ax2.XColor = RGB(2,:);

xlim(ax2,100*[-1.5 3])

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
set(ax2, 'FontSize', 18)


lamTicks=(10:-2:4); 
yticks(ax2,2*pi./lamTicks)
yticklabels(ax2,lamTicks)
ylabel(ax2,'\lambda_0,\mum')
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
