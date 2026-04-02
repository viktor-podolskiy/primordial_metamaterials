%-----------------------------
% Intrinsically nonlocal metamaterials manuscript 
% Supplementary Information
% The script calculates the figure describing 3D dispersion plot of 
% nonlocal composites (Fig.S6) and field distributions (Fig.3c) 
%
% copyright (2025) U Mass Lowell 
% Jacob LaMountain and Viktor Podolskiy 
%-----------------------------

clear 

% 3D/wavelength plots
% lamArr=(4:0.025:10);
% nxArr=(0.001:0.025:3);

% field distribution plots
lamArr=5.05; 
nxArr=0.5; 

[lamArr,nxArr]=meshgrid(lamArr,nxArr); 
omgArr=2*pi./lamArr; 

% setup permittivity and nonlocality of materials
eps1Arr=(epsDrude(lamArr,6.3,12.3));
eps2Arr=10+0*eps1Arr; 


alp1Arr=-12.3*4.1e-5*3/5/6.3^2./(1./lamArr.*(1./lamArr+0.015i)); alp1=alp1Arr(1); 
alp2=1e-4*(1+0.3i);  


dd=[0.01,0.01]; % layer thicknesses

nz1Arr=0*nxArr; 
nz2Arr=0*nxArr; 

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
    
    MM=FN{1}*(NN{1}\NN{2})*FN{2}*(NN{2}\NN{1}); 
    ML=FL{1}*(NL{1}\NL{2})*FL{2}*(NL{2}\NL{1}); 

    ne=(log(eig(MM)))/1i/omg0/sum(dd); 
    [~,ind]=sort(real(ne.^2)); 


    ne1Arr(il)=ne(ind(2)); 
    ne2Arr(il)=ne(ind(3)); 

    nl=(log(eig(ML))+0*pi*1i)/1i/omg0/sum(dd); 
    nzlArr(il)=nl(1); 
    
    if (isscalar(lamArr))

        % calculate field distribution within a unit cell using nonlocal
        % TMM
        imd=ind(3); 
        [avec,evD]=eig(MM); 
        a1=avec(:,imd); 
        di=[0 dd(2),0]; % epsPerp, epsZZ, and alpZZ exist
        epsPerpi=[epsPerp,epsPerp(1)];
        epsZZi=[epsZZ,epsZZ(1)];
        alpZZi=[alpZZ,alpZZ(1)]; 

        F1d2=diag(exp(1i*[kzi{1};-kzi{1}]*dd(2))); 
        a2=NN{2}\NN{1}*a1; 
        a3=(NN{1}*F1d2)\NN{2}*FN{2}*a2; 

        apli={[a1(1);a1(2)], [a2(1);a2(2)], [a3(1);a3(2)]}; 
        amini={[a1(3);a1(4)], [a2(3);a2(4)], [a3(3);a3(4)]};
        kzi={kzi{1}(1:2), kzi{2}(1:2), kzi{1}(1:2)}; 

        zplt=(-dd(1)/50:dd(1)/140:dd(1)+dd(2)); 

        [~,Hy]=nonlocalTMMfldZ(di,epsPerpi,kzi,apli,amini,zplt,false); 

        Hmax=max(abs(Hy)); 

        figure(31)
        clf
        plot(zplt,abs(Hy)/Hmax,'-k','LineWidth',2)

        % repeat for local calculations
        imd=2; 
        [avec,evD]=eig(ML); 
        a1=avec(:,imd); 
        F1d2=diag(exp(1i*[kzl{1},-kzl{1}]*dd(2))); 
        a2=NL{2}\NL{1}*a1; 
        a3=(NL{1}*F1d2)\NL{2}*FL{2}*a2; 

        apli={[a1(1)], [a2(1)], [a3(1)]}; 
        amini={[a1(2)], [a2(2)], [a3(2)]};
        kzi={kzl{1}, kzl{2}, kzl{1}}; 

        zplt=(-dd(1):dd(1)/5:dd(1)+dd(2)); 

        [~,Hy]=nonlocalTMMfldZ(di,epsPerpi,kzi,apli,amini,zplt,false); 
        hold on 
        Hmax=max(abs(Hy)); 
        plot(zplt,abs(Hy)/Hmax,'ko','LineWidth',2)

        legend('$|B_y|$', '$|B_y|_{\rm loc}$','Location','southeast', ...
            'interpreter','latex')

        xlabel('$z (\mu m)$','Interpreter','latex')
        ylabel('$B_y (\rm arb. units)$', 'Interpreter','latex')
        box on 
        set(gca,'FontSize',18)
        xlim([-dd(1), dd(1)+dd(2)])

        ylim([0 1.1])
        r=rectangle('Position',[0 0 dd(1) 1.1],'FaceColor',[0.8, 0.9, 0.95],'EdgeColor','none');
        uistack(r,'bottom')
        r=rectangle('Position',[dd(2) 0 sum(dd) 1.1],'FaceColor',0.75*[1 1 1],'EdgeColor','none');
        uistack(r,'bottom')
        xlim([0 sum(dd)])

    end 

end 

if isscalar(lamArr)
    return
end 


%% plotting - nonlocal plot

[~,iln]=min(abs(lamArr(:)-5.05).^2+abs(nxArr(:)-0.5).^2); 

L0=25; 
figure(22)
clf 
hold on 

surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne1Arr.^2)), 'EdgeColor','none','FaceColor','interp')
surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne2Arr.^2)),'EdgeColor','none','FaceColor','interp')

colormap jet
set(gca, 'FontSize', 18)
box on 

set(gca,'CLim', [-1 1]*10)

xlabel('$\lambda_0$', 'Interpreter','latex')
ylabel('${\tilde k}_x^2$', 'Interpreter','latex')
zlabel('${\tilde k}_{z}^2$', 'Interpreter','latex')

zt=[-1e4,-100,-1,0,1,100]*L0;
% create labels: 
ztl=cell(1,length(zt)); 
for iz=1:length(zt)
    if abs(zt(iz))<=L0
        ztl{iz}=num2str(zt(iz)); 
    else 
        pref=log10(abs(zt(iz)/L0));
        if (pref)==1
            pc='10';
        else 
            pc=['10^',num2str(pref)]; 
        end 
        if zt(iz)<0
            ztl{iz}=['-',num2str(L0),' ',pc];
        else 
            ztl{iz}=[num2str(L0),' ',pc];
        end 
    end 
end 

zticks(zScaleFun(L0,zt))
zticklabels(ztl)


view(3)
colorbar

%% local plot
figure(23)
clf 

L0=20; 
surf(lamArr,nxArr.^2,zScaleFun(L0,real(nzlArr.^2)),'EdgeColor','none','FaceColor','interp')
hold on 
plot3(lamArr(iln),nxArr(iln)^2,zScaleFun(L0,real(nzlArr(iln)^2))+2,'ko','LineWidth',2,'MarkerSize',7)

colormap jet
set(gca, 'FontSize', 18)
xlabel('$\lambda_0 (\mu m)$', 'Interpreter','latex')
ylabel('${\tilde k}_x^2$', 'Interpreter','latex')
zlabel('${\tilde k}_{z}^2$', 'Interpreter','latex')


xticks((4:2:10))
box on 
set(gca,'CLim', [-1 1]*10) 
zlim([-39 25]) 
view(3)
colorbar



