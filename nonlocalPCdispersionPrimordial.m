clear 

% % use single wavelength/single nx for field distribution plot
% lamArr=5.05; 
% nxArr=0.5; 

% use multi-wavelength/nx sweep to reproduce dispersion plot
lamArr=(4:0.025:10); % wavelength array, microns
nxArr=(0.001:0.025:3); % x-components of dimensionalized wavenumbers

% re-mesh to create 2D arrays for plotting
[lamArr,nxArr]=meshgrid(lamArr,nxArr); 
omgArr=2*pi./lamArr; % omega/c array

% local and nonlocal parts of permittivity of materials
eps1Arr=(epsDrude(lamArr,6.3,12.3));
eps2Arr=10+0*eps1Arr; 

alp1Arr=-12.3*4.1e-5*3/5/6.3^2./(1./lamArr.*(1./lamArr+0.015i)); 
alp2=1e-4*(1+0.3i); 

% layer thicknesses, in microns
dd=[0.01,0.01]; 

% setting up output arrays 
ne1Arr=0*nxArr; % nonlocal PC, nz, mode 1
ne2Arr=0*nxArr; % nonlocal PC, nz, mode 2 

nzlArr=0*nxArr; % local PC, nz


%% running

for il=1: numel(lamArr)

    %set parameters for current wavelength/nx combination

    nx=nxArr(il); 
    omg0=omgArr(il); 

    epsZZ=[eps1Arr(il), eps2Arr(il)];
    epsPerp=epsZZ; 
    alp1=alp1Arr(il);
    alpZZ=[alp1,alp2];



   
    % nonlocal material-specific matrices
    NN=cell(2,1); 
    FN=cell(2,1); 
    kzi=cell(2,1); 

    % local material-specific matrices 
    NL=cell(2,1); 
    FL=cell(2,1); 
    kzl=cell(1,1); 
    
    for ilr=[1,2] % iterate over layers in bi-layered composite
        % nonlocal calculations
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
    
        % local calculations
        nzl=sqrt(epsPerp(ilr)*(1-nx^2/epsZZ(ilr)));
        if(imag(nzl)<0)
            nzl=-nzl; 
        end 
        NL{ilr}=[1, 1; ...
            epsPerp(ilr)/nzl, -epsPerp(ilr)/nzl]; 
        FL{ilr}=diag(exp(1i*omg0*[nzl, -nzl]*dd(ilr))); 
        kzl{ilr}=omg0*nzl; 
    end 
    
    % nonlocal and local matrices of one period of the structure
    MM=FN{1}*(NN{1}\NN{2})*FN{2}*(NN{2}\NN{1}); 
    ML=FL{1}*(NL{1}\NL{2})*FL{2}*(NL{2}\NL{1}); 

    % effective indices of modes in nonlocal composite
    ne=(log(eig(MM))+0*pi*1i)/1i/omg0/sum(dd); 
    [~,ind]=sort(real(ne.^2)); 
    ne1Arr(il)=ne(ind(2)); 
    ne2Arr(il)=ne(ind(3)); 

    % effective index of mode in local composite
    nl=(log(eig(ML))+0*pi*1i)/1i/omg0/sum(dd); 
    nzlArr(il)=nl(1); 
    
    if (isscalar(lamArr)) % plotting field distributions across the unit cell
        imd=ind(3); 
        % imd=ind(2); 
        [avec,~]=eig(MM); 
        a1=avec(:,imd); 

        % extending material stack to one full period
        di=[0 dd(2),0]; % epsPerp, epsZZ, and alpZZ exist
        epsPerpi=[epsPerp,epsPerp(1)]; 
        epsZZi=[epsZZ,epsZZ(1)];
        alpZZi=[alpZZ,alpZZ(1)]; 

        F1d2=diag(exp(1i*[kzi{1};-kzi{1}]*dd(2))); 
        a2=NN{2}\NN{1}*a1; % amplitudes in material/layer 2
        a3=(NN{1}*F1d2)\NN{2}*FN{2}*a2; % amplitudes in material 1, one period away

        apli={[a1(1);a1(2)], [a2(1);a2(2)], [a3(1);a3(2)]}; 
        amini={[a1(3);a1(4)], [a2(3);a2(4)], [a3(3);a3(4)]};
        kzi={kzi{1}(1:2), kzi{2}(1:2), kzi{1}(1:2)}; 

        % positions to plot field at (must fall within period and a half)
        zplt=(-dd(1)/50:dd(1)/140:dd(1)+dd(2)); 

        [~,Hy]=nonlocalTMMfldZ(di,epsPerpi,epsZZi,alpZZi,kzi,apli,amini,zplt,false); 

        Hmax=max(abs(Hy)); 

        figure(31)
        clf
        plot(zplt,real(Hy)/Hmax, zplt,imag(Hy)/Hmax, zplt,abs(Hy)/Hmax,'-k','LineWidth',2)
    


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

        [~,Hy]=nonlocalTMMfldZ(di,epsPerpi,epsZZi,alpZZi,kzi,apli,amini,zplt,false); 
        hold on 
        Hmax=max(abs(Hy)); 
        plot(zplt,abs(Hy)/Hmax,'ko','LineWidth',2)

        legend('${\rm Re}(B_y)$','${\rm Im}(B_y)$', '$|B_y|$', '$|B_y|_{\rm loc}$','Location','southeast', ...
            'interpreter','latex')

        xlabel('$z (\mu m)$','Interpreter','latex')
        ylabel('$B_y (\rm arb. units)$', 'Interpreter','latex')
        box on 
        set(gca,'FontSize',18)
        xlim([-dd(1), dd(1)+dd(2)])

        ylim([-1.5 1])
        plot([dd(1) dd(2)], [-1.5 1], '--k','LineWidth',1.5, 'HandleVisibility','off')
        xlim([0 sum(dd)])
    
    end 


end 

if isscalar(lamArr)
    return
end 


%% plotting dispersion

L0=20; 

% nonlocal PC
figure(22)
clf 
hold on 

surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne1Arr.^2)), 'EdgeColor','none')
surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne2Arr.^2)),'EdgeColor','none')

colormap jet
set(gca, 'FontSize', 18)
box on 

set(gca,'CLim', [-1 1]*10)
zlim([-125 100])

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

xticks((4:2:10))

zticks(zScaleFun(L0,zt))
zticklabels(ztl)
view(3)
colorbar

% plot PC
figure(23)
clf 

L0=20; 
surf(lamArr,nxArr.^2,zScaleFun(L0,real(nzlArr.^2)),'EdgeColor','none')

colormap jet
set(gca, 'FontSize', 18)
xlabel('$\lambda_0 (\mu m)$', 'Interpreter','latex')
ylabel('${\tilde k}_x^2$', 'Interpreter','latex')
zlabel('${\tilde k}_{z}^2$', 'Interpreter','latex')


xticks((4:2:10))
box on 
set(gca,'CLim', [-1 1]*10) 
zlim([-35 25])

view(3)
colorbar


