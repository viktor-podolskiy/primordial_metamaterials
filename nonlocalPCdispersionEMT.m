clear 

% dispersion of the modes in periodic nonlocal bi-layered structures;
% reproduces Fig.5(c,d) of "Primordial Media: the shrouded realm of composite materials"

lamArr=(0.5:0.025:3); % decrease size to improve image quality
nxArr=(0.001:0.025:3);7.5; 

[lamArr,nxArr]=meshgrid(lamArr,nxArr); 

omgArr=2*pi./lamArr; 

eps1Arr=(epsDrudePl(lamArr,1,1));
eps2Arr=2+0*eps1Arr; 


alp1Arr=-3/5*1e-5./(1./lamArr.*(1./lamArr+0.2i)); 
alp1=1e-4*(-1+0.1i); %+0.00001i; 
alp2=0.5e-5*(1+0.2i); %+0.00001i; 


dd=1e-4*[1,1]; 


% arrays to store dispersion
ne1Arr=0*nxArr; 
ne2Arr=0*nxArr; 

nz1EMT=0*nxArr; 
nz2EMT=0*nxArr; 

%% running

for il=1: numel(lamArr)

    nx=nxArr(il); 
    omg0=omgArr(il);

    epsZZ=[eps1Arr(il), eps2Arr(il)];
    epsPerp=epsZZ; 
    alp1=alp1Arr(il);
    alpZZ=[alp1,alp2];



    %% EMT-based ini setup: 
    
    ezzEMT=sum(dd.*epsZZ)/sum(dd); % VP fix
    eperpEMT=sum(dd.*epsPerp)/sum(dd);
    alpEMT=sum(dd)/sum(dd./alpZZ); 


    DD=sqrt((alpEMT*eperpEMT+ezzEMT)^2-4*eperpEMT*alpEMT*nx^2); 
    nzz0=sqrt((alpEMT*eperpEMT-ezzEMT+[DD;-DD])/2/alpEMT);
    nz1EMT(il)=nzz0(1); 
    nz2EMT(il)=nzz0(2); 

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
    
    end 
    
    MM=FN{1}*(NN{1}\NN{2})*FN{2}*(NN{2}\NN{1}); 

    ne=(log(eig(MM))+0*pi*1i)/1i/omg0/sum(dd); 
    [~,ind]=sort(real(ne.^2)); 


    ne1Arr(il)=ne(ind(2)); 
    ne2Arr(il)=ne(ind(3)); 
   
end 



%% plotting dispersion data

L0=20; 
figure(22)
clf 
hold on 

surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne1Arr.^2)), 'EdgeColor','none','FaceColor','interp')
surf(lamArr,nxArr.^2,zScaleFun(L0,real(ne2Arr.^2)),'EdgeColor','none','FaceColor','interp')

colormap jet
set(gca, 'FontSize', 18)
box on 

xlim([0.5 3])
set(gca,'CLim', [-1 1]*3)

zlim([-120 100])

xlabel('\lambda/\lambda_p') 
ylabel('${\tilde k}_x^2$', 'Interpreter','latex')
zlabel('${\tilde k}_{z}^2$', 'Interpreter','latex')

zt=[-10000*L0,-100*L0,-L0,0,L0,100*L0];
setZTicks
view(3)
colorbar

%% nonlocal EMT plot
figure(25)
clf 
hold on 

surf(lamArr,nxArr.^2,zScaleFun(L0,real(nz1EMT.^2)), 'EdgeColor','none','FaceColor','interp')
surf(lamArr,nxArr.^2,zScaleFun(L0,real(nz2EMT.^2)),'EdgeColor','none','FaceColor','interp')

colormap jet
set(gca, 'FontSize', 18)
box on 

xlim([0.5 3])
set(gca,'CLim', [-1 1]*3)

xlabel('\lambda/\lambda_p') 
ylabel('${\tilde k}_x^2$', 'Interpreter','latex')
zlabel('${\tilde k}_{z}^2$', 'Interpreter','latex')

zlim([-120 100])

zt=[-10000*L0,-100*L0,-L0,0,L0,100*L0];

setZTicks

view(3)
colorbar



