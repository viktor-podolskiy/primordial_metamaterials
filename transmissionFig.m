%-----------------------------
% Intrinsically nonlocal metamaterials manuscript 
% Supplementary Information
% The script calculates the figures related to the period-study of the
% nonlocal composites. To select the particular composite, choose 
% appropriate values for dD and dM.
%
% copyright (2025) U Mass Lowell 
% Jacob LaMountain and Viktor Podolskiy 
%-----------------------------

clear 
lamArr=(4:0.02:10); 

numBlrs=6; %includes "cap" layer
lampD=17.2; % original material
lampM=6.3; 

eps0D=12.3; 

% --- select the stack configuration below
dB=0.005; % barrier thickness
dD=0.03; % thickness of nInAs layer
dM=0.04; % thickness of n++ layer

% dB=0.005; 
% dD=0.07; 
% dM=0.08; 

% dB=0.005; 
% dD=0.11; 
% dM=0.12; 


RArr=0*lamArr; 
TArr=0*lamArr; 

RMArr=0*lamArr; 
TMArr=0*lamArr; 
TTArr=0*lamArr; 



epsZZmArr=epsDrude(lamArr,lampM,eps0D); 
epsPerpmArr=epsZZmArr; 
epsD=epsDrude(lamArr,lampD,eps0D); 
epsB=10; 

% data from primordial manuscript
alpZZmArr=-eps0D*4.1e-5*3/5/lampM^2./(1./lamArr.*(1./lamArr+0.015i)); 
alpZZdArr=-eps0D*2.6e-5*3/5/lampD^2./(1./lamArr.*(1./lamArr+0.015i)); 
alpZZbArr=1e-5*(2.75+1.5i)+0*epsZZmArr; 


figure(10)
clf 

angArr=[1 15 30 45 60];

RAngArr=zeros(length(angArr),length(lamArr)); 
TAngArr=RAngArr; 
RMAngArr=zeros(length(angArr),length(lamArr)); 
TMAngArr=RAngArr; 

TtotAngArr=RAngArr; 


for iAng=1:length(angArr)
    
    ang0=angArr(iAng); 


    for ilam=1:length(lamArr)
        lam0=lamArr(ilam); 
        omg0=2*pi/lam0; 
        kx=omg0*sind(ang0); 
    
        di=[]; epsPerp=[]; epsZZ=[]; alpZZ=[];  
    
        for ib=1:numBlrs
            if ib==1
                di=[di,0.01,dB, dD, dB];
                epsPerp=[epsPerp, eps0D, epsB, epsD(ilam), epsB];
                epsZZ=[epsZZ,eps0D,epsB, epsD(ilam), epsB];
                alpZZ=[alpZZ,0,alpZZbArr(ilam), alpZZdArr(ilam), alpZZbArr(ilam)];
            else
                di=[di,dM, dB, dD, dB];
                epsPerp=[epsPerp, epsPerpmArr(ilam), epsB, epsD(ilam), epsB];
                epsZZ=[epsZZ,epsZZmArr(ilam),epsB, epsD(ilam), epsB];
                alpZZ=[alpZZ,alpZZmArr(ilam),alpZZbArr(ilam), alpZZdArr(ilam), alpZZbArr(ilam)];
            end
        end 
        
    
        % forward propagation
        di=[0, di, 0]; 
        epsPerp=[1 epsPerp 10]; 
        epsZZ=[1 epsZZ 10]; 
        alpZZ=[0 alpZZ 0]; 
    
    
        [Ttot,kzi,RMati,TMati] = nonlocalTMMcalcT(omg0,kx,di,epsPerp,epsZZ,alpZZ); 
    
    
        RArr(ilam)=abs(RMati{1})^2; 
        TArr(ilam)=abs(Ttot)^2*(epsPerp(end)*kzi{1}(1)/epsPerp(1)/kzi{end}(1)); 
    
        % backward propagation
        di=fliplr(di); 
        epsPerp=fliplr(epsPerp); 
        epsZZ=fliplr(epsZZ); 
        alpZZ=fliplr(alpZZ); 
        [TtotM,kziM,RMatiM] = nonlocalTMMcalcT(omg0,kx,di,epsPerp,epsZZ,alpZZ); 
    
        RMArr(ilam)=abs(RMatiM{1})^2; 
        TMArr(ilam)=abs(TtotM)^2*(epsPerp(end)*kziM{1}(1)/epsPerp(1)/kziM{end}(1)); 
    
    
        % collecting propagation to incorporate the response of the
        % optically-thick substrate
        rho=abs((kziM{1}(1)*epsPerp(end)-kziM{end}(1)*epsPerp(1))/...
            (kziM{1}(1)*epsPerp(end)+kziM{end}(1)*epsPerp(1)))^2; 
        TTArr(ilam)=TArr(ilam)*(1-rho)/(1+rho*RMArr(ilam)); 
        
    end 
    
    RAngArr(iAng,:)=RArr; 
    TAngArr(iAng,:)=TArr; 
    RMAngArr(iAng,:)=RMArr; 
    TMAngArr(iAng,:)=TMArr; 
    TtotAngArr(iAng,:)=TTArr; 
    figure(10)
    hold on 
    plot(lamArr, TTArr, '-', 'LineWidth',2)
    set(gca,'FontSize',18)
    box on 
            
end 

% finish the figure setup
xlabel('wavelength, um')
ylabel('Transmission')
ylim([0 1])
legend('0^o', '15^o', '30^o', '45^o', '60^o')