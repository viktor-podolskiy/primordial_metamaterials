clear 
% dispersion of the modes in periodic nonlocal bi-layered structures;
% reproduces Figs.2(a) and 3(b) of "Primordial Media: the shrouded realm of composite materials"


lamArr=(0.5:0.01:3);

eps1Arr=(epsDrudePl(lamArr,1,1));
eps2Arr=2+0*eps1Arr; 


alp1Arr=-3/5*1e-5./(1./lamArr.*(1./lamArr+0.2i)); 
alp1=1e-4*(-1+0.1i); %+0.00001i; 
alp2=0.5e-5*(1+0.2i); %+0.00001i; 


% plot permittivity data
figure(1)
clf
yyaxis left
hold on 
plot(lamArr,real(eps1Arr),'LineWidth',2)
plot(lamArr,real(eps2Arr),'LineWidth',2)
ylabel('\epsilon')
ylim([-10 5])
yyaxis right
hold on 
plot(lamArr,real(alp1Arr),'LineWidth',2)
plot(lamArr,real(0*alp1Arr+alp2),'LineWidth',2)
ylabel('\alpha')

box on 
xlabel('\lambda/\lambda_p')
set(gca,'FontSize',18)
legend('Re(\epsilon_p)','Re(\epsilon_d)','Re(\alpha_p)','Re(\alpha_d)', 'Location','southwest')

%% computations

TArr=0*lamArr; 


figure(10)
clf 

ang0=60; 
ddArr=[0.05 0.05 0.003 0.001 0.0001 0.0001]; %array of layer thicknesses
totThk=0.2; % total thickness of the composite


TddArr=zeros(length(ddArr),length(lamArr));

% TtotAngArr=RAngArr; 

legsArr=cell(1,length(ddArr)); % legends

for idd=1:length(ddArr)

    dd=ddArr(idd)*[1 1]; 
    
    numBlrs=totThk/2/ddArr(idd); % number of bi-layers
    % set legends
    legsArr{idd}=[num2str(ddArr(idd)),'\lambda_{pl}']; 
    if idd==1 
        legsArr{idd}=[legsArr{idd},',local'];
    elseif idd==length(ddArr)
        legsArr{idd}=[legsArr{idd},',EMT'];
    end 
    

    for ilam=1:length(lamArr)
        lam0=lamArr(ilam); 
        omg0=2*pi/lam0; 
        kx=omg0*sind(ang0); 
    
        %set material parameters
        di=[]; epsPerp=[]; epsZZ=[]; alpZZ=[];  
    
        for ib=1:numBlrs
            di=[di,dd]; 
    
            epsPerp=[epsPerp, eps1Arr(ilam), eps2Arr(ilam)];
            epsZZ=[epsZZ,eps1Arr(ilam), eps2Arr(ilam)]; 
            alpZZ=[alpZZ,alp1Arr(ilam),alp2]; 
        end 
        
        if idd==1
            alpZZ=alpZZ*0; % convert to local stack
        elseif idd==length(ddArr)
            % replace the composite with the effective-medium layer
            di=totThk; 
            epsPerp=(dd(1)*eps1Arr(ilam)+dd(2)*eps2Arr(ilam))/sum(dd); 
            epsZZ=(dd(1)*eps1Arr(ilam)+dd(2)*eps2Arr(ilam))/sum(dd);
            alpZZ=sum(dd)/(dd(1)/alp1Arr(ilam)+dd(2)/alp2);
        end 
    

        % "wrap" the composite in input/output layers
        di=[0, di, 0]; 
        epsPerp=[1 epsPerp 1]; 
        epsZZ=[1 epsZZ 1]; 
        alpZZ=[0 alpZZ 0]; 
    
    
        [Ttot,kzi] = nonlocalTMMcalcT(omg0,kx,di,epsPerp,epsZZ,alpZZ); % calculate transmission
    
    
        % RArr(ilam)=abs(RMati{1})^2; 
        TArr(ilam)=abs(Ttot)^2*(epsPerp(end)*kzi{1}(1)/epsPerp(1)/kzi{end}(1)); 
    
    end 

    % RAngArr(idd,:)=RArr; 
    TddArr(idd,:)=TArr; 
    
    %% plot
    figure(10)
    hold on 
    if idd==1 
        plot(lamArr([5:10:end]), TArr(([5:10:end])), 'o', 'LineWidth',2)
        set(gca,'ColorOrderIndex',1)
    elseif idd==length(ddArr)
        set(gca,'ColorOrderIndex',length(ddArr)-2)
        plot(lamArr([5:10:end]), TArr(([5:10:end])), '^', 'LineWidth',2)
    else 
        plot(lamArr, TArr, '-', 'LineWidth',2)
    end 
    set(gca,'FontSize',18)
    box on 
                
end 

%% add legends
figure(10)
xlabel('\lambda/\lambda_{pl}')
ylabel('Transmission')
legend(legsArr,'Location','northeastoutside')

