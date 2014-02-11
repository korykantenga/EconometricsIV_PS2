%=========================================================================
%      IMPULSE RESPONSE FUNCTIONS USING UHLIG PROCEDURE
%=========================================================================
clc

H      = 24;
restr  = 2;             %Integer indicating number of horizons to restrict
PPhi   = zeros(n*p+1,n);
SSigma = zeros(n,n);
IRFunc = zeros(n,H,nsim);
IRFunc2= zeros(n,H,nsim);
IRFuncV= zeros(nsim,1);
meanIRF= zeros(n,H); %#ok<NASGU>
lowerB = zeros(n,H);
upperB = zeros(n,H);
storeB = zeros(2,n,H);

for j=1:nsim
    SSigma(:,:)   = Sigmap(j,:,:);
    PPhi(:,:)     = Phip(j,:,:);
    IRFunc(:,:,j) = irf_Uhlig(PPhi,SSigma,H,restr);
end


for j = 1:size(IRFunc,3)
    if IRFunc(2,1:restr,j)>0 | IRFunc(3,1:restr,j)<0 | ...
            IRFunc(4,1:restr,j)>0
        IRFuncV(j)=1;
    end
end

for j = 1:size(IRFunc,3)
    if IRFuncV(j)==0
        IRFunc2(:,:,j)=IRFunc(:,:,j);
    else
        IRFunc2(:,:,j)=NaN;
    end
end

tempIRF = IRFunc2(isnan(IRFunc2)==0);
IRFunc2 = reshape(tempIRF,n,H,length(tempIRF)/(n*H));

meanIRF = mean(IRFunc2,3);

for t=1:n
    for s = 1:H
        storeB(:,t,s) = hpdi(IRFunc2(t,s,:),90)';
        lowerB(t,s)   = storeB(1,t,s);
        upperB(t,s)   = storeB(2,t,s);
    end
end

% Plots

figure;
subplot(2,2,1)
hold on
plot(1:H,meanIRF(1,:),'-m')
plot(1:H,lowerB(1,:),':')
plot(1:H,upperB(1,:),':')
title('ln(GDP per Capita)')
hline = refline(0,0);
xlim([0 H])
set(hline,'Color','k','LineStyle','--')

subplot(2,2,2)
hold on
plot(1:H,meanIRF(2,:),'-m')
plot(1:H,lowerB(2,:),':')
plot(1:H,upperB(2,:),':')
title('Inflation')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

subplot(2,2,3)
hold on
plot(1:H,meanIRF(3,:),'-m')
plot(1:H,lowerB(3,:),':')
plot(1:H,upperB(3,:),':')
title('Interest Rate')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

subplot(2,2,4)
hold on
plot(1:H,meanIRF(4,:),'-m')
plot(1:H,lowerB(4,:),':')
plot(1:H,upperB(4,:),':')
title('Money Balances')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

print -depsc2 VAR_IRF_Q2_5.eps
