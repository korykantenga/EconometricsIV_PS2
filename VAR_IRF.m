%=========================================================================
%      IMPULSE RESPONSE FUNCTIONS WITH PURE SIGN RESTRICTIONS
%=========================================================================

H      = 24;
PPhi   = zeros(n*p+1,n);
SSigma = zeros(n,n);
IRFunc = zeros(n,H,nsim);
meanIRF= zeros(n,H); %#ok<NASGU>
lowerB = zeros(n,H);
upperB = zeros(n,H);
storeB = zeros(2,n,H);



for j=1:nsim
    SSigma(:,:)   = Sigmap(j,:,:);
    PPhi(:,:)     = Phip(j,:,:);
    IRFunc(:,:,j) = irf(PPhi,SSigma,H,2);
end

meanIRF = mean(IRFunc,3);

for t=1:n
    for s = 1:H
        storeB(:,t,s) = hpdi(IRFunc(t,s,:),90)';
        lowerB(t,s)   = storeB(1,t,s);
        upperB(t,s)   = storeB(2,t,s);
    end
end

% Plots

figure;
subplot(2,2,1)
hold on
plot(1:H,meanIRF(1,:),'-r')
plot(1:H,lowerB(1,:),':')
plot(1:H,upperB(1,:),':')
title('ln(GDP per Capita)')
hline = refline(0,0);
xlim([0 H])
set(hline,'Color','k','LineStyle','--')

subplot(2,2,2)
hold on
plot(1:H,meanIRF(2,:),'-r')
plot(1:H,lowerB(2,:),':')
plot(1:H,upperB(2,:),':')
title('Inflation')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

subplot(2,2,3)
hold on
plot(1:H,meanIRF(3,:),'-r')
plot(1:H,lowerB(3,:),':')
plot(1:H,upperB(3,:),':')
title('Interest Rate')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

subplot(2,2,4)
hold on
plot(1:H,meanIRF(4,:),'-r')
plot(1:H,lowerB(4,:),':')
plot(1:H,upperB(4,:),':')
title('Money Balances')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

print -depsc2 VAR_IRF_Q2_4.eps
