
%=========================================================================
%            BAYESIAN POINT FORECASTS & CREDIBLE BANDS
%=========================================================================

% Length of Forecast
stepH    = 12;
forecast = zeros(stepH,4,nsim);

for j=1:nsim
    
    %Initialisation
    simYY    = Y(1:p,:);
    
    %Extract Draw of Phi
    Phi             = zeros(size(Phip,2),size(Phip,3));
    Phi(:,:)        = Phip(j,:,:);
    F(1:n,1:n*p)    = Phi(1:n*p,:)';
    
    % Update X and Y matrices
    simXX           = simYY(1:p,:);
    simYY           = [0 0 0 0; simYY(1:p-1,:)];
    PhiC            = zeros(size(simYY,1),size(simYY,2));
    PhiC(1,:)       = Phi(end,:);
    simXX           = reshape(simXX',n*p,1);
    simYY           = reshape(simYY',n*p,1); %#ok<NASGU>
    PhiC            = reshape(PhiC',n*p,1);
    
    % Forecast for t=1
    simYY           = F*simXX;
    simYY(1:4)      = PhiC(1:4,:)+simYY(1:4);
    forecast(1,:,j)   = simYY(1:4)';
    
    %Update Forecast at time t>1
    t=2;
    while (t<stepH+1)
        simXX           = simYY;
        simYY           = [0; 0; 0; 0; simYY(1:end-4)];
        simYY           = F*simXX;
        simYY(1:4)      = PhiC(1:4,:)+simYY(1:4);
        forecast(t,:,j)   = simYY(1:4)';
        t=t+1;
    end
    
end

meanForecast = mean(forecast,3);

%Credible Sets
CSets1 = zeros(stepH,4);
for s = 1:stepH
    CSets1(s,:) = [s; meanForecast(s,1); hpdi(forecast(s,1,:),90)];
end

CSets2 = zeros(stepH,4);
for s = 1:stepH
    CSets2(s,:) = [s; meanForecast(s,2); hpdi(forecast(s,2,:),90)];
end

CSets3 = zeros(stepH,4);
for s = 1:stepH
    CSets3(s,:) = [s; meanForecast(s,3); hpdi(forecast(s,3,:),90)];
end

CSets4 = zeros(stepH,4);
for s = 1:stepH
    CSets4(s,:) = [s; meanForecast(s,4); hpdi(forecast(s,4,:),90)];
end


%Plots of Mean Forecasts and 90% Error Bands

figure;
subplot(2,2,1)
hold on
plot(CSets1(:,1),CSets1(:,2),'-r')
plot(CSets1(:,1),CSets1(:,3),':')
plot(CSets1(:,1),CSets1(:,4),':')
xlabel('Quarters Ahead Forecast')
ylabel('Residual Log(GDP per Capita)')
title('ln(GDP per Capita)')

subplot(2,2,2)
hold on
plot(CSets2(:,1),CSets2(:,2),'-r')
plot(CSets2(:,1),CSets2(:,3),':')
plot(CSets2(:,1),CSets2(:,4),':')
xlabel('Quarters Ahead Forecast')
ylabel('\Delta GDP Deflator')
title('Inflation')

subplot(2,2,3)
hold on
plot(CSets3(:,1),CSets3(:,2),'-r')
plot(CSets3(:,1),CSets3(:,3),':')
plot(CSets3(:,1),CSets3(:,4),':')
xlabel('Quarters Ahead Forecast')
ylabel('Federal Funds Rate')
title('Interest Rate')

subplot(2,2,4)
hold on
plot(CSets4(:,1),CSets4(:,2),'-r')
plot(CSets4(:,1),CSets4(:,3),':')
plot(CSets4(:,1),CSets4(:,4),':')
xlabel('Quarters Ahead Forecast')
title('Real Money Balances')

print -depsc2 VARForecast_Q1_6.eps
