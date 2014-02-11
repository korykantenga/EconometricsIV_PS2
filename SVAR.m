%==========================================================================
%                   IRFs TO A MONETARY POLICY SHOCK                      
%                                                                 
%                             
%                     See README_SVAR for details
%
%
% Author: Luigi Bocola    lbocola@sas.upenn.edu
% Date  : 06/20/2010
%==========================================================================


%=========================================================================
%                              Housekeeping
%=========================================================================

tic
close all
clear all
clc


l = path;
path('Impulse Responses',path);
path('Minnesota Prior',path);
path('Data',path);


%=========================================================================
%             Load Data and Dummy Observations (Minnesota Prior) 
%=========================================================================
 
load data

%=========================================================================
%     Definition of Data, Lag Structure of VAR and Posterior Simulation
%=========================================================================

[Tdummy,n] = size(YYdum);
[Tobs,n]   = size(YYact);
X          = [XXact; XXdum];
Y          = [YYact; YYdum];
n          = n;     % Number of variables in the VAR
p          = 4;     % Number of lags in the VAR
T          = Tobs+Tdummy;
H          = 40;    % Number of Periods for the computation of IRFs
nsim       = 1000;  % Number of draws from Posterior Density


%=========================================================================
%                OLS Estimator For Phi and SSR (Sigma)
%=========================================================================

Phi_tilde = inv(X'*X)*X'*Y;
Sigma     = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);

% Matrices for collecting draws from Posterior Density
Sigmap = zeros(nsim,n,n);
Phip   = zeros(nsim,n*p+1,n);
IRF    = zeros(nsim,n,H+1);
counter= 0;
%=========================================================================
%       Draws from the Posterior Density using Direct Sampling
%=========================================================================

disp('                                                                  ');
disp('       SIGN RESTRICTIONS FOR IRFs : ACCEPTANCE SAMPLING...        ');
disp('                                                                  ');

for j=1:nsim
    
    % Draws from the density Sigma | Y
    
    sigma   = iwishrnd(Sigma,T-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv(X'*X)));
    
    % Rearrange vec(Phi) into Phi
    
    Phi     = reshape(phi_new,n*p+1,n);
    
    % Computation of Impulse Responce using sign restrictions for the first
    % two periods; 
    
   IR            = irf(Phi,sigma,H+1,2);
    
    % Store the draws from the posterior density
   Sigmap(j,:,:) = sigma;
   Phip(j,:,:)   = Phi;
   IRF(j,:,:)    = IR;
   counter       = counter+1;
   
   if counter==300
disp(['         DRAW NUMBER:   ', num2str(j)]);
disp('                                                                  ');
disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
disp('                                                                  ');
      counter=0;
      
   end
end

%=========================================================================
% Compute mean and 90% pointwise credible sets (HPD) for Impulse Responces
%=========================================================================

[irmean,ir05,ir95] = impulsemoment(IRF);


%=========================================================================
%           Figure 1: Impulse Responce Functions (Posterior means and 90%
%                     Credible Sets)
%=========================================================================

pnames = strvcat( 'Output','Inflation', 'Interest Rate','Real Money');

figure('Position',[20,20,900,600],'Name',...
    'IRFs to a Contractionary Monetary Policy Shock','Color','w')


for i=1:4
subplot(2,2,i), plot(irmean(i,:),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
xlim([1 40])    
plot(ir05(i,:),'LineStyle',...
    ':','Color',...
    'r','LineWidth',2.5), hold on
plot([1 40],[0 0], 'LineStyle','--','Color','black',...
        'LineWidth',2.5)
%if i==1 || i==3    
%h2 = legend('Posterior Mean','90% Credible Set ',1);
%else
%h2 = legend('Posterior Mean','90% Credible Set',4);
%end    
%set(h2,'Interpreter','none');
plot(ir95(i,:),'LineStyle',...
    ':','Color',...b
    'r','LineWidth',2.5), hold off
title(pnames(i,:),'FontSize',13,'FontWeight','bold');
end


path = l;

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;
