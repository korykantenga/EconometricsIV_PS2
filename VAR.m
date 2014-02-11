%==========================================================================
%                       VAR WITH MINNESOTA PRIOR                      
%
%                      See README_VAR for details
%
%
% Author: Luigi Bocola     lbocola@sas.upenn.edu
% Date  : 06/20/2010
%
% Modified: Kory Kantenga
% Date  : 02/10/2014
%
%==========================================================================


%=========================================================================
%                             HOUSEKEEPING
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
%         GENERATE DUMMY OBSERVATIONS FROM MINNESOTA PRIOR 
%=========================================================================
 
vm_dummy

%=========================================================================
%     DEFINITION OF DATA, LAG STRUCTURE AND POSTERIOR SIMULATION
%=========================================================================

[Tdummy,n] = size(YYdum);
[Tobs,n]   = size(YYact);
X          = [XXact; XXdum];
Y          = [YYact; YYdum];
n          = n;                 % Number of variables in the VAR
p          = 4;                 % Number of lags in the VAR
T          = Tobs+Tdummy;
nsim       = 10000;             % Number of draws from Posterior Density
nburn      = 0.2*nsim;          % Number of draws to discart
F          = zeros(n*p,n*p);    % Matrix for Companion Form
I          = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end


%=========================================================================
%               OLS ESTIMATOR FOR PHI AND SSR (SIGMA)
%=========================================================================

Phi_tilde = inv(X'*X)*X'*Y;
Sigma     = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);

% Matrices for collecting draws from Posterior Density

Sigmap    = zeros(nsim,n,n);
Phip      = zeros(nsim,n*p+1,n);
largeeig  = zeros(nsim,1);
counter   = 0;

%=========================================================================
%            DRAWS FROM POSTERIOR DENSITY (DIRECT SAMPLING)
%=========================================================================
disp('                                                                  ');
disp('        BAYESIAN ESTIMATION OF VAR: DIRECT SAMPLING...            ');
disp('                                                                  ');

for j=1:nsim

    
    % Draws from the density Sigma | Y
    
    sigma   = iwishrnd(Sigma,T-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv(X'*X)));
    
    % Rearrange vec(Phi) into Phi
    
    Phi     = reshape(phi_new,n*p+1,n);
        
    Sigmap(j,:,:) = sigma;
    Phip(j,:,:)   = Phi;
    
    Phi = Phi(1:n*p,:);

    
    % Calculating largest eigenvalue of Companion form
    
     F(1:n,1:n*p)    = Phi';

     eigen           = eig(F);
     eigen           = max(eigen);
     largeeig(j)     = abs(eigen);
     counter         = counter +1; 
     
     if counter==2000
disp(['         DRAW NUMBER:   ', num2str(j)]);
disp('                                                                  ');
disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
disp('                                                                  ');

     counter = 0;
     
     end
     
end

%=========================================================================
%                        MARGINAL DATA DENSITY
%=========================================================================

vm_mdd

mdd = lnpYY;               % Marginal Data Density


%=========================================================================
%           FIGURE 1: LARGEST EIGENVALUE (Companion Form)
%=========================================================================

% 
% pnames = strvcat( 'Largest Eigenvalue (Recursive Average)',...
%     'Largest Eigenvalue (Posterior Marginal Density)');
% 
% figure('Position',[20,20,900,600],'Name',...
%     'Largest Eigenvalue (Companion Form)','Color','w')
% 
% rmean = zeros(nsim,1);
% 
% for i=1:nsim
%     rmean(i) = mean(largeeig(1:i));
% end
% 
% subplot(1,2,1), plot(rmean,'LineStyle','-','Color','b',...
%         'LineWidth',2.5), hold on
% title(pnames(1,:),'FontSize',13,'FontWeight','bold');
% 
% [density,x]  = ksdensity(largeeig(nburn:end));
% 
% subplot(1,2,2), plot(x,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5), hold on
% 
% title(pnames(2,:),'FontSize',13,'FontWeight','bold');
% 
% path=l;
% 
% disp(['         ELAPSED TIME:   ', num2str(toc)]);


elapsedtime=toc;


