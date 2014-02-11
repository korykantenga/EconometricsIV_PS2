%********************************************************************
%  Generate Table 1 for Hyperparameter Selection
%
%*******************************************************************/

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
%                             PRIOR SELECTION LOOP
%=========================================================================


mdd = zeros(4,5);

for nprior=1:5
    mprior   = nprior;          % choose prior */
    subT     = 1;               % choose subsample */
    %dataselstr = "d1";         % choose data set */
    nlags_   = 4;               % number of lags   */
    nirf_    = 40;              % choose number of responses */
    T0       = 4;               % size of pre-sample */
    nex_     = 1;               % number of exogenous vars; 1 means intercept only */
    
    %Setting up the values for hyperparameters for the Minnesota prior **/
    
    if (mprior == 1);
        tau     =   0.01;            % 0.01 overall tightness */
        d       =   4;               % scaling down the variance for the coefficients of a distant lag */
        w       =   1;               % number of observations used for obtaining the prior for the covariance matrix of error terms*/
        lambda  =   1;               % 2.3 tuning parameter for coefficients for constant*/
        mu      =   1;               % 12 tuning parameter for the covariance between coefficients*/
    elseif (mprior == 2);
        tau     =   0.1;             % 0.1 overall tightness */
        d       =   4;               % scaling down the variance for the coefficients of a distant lag */
        w       =   1;               % number of observations used for obtaining the prior for the covariance matrix of error terms*/
        lambda  =   1;               % 2.3 tuning parameter for coefficients for constant*/
        mu      =   1;               % 12 tuning parameter for the covariance between coefficients*/
    elseif (mprior == 3);
        tau     =   0.5;            % 0.5 overall tightness */
        d       =   4;              % scaling down the variance for the coefficients of a distant lag */
        w       =   1;              % number of observations used for obtaining the prior for the covariance matrix of error terms*/
        lambda  =   1;              % 2.3 tuning parameter for coefficients for constant*/
        mu      =   1;              % 12 tuning parameter for the covariance between coefficients*/
    elseif (mprior == 4);
        tau     =   1.0;             % 1.0 overall tightness */
        d       =   4;               % scaling down the variance for the coefficients of a distant lag */
        w       =   1;               % number of observations used for obtaining the prior for the covariance matrix of error terms*/
        lambda  =   1;               % 2.3 tuning parameter for coefficients for constant*/
        mu      =   1;               % 12 tuning parameter for the covariance between coefficients*/
    elseif (mprior == 5);
        tau     =  2.0;              % 2.0 overall tightness */
        d       =   4;               % scaling down the variance for the coefficients of a distant lag */
        w       =   1;               % number of observations used for obtaining the prior for the covariance matrix of error terms*/
        lambda  =   1;               % 2.3 tuning parameter for coefficients for constant*/
        mu      =   1;               % 12 tuning parameter for the covariance between coefficients*/
    end
    
    
    %********************************************************
    % Import data series                                    *
    %*******************************************************/
    vm_loaddata
    
    nv      = size(YY,2);     %* number of variables */
    nobs    = size(YY,1)-T0;  %* number of observations */
    
    %********************************************************
    % Dummy Observations                                    *
    %*******************************************************/
    
    %** Obtain mean and standard deviation from expandend pre-sample data
    
    YY0     =   YY(1:T0+16,:);
    ybar    =   mean(YY0)';
    sbar    =   std(YY0)';
    premom  =   [ybar sbar];
    
    % Generate matrices with dummy observations
    hyp = [tau; d; w; lambda; mu];
    [YYdum, XXdum, breakss] = varprior_h(nv,nlags_,nex_,hyp,premom);
    
    % Actual observations
    
    YYact = YY(T0+1:T0+nobs,:);
    ti = ti(T0+1:T0+nobs,:);
    XXact = zeros(nobs,nv*nlags_);
    i = 1;
    
    while (i <= nlags_)
        XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
        i = i+1;
    end
    
    % last column of XXact = constant
    XXact = [XXact ones(nobs,1)];
    
    % store observations
    save data YYdum XXdum YYact XXact
    
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
    fprintf('        THE MINNESOTA PRIOR TIGHTNESS IS %1.2f',tau);
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
    
    
    %********************************************************
    % Import data series to determine length of sample      *
    %*******************************************************/
    vm_loaddata
    
    nobs = size(YY,1)-T0;
    nv   = size(YY,2);
    
    
    %********************************************************
    %* open output generated by vm_dummy.g
    %* dummy  YYdum XXdum
    %* actual YYact XXact
    %*******************************************************/
    
    load data
    
    YY=[YYdum' YYact']';
    XX=[XXdum' XXact']';
    YYdum=YYdum;
    XXdum=XXdum;
    
    n_total = size(YY,1);
    n_dummy = n_total-nobs;
    nv     = size(YY,2);
    k      = size(XX,2);
    
    
    %********************************************************
    %* Compute the log marginal data density for the VAR model
    %*******************************************************/
    
    Phi0   = (inv(XXdum'*XXdum))*(XXdum'*YYdum);
    S0     = (YYdum'*YYdum)-YYdum'*XXdum*(inv(XXdum'*XXdum))*XXdum'*YYdum;
    
    Phi1   = (inv(XX'*XX))*(XX'*YY);
    S1     = (YY'*YY)-YY'*XX*(inv(XX'*XX))*XX'*YY;
    
    %* compute constants for integrals
    
    i=1;
    gam0=0;
    gam1=0;
    
    while i <= nv;
        gam0=gam0+log(gamma(0.5*(n_dummy-k+1-i)));
        gam1=gam1+log(gamma(0.5*(n_total-k+1-i)));
        i=i+1;
    end;
    
    %** dummy observation
    
    lnpY0 = -nv*(n_dummy-k)*0.5*log(pi)-(nv/2)*log(abs(det(XXdum'*XXdum)))-(n_dummy-k)*0.5*log(abs(det(S0)))+nv*(nv-1)*0.25*log(pi)+gam0;
    
    %** dummy and actual observations
    
    lnpY1 = -nv*(n_total-k)*0.5*log(pi)-(nv/2)*log(abs(det(XX'*XX)))-(n_total-k)*0.5*log(abs(det(S1)))+nv*(nv-1)*0.25*log(pi)+gam1;
    lnpYY = lnpY1-lnpY0;
    tau;
    save lnpYY
    
    mdd(1,nprior) = tau;                 % Prior for Marginal Data Density
    mdd(2,nprior) = 0.2;                 % Prior Probability
    mdd(3,nprior) = lnpYY;               % Marginal Data Density
    mdd(4,nprior) = 0;                   % Posterior Probability
    
end

[temp1,temp2] = max(mdd(3,:));

mdd(4,temp2) = 1;

clc
disp('Latex Code for Table 1');
for j=1:size(mdd,1)
    fprintf('%8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',...
        mdd(j,1), mdd(j,2), mdd(j,3), mdd(j,4), mdd(j,5))
end