%********************************************************************
%  Prior, VAR,  Data Specification and Global Setting
%
%*******************************************************************/

subT     = 1;               % choose subsample */ 
mprior   = 2;               % choose prior */
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
  tau     =  5.0;              % 2.0 overall tightness */ 
  d       =   4;               % scaling down the variance for the coefficients of a distant lag */ 
  w       =   1;               % number of observations used for obtaining the prior for the covariance matrix of error terms*/ 
  lambda  =   0.1;             % 2.3 tuning parameter for coefficients for constant*/ 
  mu      =   0.1;             % 12 tuning parameter for the covariance between coefficients*/

else
   tau     =   0.5;            % 0.1 overall tightness */ 
   d       =   4;              % scaling down the variance for the coefficients of a distant lag */ 
   w       =   1;              % number of observations used for obtaining the prior for the covariance matrix of error terms*/ 
   lambda  =   2.3;            % 2.3 tuning parameter for coefficients for constant*/ 
   mu      =   12;             % 1 tuning parameter for the covariance between coefficients*/
   
end
