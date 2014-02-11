%***************************************************************/
%*                                                             */
%*     VAR analysis: compute marginal data density             */
%*                                                             */
%***************************************************************/


vm_spec 


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


