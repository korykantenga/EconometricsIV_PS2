%Data set:
%Output, Inflation, Interest Rates, InvVelocity  

load datampshock.txt
YY=datampshock;
ti=linspace(1964,0.2, size(YY,1))';
nobs=size(YY,1);

%Convert velocity into real money balances

YY(:,4)=YY(:,4)+YY(:,1);
