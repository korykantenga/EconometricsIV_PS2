function [y] = irf(Phi,Sigma,H,j)

% This function draw the IRFs from the posterior density. The matrix Omega
% is drawn from the unit ball until the sign restrictions on the impulse
% response of Inflation(y(2,:)), the Interest Rate (y(3,:)) and Real Money
% (y(4,:)) are satisfied for j periods.

[m,n] = size(Sigma); %#ok<*ASGLU>
y = ones(n,H);
d=0;
r = zeros(1,j);

while d==0
    e = randn(n,1);
    
    Omega = (e./norm(e));
    
    for i=1:H
        y(:,i) = impulseresponce(Phi,Sigma,Omega,i);
    end
    
    if y(2,1:j)<=r & y(3,1:j)>=r & y(4,1:j)<=r;
        d=1;
    else
        d=0;
    end
end

end
