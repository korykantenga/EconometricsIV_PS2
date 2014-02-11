function [y] = irf_Uhlig(Phi,Sigma,H,j)

% This function draws the IRFs from the posterior density. The matrix Omega
% is drawn from the unit ball.

[m,n] = size(Sigma); %#ok<*ASGLU>
y = ones(n,H);

r = zeros(1,j);
e = randn(n,1);

Omega = (e./norm(e));

for i=1:H
    y(:,i) = impulseresponce(Phi,Sigma,Omega,i);
end

end