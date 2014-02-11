function y = impulseresponce(Phi0,Sigma,Omega,H)

[m , n]= size(Phi0);
p = (m-1)/n;

Phi = Phi0(1:n*p,:);
F = zeros(n*p,n*p);
F(1:n,1:n*p) = Phi';
I = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end

M = zeros(n,n*p);
M(1:n,1:n) = I;
A = zeros(n*p,n*p);
A(1:n,1:n) = chol(Sigma)';
Omega1 = [Omega;zeros((n-1)*p,1)];
y = M*(F^(H-1))*A*Omega1;



