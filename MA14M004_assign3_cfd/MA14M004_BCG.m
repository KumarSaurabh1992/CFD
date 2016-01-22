%% Author: Kumar Saurabh
% Roll No. MA14M004
% Bi-Conjugate Gradient Stabilized Method to solve AX = B

%% Function:
function [ X,n ] = MA14M004_BCG( A, B, X,n )
len = size(B,1);
%X = ones(len,1);
X1 = ones(len,1);
Q = MA14M004_mul(A,X);
r = (B - Q);
r_star = r;
alpha = 1;
omega = 1;
rho = 1;
v = 0;
p1 = 0;
error = 199;

while (error > 10^(-6))
    rho_1 = r_star'*r;
    beeta = (rho_1/rho)*(alpha/omega);
    p1 = r + beeta*(p1 - omega*v);
    v = MA14M004_mul(A,p1);
    alpha = rho_1/(r_star'*v);
    s = r - alpha*v;
    t = MA14M004_mul(A,s);
    omega = (t'*s)/(t'*t);
    X = X + omega*s + alpha*p1;
    r = s - omega*t;
    rho = rho_1;
    n = n + 1;
    
    error = 0;
    for i = 1:len
        error = error + abs((X1(i) - X(i))/X1(i));
    end
    X1 = X;
    %disp(error);
end

end
