%% Author: Kumar Saurabh
% Roll No. MA14M004
% Conjugate Gradient Method to solve AX = B

%% Function:
function [ X,n ] = MA14M004_CG( A, B, X, n )

len = size(B,1);
%X = ones(len,1);
X1 = ones(len,1);
Q = MA14M004_mul(A,X);
r = (B - Q);
d = r;
rnew = r'*r;

error = 199;

while (error > 10^(-6))
    q = MA14M004_mul(A,d);
    alpha = rnew/(d'*q);
    X = X + alpha*d;
    r = r - alpha*q;
    rold = rnew;
    rnew = r'*r;
    beeta = rnew/rold;
    d = r + beeta*d;
    error = 0;
    for i = 1:len
        error = error + abs((X1(i) - X(i))/X1(i));
    end
    X1 = X;
    n = n + 1;
end
end