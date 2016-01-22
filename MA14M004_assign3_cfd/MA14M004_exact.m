%% Author: Kumar Saurabh
% Roll No. MA14M004
% Computes the exact solution of the given problem.

%% Function
function [D1] = MA14M004_exact(t1)

    N = 100;
    D1 = zeros(N,N);
    delta_x = 1/N;
    delta_y = 1/N;

    for i = 1:N
        for j = 1:N
            D1(i,j) = exp(-2*pi*pi*t1)*cos(pi*((delta_x/2) + (i - 1)*delta_x - 0.5))*cos(pi*((delta_y/2) + (j - 1)*delta_y - 0.5));

        end
    end
    D1 = flipdim(D1,1);
end
