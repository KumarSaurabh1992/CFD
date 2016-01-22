%% Individual Details
% Author: Kumar Saurabh
% Roll No. MA14M004
% Gauss siedel Method to solve system of linear equation

%% Function :
function [  X, n ] = MA14M004_Gauss( A, B, X, n, relax)

len = size(B,1);
error = 1999;
X1 = zeros(len,1);
%X = zeros(len,1);
N = sqrt(len);

%% Traditional Gauss Siedel %%
% while(error > 10^(-6))
%     for j = 1:len
%         sum = 0;
%         for k = 1:len
%             if j ~= k
%                 sum = sum + A(j,k)*X(k);
%             end
%         end
%         X(j) =  X1(j) + relax*((B(j) - sum)/A(j,j) - X1(j));
%     end
%     alpha = zeros(len,1);
%     for i = 1:len
%             alpha(i) = abs((X(i) - X1(i))/X(i));
%     end
%     error = max(alpha);
%
%     disp(error);
%     clear X1;
%     X1 = X;
%     n = n + 1;
%
% end

%% Fast Gauss Siedel specific to this case %%
while(error > 10^(-6))
    for k = 1:len
        row = floor((k-1)/N) + 1;
        column = mod((k - 1),N) + 1;
        index = (row - 1)*N + column;
        left = index - 1;
        right = index + 1;
        top = index + N;
        bottom = index - N;
        
        if (row == 1)
            if (column == 1)
                sum = A(k,5)*X(top) + A(k,4)*X(right);
            elseif (column == N)
                sum = A(k,5)*X(top) + A(k,2)*X(left);
            else
                sum = A(k,5)*X(top) + A(k,2)*X(left) + A(k,4)*X(right);
            end
        elseif (row == N)
            if (column == 1)
                sum = A(k,1)*X(bottom) + A(k,4)*X(right);
            elseif (column == N)
                sum = A(k,1)*X(bottom) + A(k,2)*X(left);
            else
                sum = A(k,1)*X(bottom) + A(k,2)* X(left) + A(k,4)*X(right);
            end
        else
            if(column == 1)
                sum = A(k,1)*X(bottom) + A(k,4)*X(right) + A(k,5)*X(top);
            elseif (column == N)
                sum = A(k,1)*X(bottom) + A(k,2)*X(left) + A(k,5)*X(top);
            else
                sum = A(k,1)*X(bottom) + A(k,2)*X(left) + A(k,5)*X(top) + A(k,4)*X(right);
            end
        end
        X(k) =  X1(k) + relax*((B(k) - sum)/A(k,3) - X1(k));
    end
    
    % Calculation of error
    error = 0;
    for i = 1:len
        error = error + abs((X(i) - X1(i))/X(i));
    end
    
    
    %disp(error);
    clear X1;
    X1 = X;
    n = n + 1;
    
end



%disp(n);
end

