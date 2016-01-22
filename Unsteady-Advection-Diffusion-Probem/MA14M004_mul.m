%% Author: Kumar Saurabh
% Roll No. MA14M004
% Matrix Multiplication:

%% Function:
function [ C ] = MA14M004_mul( A,B )
len = size(B,1);
C = zeros(len,1);
N = sqrt(len);
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
            C(k,1) = A(k,5)*B(top) + A(k,4)*B(right) + A(k,3)*B(index);
        elseif (column == N)
            C(k,1) = A(k,5)*B(top) + A(k,2)*B(left) + A(k,3)*B(index);
        else
            C(k,1) = A(k,5)*B(top) + A(k,2)*B(left) + A(k,4)*B(right) + A(k,3)*B(index);
        end
    elseif (row == N)
        if (column == 1)
            C(k,1) = A(k,1)*B(bottom) + A(k,4)*B(right) + A(k,3)*B(index);
        elseif (column == N)
            C(k,1) = A(k,1)*B(bottom) + A(k,2)*B(left)+ A(k,3)*B(index);
        else
            C(k,1) = A(k,1)*B(bottom) + A(k,2)* B(left) + A(k,4)*B(right) + A(k,3)*B(index);
        end
    else
        if(column == 1)
            C(k,1) = A(k,1)*B(bottom) + A(k,4)*B(right) + A(k,5)*B(top) + A(k,3)*B(index);
        elseif (column == N)
            C(k,1) = A(k,1)*B(bottom) + A(k,2)*B(left) + A(k,5)*B(top) + A(k,3)*B(index);
        else
            C(k,1) = A(k,1)*B(bottom) + A(k,2)*B(left) + A(k,5)*B(top) + A(k,4)*B(right) + A(k,3)*B(index);
        end
        
    end
    
end

end

