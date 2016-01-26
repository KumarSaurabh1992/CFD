%% Individual Details
% Author: Kumar Saurabh
% Roll No. MA14M004
% Main File

close all;
clear;
clc;
%% Reading data from Input file %%
Biot = xlsread('inputMA14M004_assign2.xlsx','B2:B6');
partition = xlsread('inputMA14M004_assign2.xlsx','A2:A3');

%% Computation 
V = 1;
Max = zeros(2,5);
Max1 = zeros(2,5);

for k = 1:5
    
    Bi = Biot(k);
    for p = 1:2
        N = partition(p);
        deltax = 0.5/(N - 1);
        deltay = 0.5/(N - 1);
        
        B = zeros(N*N,N*N);
        C = zeros(N*N,1);
        
        n = 1;
        %% Computation of coeffecient matrix %%
        for j = 1:N
            for i = 1:N
                if (i == 1)
                    
                    if(j == 1)
                        B(n,((j - 1)*N + i+1)) = 2;
                        B(n,(j - 1)*N + i) = -4;
                        B(n,((j*N) + i)) = 2;
                    elseif (j == N)
                        B(n,((j-1)*N + i)) = -(4 + 2*Bi*deltax);
                        B(n,((j - 2)*N + i)) = 2;
                        B(n,((j - 1)*N+i+1)) = 2;
                    else
                        B(n,((j - 1) * N + i)) = -4;
                        B(n,((j - 2) * N + i)) = 1;
                        B(n,((j)*N + i)) = 1;
                        B(n,((j - 1) * N + i + 1)) = 2;
                    end
                elseif(i == N)
                    if(j == 1)
                        B(n,((j-1)*N + i)) = -(4 + 2*Bi*deltax);
                        B(n,((j - 1)*N + i - 1)) = 2;
                        B(n,((j)*N + i)) = 2;
                        
                    elseif(j == N)
                        B(n,((j - 1)*N + i)) = - (4 + (4*Bi*deltax));
                        B(n,((j - 1)*N + i - 1)) = 2;
                        B(n,((j - 2)*N + i)) = 2;
                    else
                        B(n,((j - 1)*N + i)) = -(4 + (2*Bi*deltax));
                        B(n,((j - 1)*N + i - 1)) = 2;
                        B(n,((j - 2)*N + i)) = 1;
                        B(n,((j)*N + i)) = 1;
                    end
                elseif(j == 1)
                    B(n,((j - 1)*N + i)) = -4;
                    B(n,((j - 1)*N + i - 1)) = 1;
                    B(n,((j - 1)*N + i + 1)) = 1;
                    B(n,((j)*N + i)) = 2;
                elseif(j == N)
                    B(n,((j - 1)*N + i)) = -(4 + 2*Bi*deltay);
                    B(n,((j - 1)*N + i - 1)) = 1;
                    B(n,((j - 1)*N + i + 1)) = 1;
                    B(n,((j - 2)*N + i)) = 2;
                else
                    B(n,((j - 1)*N + i)) = -4;
                    B(n,((j - 1)*N + i - 1)) = 1;
                    B(n,((j - 1)*N + i + 1)) = 1;
                    B(n,((j)*N + i)) = 1;
                    B(n,((j - 2)*N + i)) = 1;
                end
                C(n) = -(deltax)^2;
                n = n + 1;
            end
        end
        
        %% Solution of the system of Equation       
        X = MA14M004_Gauss(B,C);  % Solution by Gauss Siedel Method
        X1 = B\C; % Solution by Direct Method
        
        %% Updating the temperature to corresponding temperature of the slab
        D = zeros(N,N);
        D1 = zeros(N,N);
        l = 1;
        for j = 1:N
            for i = 1:N
                D(i,j) = X(l);
                D1(i,j) = X1(l);
                l = l  + 1;
            end
        end
        
        %% Calculation of maximum temperature
        Max(p,k) = max(max(D)); % from Gauss Siedel Method
        Max1(p,k) = max(max(D1)); %from Direct Method
        
        %% Plotting of figure
        figure(V);
        subplot(2,2,p*2 - 1);
        contour(D,50);
        colorbar;
        title(['Contour plot for ',num2str(N), ' X ',num2str(N),'grid']);
        subplot(2,2,p*2);
        length = 0:deltax:0.5;
        plot(length,D(:,1));
        axis tight;
        xlabel('x');
        ylabel(['{\theta} (\circC)']);
        title(['Tempearture distribution along x = 0 for ',num2str(N), ' X ',num2str(N),'grid']);
        
        
        
    end
    suptitle(['Biot Number =  ',num2str(Bi)]);
    V = V+1;
end

%% Writing the final result to Excel file
xlswrite('outputMA14M004_asign2.xlsx',Max,'C2:G3');
xlswrite('outputMA14M004_asign2.xlsx',Max1,'C6:G7');