close all;
clear all;
clc;

%% Details 
% Author Name: Kumar Saurabh
% Roll No. MA14M004

%% Reading the data from input file
Time_interval = xlsread('inputMA14M004_assign1.xlsx','Sheet1','A2:A4');
alpha = xlsread('inputMA14M004_assign1.xlsx','Sheet1','B2');
len = xlsread('inputMA14M004_assign1.xlsx','Sheet1','C2');
delta_x = xlsread('inputMA14M004_assign1.xlsx','Sheet1','D2');
for i = 1:3
    
    %% Initializing the input
    delta_t = Time_interval(i);
    X = 0:delta_x:len;
    time = 0;
    figure(i);
    
    %% Calculations 
    x_size = (len/delta_x) + 1;
    gamma = alpha*delta_t/(delta_x^2);
    T_old = zeros(x_size,1);
    T_old(1) = 100; %Boundary Condition
    
    %% At time t = 0 s
    T_0 = T_old;
    subplot(2,3,1);
    plot(X,T_0);
    xlabel('X (m)');
    ylabel('Temperature (\circC)');
    title('time = 0 s');
    %% Calculations for various times
    while(time < 10)
        T_new = zeros(x_size,1);
        T_new(1) = 100;
        time = time + delta_t;
        
        for j = 2:(x_size - 1)
            T_new(j) = gamma*T_old(j+1) + (1 - (2*gamma))*T_old(j) + gamma*T_old(j - 1);
        end
        
        %% At time t = 0.5 s
        if(abs(time - 0.5)<(10^(-6)))
            T_point5 = T_new;
            subplot(2,3,2);
            plot(X,T_point5);
            title('time = 0.5 s');
            xlabel('X (m)');
            ylabel('Temperature (\circC)');
        end
        
        %% At time t = 1 s
        if(abs(time - 1) < (10^(-6)))
            T_1 = T_new;
            subplot(2,3,3);
            plot(X,T_1);
            title('time = 1 s');
            xlabel('X (m)');
            ylabel('Temperature (\circC)');
        end
        
        %% At time t = 2 s
        if(abs(time - 2.000) < (10^(-6)))
            T_2 = T_new;
            subplot(2,3,4);
            plot(X,T_2);
            title('time = 2 s');
            xlabel('X (m)');
            ylabel('Temperature (\circC)');
        end
        
        %% At time t = 5 s
        if(abs(time - 5.000) < (10^(-6)))
            T_5 = T_new;
            subplot(2,3,5);
            plot(X,T_5);
            title('time = 5 s');
            xlabel('X (m)');
            ylabel('Temperature (\circC)');
        end
        
        %% At time t = 10 s
        if(abs(time - 10.000) < (10^(-6)))
            T_10 = T_new;
            subplot(2,3,6);
            plot(X,T_10);
            title('time = 10 s');
            xlabel('X (m)');
            ylabel('Temperature (\circC)');
        end
        
        %% Updation of Temperature
        T_old = T_new;
        
        
    end
    
    suptitle(['Temperature distribution with {\Delta t} =  ',num2str(delta_t)])
end