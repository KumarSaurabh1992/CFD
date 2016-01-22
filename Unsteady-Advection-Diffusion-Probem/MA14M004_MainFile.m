%% Author: Kumar Saurabh
% Roll No. MA14M004
% Main File: Run this file
close all;


%% Code for Convection term with Central Difference Scheme
clear all;
clc;

% Gauss Seidel Iteration without relaxation
[X1, time, iteration, t1, t2] = MA14M004_CDS(1);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B4:B4');
xlswrite('outputMA14M004_asign3.xlsx',time,'C4:C4');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D4:D4');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E4:E4');

% Gauss Seidel Iteration with Relaxation
[X2, time, iteration, t1, t2] = MA14M004_CDS(2);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B6:B6');
xlswrite('outputMA14M004_asign3.xlsx',time,'C6:C6');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D6:D6');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E6:E6');

% Conjugate Gradient Method
[X3, time, iteration, t1, t2] = MA14M004_CDS(3);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B8:B8');
xlswrite('outputMA14M004_asign3.xlsx',time,'C8:C8');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D8:D8');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E8:E8');

% Bi-Conjugate Gradient Stabilized Method
[X4, time, iteration, t1, t2] = MA14M004_CDS(4);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B10:B10');
xlswrite('outputMA14M004_asign3.xlsx',time,'C10:C10');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D10:D10');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E10:E10');

% Code for exact Solution:
D = MA14M004_exact(time);

% plotting the results:
figure(1);
subplot(2,2,1);
contourf(X1);
colorbar;
title('Gauss Seidel Method without relaxation');
subplot(2,2,2);
contourf(X2);
colorbar;
title('Gauss Seidel Method with relaxation');
subplot(2,2,3);
contourf(X3);
colorbar;
title('Conjugate Gradient Method');
subplot(2,2,4);
contourf(X4);
colorbar;
title('Biconjugate Gradient Stabilized Method');
suptitle(['Convection term with Central Difference scheme at time = ',num2str(time), 's']);

% Comparison of Results:
X = 1:100;
figure(2);
subplot(2,2,1);
plot(X,D(:,50),X,X1(:,50));
legend('Exact Solution','Gauss Seidel without relaxation');
title('Gauss Seidel without relaxation');
subplot(2,2,2);
plot(X,D(:,50),X,X2(:,50));
legend('Exact Solution','Gauss Seidel with relaxation');
title('Gauss Seidel with relaxation');
subplot(2,2,3);
plot(X,D(:,50),X,X3(:,50));
legend('Exact Solution','Conjugate Gradient Method');
title('Conjugate Gradient Method');
subplot(2,2,4);
plot(X,D(:,50),X,X4(:,50));
legend('Exact Solution','Bi-Conjugate Stabilized Method');
title('Bi-Conjugate Stabilized Method');
suptitle('Comparison of the results of Central Difference Scheme with exact solution along the central line');

% %% Code for Convection term with First Order Upwind.
clear all;
clc;

% Gauss Seidel Iteration without relaxation
[X1, time, iteration, t1, t2] = MA14M004_upwind(1);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B16:B16');
xlswrite('outputMA14M004_asign3.xlsx',time,'C16:C16');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D16:D16');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E16:E16');

% Gauss Seidel Iteration with relaxation
[X2, time, iteration, t1, t2] = MA14M004_upwind(2);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B18:B18');
xlswrite('outputMA14M004_asign3.xlsx',time,'C18:C18');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D18:D18');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E18:E18');

% Conjugate Gradient Method
[X3, time, iteration, t1, t2] = MA14M004_upwind(3);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B20:B20');
xlswrite('outputMA14M004_asign3.xlsx',time,'C20:C20');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D20:D20');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E20:E20');

% Bi-Conjugate Gradient Stabilized Method
[X4, time, iteration, t1, t2] = MA14M004_upwind(4);
xlswrite('outputMA14M004_asign3.xlsx',iteration,'B22:B22');
xlswrite('outputMA14M004_asign3.xlsx',time,'C22:C22');
xlswrite('outputMA14M004_asign3.xlsx',t1,'D22:D22');
xlswrite('outputMA14M004_asign3.xlsx',t2,'E22:E22');

% Code for exact Solution:
D = MA14M004_exact(time);

% plotting the results
figure(3);
subplot(2,2,1);
contourf(X1);
colorbar;
title('Gauss Seidel Method without relaxation');
subplot(2,2,2);
contourf(X2);
colorbar;
title('Gauss Seidel Method with relaxation');
subplot(2,2,3);
contourf(X3);
colorbar;
title('Conjugate Gradient Method');
subplot(2,2,4);
contourf(X4);
colorbar;
title('Biconjugate Gradient Stabilized Method');
suptitle(['Convection term with First Order Upwind at time = ',num2str(time),'s']);

%comparison of results
X = 1:100;
figure(4);
subplot(2,2,1);
plot(X,D(:,50),X,X1(:,50));
legend('Exact Solution','Gauss Seidel without relaxation');
title('Gauss Seidel without relaxation');
subplot(2,2,2);
plot(X,D(:,50),X,X2(:,50));
legend('Exact Solution','Gauss Seidel with relaxation');
title('Gauss Seidel with relaxation');
subplot(2,2,3);
plot(X,D(:,50),X,X3(:,50));
legend('Exact Solution','Conjugate Gradient Method');
title('Conjugate Gradient Method');
subplot(2,2,4);
plot(X,D(:,50),X,X4(:,50));
legend('Exact Solution','Bi-Conjugate Stabilized Method');
title('Bi-Conjugate Stabilized Method');
suptitle('Comparison of the results of First Order Upwind with exact solution along the central line');