clear;
clc;
close all;

H = 1; % m
L = 14; % m
nx = 281;  % number of x nodes
ny = 21;  % number of y nodes
dx = L/(nx-1);
dy = H/(ny-1);
t = 0;
delta_t = 0.001;
Re = 100;
% nondimensionalization
dX = dx/H;
dY = dy/H;
PSI = zeros(ny,nx);
vort = zeros(ny,nx);
Pr = 0.7;
%% Boundary conditions
PSI(1,:) = 0; % bottom
PSI(ny,:) = 1; % top
w = zeros(ny,nx);
u_prime = zeros(ny,nx);
u = zeros(ny,nx);
v = zeros(ny,nx);
thetha= zeros(ny,nx);
relax = 1.8;
for j = 1:ny
    PSI(j,1) = (j-1)*dY; % inlet
end
diff = 20;
while (diff > 10^(-6))
    W = PSI;
    u_prime = u;
    err = 10;
    n = 1;
    while (err > 10^(-6))
        tempPSI = PSI;
        for i = 2:nx - 1
            for j = 2:ny - 1
                PSI(j,i) = relax*(1/(2/dX^2+2/dY^2)*(((PSI(j,i+1)+PSI(j,i-1))/dX^2+...
                    (PSI(j+1,i)+PSI(j-1,i))/dY^2) + vort(j,i))) + (1 - relax)*PSI(j,i);
            end
        end
        PSI(2:ny,nx) = 2*PSI(2:ny,(nx - 1)) - PSI(2:ny,(nx - 2)); 
        err = 0;
        for i = 2:nx - 1
            for j = 2:ny - 1
                err = err+ (abs(PSI(j,i) - tempPSI(j,i))/PSI(j,i));
            end
        end
        
        n = n+1;
    end
    
    vort(ny,2:nx) = 2.0*(PSI(ny,2:nx) - PSI(ny - 1,2:nx))/(dy^2); % top wall
    vort(1,2:nx) = -2.0*(PSI(2,2:nx))/(dy^2); % bottom wall
    vort(2:ny-1,nx) = vort(2:ny -1,nx - 1);
    for i = 2:nx - 1
        for j = 2:ny - 1
        w(j,i) = vort(j,i) + delta_t*(-((PSI(j+1,i) - PSI(j - 1,i))/(2*dy))*((vort(j,i+1) - vort(j,i - 1))/(2*dx))+...
            +((PSI(j,i+1) - PSI(j,i - 1))/(2*dx))*((vort(j+1,i) - vort(j - 1,i))/(2*dy))+...
            (1/Re)*(((vort(j,i+1)-2*vort(j,i)+vort(j,i-1))/(dx^2)) + (vort(j+1,i)-2*vort(j,i)+vort(j - 1,i))/dy^2));

        end
    end
    vort = w;
    t = t + delta_t;
    
% contourf(PSI);
% pause(0.5)


psi = PSI*1*1; % dimentionalization

for j = 2:ny-1
    u(j,:) = (psi(j+1,:)-psi(j-1,:))/(2*dy);
end
for i = 2:nx-1
    v(:,i) = -(psi(:,i+1)-psi(:,i-1))/(2*dx);
end
diff = max(max(abs((W(:,:) - PSI(:,:)))));
u(:,1) = 1;                 % inlet
v(:,1) = 0;                 % inlet
u(1,:) = 0;                 % bottom wall
u(ny,:) = 0;                % top wall
u(2:ny,nx) = u(2:ny,nx - 1); % outflow



delta_t1 = 1/(2*(((max(max(abs(u))))/dx) + (max(max(abs(v))))));
delta_t2 = (1/(2*((1/dx^2) + (1/dy^2))))/(1.3*10^(-5));
delta_t3 = min(delta_t2,delta_t1);
delta_t = min(delta_t3,dx);
%disp(diff);
err = 20;
%while (err < 10^(-6))
thetha_prime = thetha;
thetha(1,:) = 1;
thetha(ny,:) = 1;
thetha(2:ny - 1,1) = 0;
thetha(2:ny - 1,nx) = thetha(2:ny - 1,nx -1);
for i = 2:nx - 1
    for j = 2: ny - 1
        Conv1 = -u(j,i)*(thetha_prime(j,i+1) - thetha_prime(j,i - 1))/(2*dx);
        Conv2 = -v(j,i)*(thetha_prime(j+1,i) - thetha_prime(j-1,i))/(2*dy);
        diff1 = (1/(Re*Pr))*(thetha_prime(j,i-1) + thetha_prime(j,i+1) - 2*thetha_prime(j,i))/dx^2;
        diff2 = (1/(Re*Pr))*(thetha_prime(j+1,i) + thetha_prime(j-1,i) - 2*thetha_prime(j,i))/dy^2;
        thetha(j,i) = thetha_prime(j,i) + delta_t*(Conv1 + Conv2 + diff1 + diff2);
    end
end


end
Sum = 0;
K = zeros(nx,1);
for i = 1:nx
        K(i) =  (thetha(1,i) - thetha(2,i))/dy;
end
Nu = sum(K)/nx;

%% Printing the result
fprintf('Steady state reached \n');
fprintf('Local Nusselt Number near the walls = %.3f\n', Nu);
%% Zooming at inlet:
U = zeros(21,50);
V = zeros(21,50);

for i = 1:50
    for j = 1:21
        U(j,i) = u(j,i);
        V(j,i) = v(j,i);
    end
end

temp = zeros(21,80);


for i = 1:80
    for j = 1:21
        temp(j,i) = thetha(j,i);
    end
end
%% Plotting the result:

% velocity contours
figure(1);
subplot(2,1,1);
contourf(u);
colorbar;
title('Velocity contours across channel');
subplot(2,1,2);
contourf(U);
colorbar;
title('Velocity contours zoomed at inlet');
suptitle('Velocity Contours');

% velocity profile
y = 0:dy:1;
figure(2);
plot(y,u(:,250));
xlabel('y');
ylabel('Velocity (in non-dimensionalized form)');
suptitle('Velocity profile');

% temperature contours
figure(3);
subplot(2,1,1);
contourf(thetha);
colorbar;
title('Tempearture contours across channel');
subplot(2,1,2);
contourf(temp);
colorbar;
title('Temperature contours zoomed at inlet');
suptitle('Temperature Contours');

% velocity profile
y = 0:dy:1;
figure(4);
plot(y,thetha(:,250));
xlabel('y');
ylabel('Temperature (in non-dimensionalized form)');
suptitle('Temperature profile');

% Nusselt number
x = 0:dx:14;
figure(5);
plot(x,K);
xlabel('Length along the channel');
ylabel('Nusselt number');
suptitle('Variation of Nusselt Number along the length of channel');
