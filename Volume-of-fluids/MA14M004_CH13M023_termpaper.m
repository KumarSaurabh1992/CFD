clear all;
clc;

ny = 51;
nx = 51;
L = 1;
H = 1;
dx = L/(nx - 1);
dy = H/(ny - 1);
t = 0;
u = zeros(ny,nx);
v = zeros(ny,nx);
F0 = zeros(ny,nx);
F1 = zeros(ny,nx);
H_V = zeros(ny,nx);
yi = zeros(ny,nx);
xi = zeros(ny,nx);
Dy_Dx = zeros(ny,nx);
Dx_Dy = zeros(ny,nx);

epsilon = 10^(-6);

for j = 1:ny
    for i = 1:nx
        v(j,i) = 0.5; %-0.5*sin(pi*((i - 1)*dx-0.5))*cos(pi*((j - 1)*dy) - 0.5);
        u(j,i) = 0.5; %0.5*cos(pi*((i - 1)*dx-0.5))*sin(pi*((j - 1)*dy) - 0.5);
    end
end


for j = 1:ny
    for i = 1:nx
        if((j > 2) && (j < 10) && (i > 2) && (i < 10))
            F0(j,i) = 1;
        end
    end
end
time = 1.5;
dt = 0.001;
while (t < time)
    F1 = F0;
    %% Deciding the orientation
    for j = 2:ny - 1
        for i = 2:nx - 1
            yi(j,i) = F0(j - 1,i)*dy + F0(j,i)*dy + F0(j + 1,i)*dy;
            xi(j,i) = F0(j - 1,i)*dx + F0(j,i)*dx + F0(j ,i+1)*dx;
        end
    end
    
    for j = 2:ny - 1
        for i = 2:nx - 1
            Dy_Dx = 2*(yi(j,i + 1) - yi(j,i - 1))/(4*dy);
            Dx_Dy = 2*(xi(j+1,i) - xi(j - 1,i))/(4*dx);
            if (abs(Dy_Dx) < abs(Dx_Dy))      %horizontal
                H_V(j,i) = 1;       
            else
                H_V(j,i) = 0;
            end
        end
    end
    %% Deciding the coordinates of acceptor and donor cell
    for j = 1:ny - 1
        for i = 1: nx - 1
            Vx = abs(u(j,i)*dt);
            if(u(j,i) > 0)
                donor_x = i;
                acceptor_x = i + 1;
            else
                donor_x = i + 1;
                acceptor_x = i;
            end
            if (H_V(j,donor_x) == 1) % convection is more parallel so donor cell is used
                Acc_Don_x = donor_x;
            else % convection is more normal so acceptor cell is used
                Acc_Don_x = acceptor_x;
            end
            if(F0(j,acceptor_x) < epsilon)
                Acc_Don_x = acceptor_x;
            end
            CF = max((1 - F0(j,Acc_Don_x))*Vx - (1 - F0(j,donor_x))*dx,0);
            df = min(F0(j,Acc_Don_x)*Vx + CF, F0(j,donor_x)*dx);
            F1(j,donor_x) = F1(j,donor_x) - df/dx;
            F1(j,acceptor_x) = F1(j,acceptor_x) + df/dx;
            %% y velocity
            Vy = abs(v(j,i)*dt);
            if(v(j,i) > 0)
                donor_y = j;
                acceptor_y = j + 1;
            else
                donor_y = j + 1;
                acceptor_y = j;
            end
            if (H_V(donor_y,i) == 0) % convection is more normal so acceptor cell is used
                Acc_Don_y = donor_y;
            else % convection is more parallel so donor cell is used
                Acc_Don_y = acceptor_y;
            end
            if(F0(acceptor_y,i) < epsilon)
                Acc_Don_y = acceptor_y;
            end
            CF = max((1 - F0(Acc_Don_y,i))*Vy - (1 - F0(donor_y,i))*dy,0);
            df = min(F0(Acc_Don_y,i)*Vy + CF, F0(donor_y,i)*dy);
            F1(donor_y,i) = F1(donor_y,i) - df/dy;
            F1(acceptor_y,i) = F1(acceptor_y,i) + df/dy;
        end
    end
    %% Bookkeeping adjustment
    for i = 2:nx-1
        for j = 2:ny-1
            if (F1(j,i) > 1 - epsilon)
                F1(j,i) = 1;
            end
            if(F1(j,i) < epsilon)
                F1(j,i) = 0;
                if (F1(j,i+1) > 1.0 - epsilon)
                    F1(j,i+1) = F1(j,i+1) - 1.1*epsilon;
                end
                if (F1(j,i-1) > 1.0 - epsilon)
                    F1(j,i-1) = F1(j,i-1) - 1.1*epsilon;
                end
                if (F1(j+1,i) > 1.0 - epsilon)
                    F1(j+1,i) = F1(j+1,i) - 1.1*epsilon;
                end
                if (F1(j-1,i) > 1.0 - epsilon)
                    F1(j-1,i) = F1(j-1,i) - 1.1*epsilon;
                end
            end
            
             F0(j,i) = min(1,max(F1(j,i),0));
        end
           
    end
    t = t + dt;
    contourf(F0);
    suptitle(['time = ',num2str(t),'s']);
    pause(0.01);

end