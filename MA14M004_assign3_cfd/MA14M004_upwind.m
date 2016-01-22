%% Author: Kumar Saurabh
% Roll No. MA14M004
% This file computes the results when convection term is discretized by
% First Order Upwind.

%% Function:
function [D, t, it, t1, t2] = MA14M004_upwind(option)
tic;
t2 = cputime;

thetha = 0.5; % Crank Nicholson Scheme
N = 100;
gama = 1;
len = 1;
delta_x = len/N;
delta_y = len/N;
delta_t = 0.01;
ap_not = delta_x*delta_y/delta_t;
B = zeros(N*N,5);
X = zeros(N*N,1);
coeff = zeros(N*N,5);
n = 1;
De = 1*delta_y/delta_x;
Dw = 1*delta_y/delta_x;
Dn = 1*delta_x/delta_y;
Ds = 1*delta_x/delta_y;
u = @(x,y) gama*cos(pi*(x - 0.5))*sin(pi*(y - 0.5));
v = @(x,y) -gama*sin(pi*(x - 0.5))*cos(pi*(y - 0.5));
phi = @(x,y) cos(pi*(x - 0.5))*cos(pi*(y - 0.5));

%% Computation of Coefficient Matrix
for j = 1:N
    for i = 1:N
        east = (i)*delta_x;
        East = delta_x/2 + (i)*delta_x;
        west = (i - 1)*delta_x;
        West = delta_x/2  + (i - 2)*delta_x;
        north = j*delta_y;
        North = delta_y/2 + (j)*delta_y;
        south = (j - 1)*delta_y;
        South = delta_y/2 + (j - 2)*delta_y;
        Point_x = delta_x/2 + (i - 1)*delta_x;
        Point_y = delta_y/2 + (j - 1)*delta_y;
        %% left face
        Fe = (u(east,Point_y))*delta_y;
        Fw = (u(west,Point_y))*delta_y;
        Fn = (v(Point_x,north))*delta_x;
        Fs = (v(Point_x,south))*delta_x;
        
        if (i == 1)
            if (j == 1) %left bottom corner
                aw = 0;
                an = Dn + max(-Fn,0);
                as = 0;
                ae = De + max(-Fe,0);
                Sp = -(2*Dw + 2*Ds);
                Fw = 0;
                Fs = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                B(n,4) = -thetha *ae;
                B(n,3) = ap;
                B(n,5) = -thetha * an;
            elseif (j == N) % top left corner
                ae = De + max(-Fe,0);
                aw = 0;
                an = 0;
                as = Ds + max(Fs,0);
                Sp = -(2*Dw + 2*Dn);
                Fn = 0;
                Fw = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs)- Sp);
                B(n,3) = ap;
                B(n,1) = -thetha * as;
                B(n,4) = -thetha*ae;
            else % left face except at the corners
                ae = De + max(-Fe,0);
                aw = 0;
                an = Dn + max(-Fn,0);
                as = Ds + max(Fs,0);
                Sp = -(2*Dw);
                Fw = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                B(n,3) = ap;
                B(n,1) = -thetha*as;
                B(n,5) = -thetha * an;
                B(n,4) = -thetha*ae;
            end
            
        elseif(i == N)
            % right face
            if(j == 1) % bottom right corner
                ae = 0;
                aw = Dw + max(Fw,0);
                an = Dn + max(-Fn,0);
                as = 0;
                Sp = -(2*De + 2*Ds);
                Fe = 0;
                Fs = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                B(n,3) = ap;
                B(n,2) = -thetha*aw;
                B(n,5) = -thetha * an;
            elseif(j == N) % top right corner
                ae = 0;
                aw = Dw + max(Fw,0);
                an = 0;
                as = Ds + max(Fs,0);
                Sp = -(2*De + 2*Dn);
                Fe = 0;
                Fn = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                B(n,3) = ap;
                B(n,2) = -thetha*aw;
                B(n,1) = -thetha*as;
            else % right face except at the corners
                ae = 0;
                aw = Dw + max(Fw,0);
                an = Dn + max(-Fn,0);
                as = Ds + max(Fs,0);
                Sp = -(2*De);
                Fe = 0;
                ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                B(n,3) = ap;
                B(n,2) = -thetha*aw;
                B(n,1) = -thetha*as;
                B(n,5) = -thetha*an;
            end
            
        elseif(j == 1) % bottom face except at the corners
            ae = De + max(-Fe,0);
            aw = Dw + max(Fw,0);
            an = Dn + max(-Fn,0);
            as = 0;
            Sp = -(2*Ds);
            Fs = 0;
            ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
            B(n,3) = ap;
            B(n,2) = -thetha*aw;
            B(n,4) = -thetha*ae;
            B(n,5) = -thetha*an;
        elseif(j == N) % top face except at the corners
            ae = De + max(-Fe,0);
            aw = Dw + max(Fw,0);
            an = 0;
            as = Ds + max(Fs,0);
            Sp = -(2*Dn);
            Fn = 0;
            ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
            B(n,3) = ap;
            B(n,2) = -thetha*aw;
            B(n,4) = -thetha*ae;
            B(n,1) = -thetha*as;
        else % interior points
            ae = De + max(-Fe,0);
            aw = Dw + max(Fw,0);
            an = Dn + max(-Fn,0);
            as = Ds + max(Fs,0);
            Sp = 0;
            ap = ap_not + thetha*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
            B(n,3) = ap;
            B(n,2) = -thetha*aw;
            B(n,4) = -thetha*ae;
            B(n,5) = -thetha * an;
            B(n,1) = -thetha*as;
        end
        X(((j - 1)*N+i)) = phi(Point_x,Point_y);
        coeff(n,1) = as;
        coeff(n,2) = aw;
        coeff(n,3) = ap;
        coeff(n,4) = ae;
        coeff(n,5) = an;
        n = n + 1;
    end
    
end


%% Iteration begins:

Y2 = X;
diff = 20;
t = 0;
it = 0;

while (diff > 10^(-6))
    C1 = zeros(N*N,1);
    n = 1;
    
    for j = 1:N
        for i = 1:N
            %% Compuatation of the right hand Side.
            east = (i)*delta_x;
            East = delta_x/2 + (i)*delta_x;
            west = (i - 1)*delta_x;
            West = delta_x/2  + (i - 2)*delta_x;
            north = j*delta_y;
            North = delta_y/2 + (j)*delta_y;
            south = (j - 1)*delta_y;
            South = delta_y/2 + (j - 2)*delta_y;
            Point_x = delta_x/2 + (i - 1)*delta_x;
            Point_y = delta_y/2 + (j - 1)*delta_y;
            %% left face
            Fe = (u(east,Point_y))*delta_y;
            Fw = (u(west,Point_y))*delta_y;
            Fn = (v(Point_x,north))*delta_x;
            Fs = (v(Point_x,south))*delta_x;
            left = (j - 1)*N + i - 1;
            right = (j - 1)*N + i + 1;
            top = j*N + i;
            point = (j - 1)*N + i;
            bottom = (j - 2)*N + i;
            if (i == 1)
                if (j == 1) %left bottom corner
                    aw = 0;
                    an = Dn + max(-Fn,0);
                    as = 0;
                    ae = De + max(-Fe,0);
                    Sp = -(2*Dw + 2*Ds);
                    Fw = 0;
                    Fs = 0;
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*ae*X(right) +(1 - thetha)*an*X(top) + X(point)*ap_prime;
                elseif (j == N) % top left corner
                    ae = De + max(-Fe,0);
                    aw = 0;
                    an = 0;
                    as = Ds + max(Fs,0);
                    Sp = -(2*Dw + 2*Dn);
                    Fn = 0;
                    Fw = 0;
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*ae*X(right) + (1 - thetha)*as*X(bottom)+ X(point)*ap_prime;
                else % left face except at the corners
                    ae = De + max(-Fe,0);
                    aw = 0;
                    an = Dn + max(-Fn,0);
                    as = Ds + max(Fs,0);
                    Sp = -(2*Dw);
                    Fw = 0;
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*ae*X(right) +(1 - thetha)*an*X(top) + (1 - thetha)*as*X(bottom)+ X(point)*ap_prime;
                end
                
            elseif(i == N)
                %% right face
                if(j == 1) % bottom right corner
                    ae = 0;
                    aw = Dw + max(Fw,0);
                    an = Dn + max(-Fn,0);
                    as = 0;
                    Sp = -(2*De + 2*Ds);
                    Fe = 0;
                    Fs = 0;
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*aw*X(left) +(1 - thetha)*an*X(top)+ X(point)*ap_prime;
                elseif(j == N) % top right corner
                    ae = 0;
                    aw = Dw + max(Fw,0);
                    an = 0;
                    as = Ds + max(Fs,0);
                    Sp = -(2*De + 2*Dn);
                    Fe = 0;
                    Fn = 0;
                    
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*aw*X(left)+ (1 - thetha)*as*X(bottom)+ X(point)*ap_prime;
                else % right face except at the corners
                    ae = 0;
                    aw = Dw + max(Fw,0);
                    an = Dn + max(-Fn,0);
                    as = Ds + max(Fs,0);
                    Sp = -(2*De);
                    Fe = 0;
                    ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                    C1(n,1) = (1 - thetha)*aw*X(left)+(1 - thetha)*an*X(top) + (1 - thetha)*as*X(bottom) + X(point)*ap_prime;
                end
                
            elseif(j == 1) % bottom face except at the corners
                ae = De +max(-Fe,0);
                aw = Dw + max(Fw,0);
                an = Dn + max(-Fn,0);
                as = 0;
                Sp = -(2*Ds);
                Fs = 0;
                ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                C1(n,1) = (1 - thetha)*ae*X(right) +(1 - thetha)*aw*X(left)+(1 - thetha)*an*X(top) + X(point)*ap_prime;
            elseif(j == N) % top face except at the corners
                ae = De + max(-Fe,0);
                aw = Dw + max(Fw,0);
                an = 0;
                as = Ds + max(Fs,0);
                Sp = -(2*Dn);
                Fn = 0;
                ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                C1(n,1) = (1 - thetha)*ae*X(right) +(1 - thetha)*aw*X(left)+ (1 - thetha)*as*X(bottom)+ X(point)*ap_prime;
            else % interior points
                ae = De + max(-Fe,0);
                aw = Dw + max(Fw,0);
                an = Dn + max(-Fn,0);
                as = Ds + max(Fs,0);
                Sp = 0;
                ap_prime = ap_not - (1 - thetha)*(aw + an + as + ae + (Fe - Fw) + (Fn - Fs) - Sp);
                C1(n,1) = (1 - thetha)*ae*X(right) +(1 - thetha)*aw*X(left)+(1 - thetha)*an*X(top) + (1 - thetha)*as*X(bottom) + X(point)*ap_prime;
            end
            coeff(n,1) = as;
            coeff(n,2) = aw;
            coeff(n,3) = ap;
            coeff(n,4) = ae;
            coeff(n,5) = an;
            
            
            n = n + 1;
        end
    end
    %% Gauss Seidel Method without relaxation
    if (option == 1)
        [Y2,it] = MA14M004_Gauss(B,C1,Y2,it,1);
    end
    %% Gauss Seidel Method with relaxation
    if (option == 2)
        [Y2,it] = MA14M004_Gauss(B,C1,Y2,it,1.8);
    end
    %% Conjugate Gradient Method
    if (option == 3)
        [Y2,it] = MA14M004_CG(B,C1,Y2,it);
    end
    %% Biconjugate Gradient Stabilized Method
    if (option == 4)
        [Y2,it] = MA14M004_BCG(B,C1,Y2,it);
    end
    
    diff = (abs(X(49*100 + 50) - Y2(49*100 + 50)));
    
    t = t + delta_t;
    X = Y2;
end

l = 1;
%% Final result
D = zeros(N,N);
for j = 1:N
    for i = 1:N
        D(i,j) = Y2(l);
        
        l = l  + 1;
    end
end

D = flipdim(D,1);
t1 = toc;
t2 = cputime - t2;
end


