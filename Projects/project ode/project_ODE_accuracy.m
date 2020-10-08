clear,clc

h = 1e-1;

t = 0:h:120;
N = length(t);

%%%%% Constants
Gamma = 1.36;
My = 1.36e-3;
Tau = 0.2;
Beta = 0.00027;
Rho = 0.1;
Alpha = 3.6e-2;
Sigma = 2;
Delta = 0.33;
Pi = 100;

%%%%% Initial value conditions
R = 2e2;
L = 0;
E = 0;
V = 4e-7;

% runge-kutta 4 algorithm
sR = zeros(1,4); 
sL = zeros(1,4);  
sE = zeros(1,4); 
sV = zeros(1,4);  
rk4 = [1 2 2 1];   

% System of ODE
dR = @(R,L,E,V) Gamma*Tau - My*R - Beta*R*V;
dL = @(R,L,E,V) Rho*Beta*R*V - My*L - Alpha*L;
dE = @(R,L,E,V) (1-Rho)*Beta*R*V + Alpha*L - Delta*E;
dV = @(R,L,E,V) Pi*E - Sigma*V;

for i = 1:(N-1)
    sR = dR( R, L, E, V);
    sL = dL( R, L, E, V);
    sE = dE( R, L, E, V);
    sV = dV( R, L, E, V);

    sR(2) = dR( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
    sL(2) = dL( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
    sE(2) = dE( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
    sV(2) = dV( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));

    sR(3) = dR( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
    sL(3) = dL( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
    sE(3) = dE( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
    sV(3) = dE( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));

    sR(4) = dR( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
    sL(4) = dL( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
    sE(4) = dE( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
    sV(4) = dE( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));

    R = R + (h/6)*sum(rk4.*sR);       
    L = L + (h/6)*sum(rk4.*sL);       
    E = E + (h/6)*sum(rk4.*sE);
    V = V + (h/6)*sum(rk4.*sV);
    LC = 1000*(1-Tau)+R+L+E;   
end


% % %comparison fine grid with halving step-length
% L1 = LC;
% LC = 0;
% h = 1;
% j = 1;
% while(true)
%     L2(j) = gridapprox(h);
%     h = h/2;  
%     
%     LFG(j) = abs(L1-L2(j));
%     
%     if LFG(j) < 1e-5
%         fine_grid = LFG;
%         break
%     end
%     
%     j = j + 1
% end
% fine_grid'



% %comparison n with 2n
% h = 1;
% j = 1;
% while(true)
%     L2(j) = gridapprox(h);
%     h = h/2;  
%     
%     if j > 1
%         LFG(j-1) = abs(L2(j)-L2(j-1));
%     end
%     
%     if j > 1 && LFG(j-1) < 1e-5
%         break
%     end
%     
%     j = j + 1
% end
% LFG'


% % comparison with ode
% function
% f = @(t,y) [Gamma*Tau - My*y(1) - Beta*y(1)*y(4); Rho*Beta*y(1)*y(4) - My*y(2) - Alpha*y(2) ; (1-Rho)*Beta*y(1)*y(4) + Alpha*y(2) - Delta*y(3) ; Pi*y(3) - Sigma*y(4)];
% setting tolerance value
% options = odeset('RelTol',1e-5);
% 
% [t1,xa1] = ode113(f,[0 120],[200 0 0 4e-7]);
% LC_ode = 1000*(1-Tau)+xa1(end,1)+xa1(end,2)+xa1(end,3);
% h = 1;
% while(true)
%     uh = gridapprox(h);
%     h = h/10
%     diff = abs(LC_ode-uh)
%     if diff < 1e-5
%         break
%     end
% end

function LC_TOT = gridapprox(h)
t = 0:h:120;
N = length(t);

Gamma = 1.36;
My = 1.36e-3;
Tau = 0.2;
Beta = 0.00027;
Rho = 0.1;
Alpha = 3.6e-2;
Sigma = 2;
Delta = 0.33;
Pi = 100;

%%%%% Initial value conditions
R = 2e2;
L = 0;
E = 0;
V = 4e-7;
LC = 1000*(1-Tau)+R+L+E;
        

% runge-kutta 4 algorithm
sR = zeros(1,4);
sL = zeros(1,4);
sE = zeros(1,4);
sV = zeros(1,4);
rk4 = [1 2 2 1];
    
% System of ODE
dR = @(R,L,E,V) Gamma*Tau - My*R - Beta*R*V;
dL = @(R,L,E,V) Rho*Beta*R*V - My*L - Alpha*L;
dE = @(R,L,E,V) (1-Rho)*Beta*R*V + Alpha*L - Delta*E;
dV = @(R,L,E,V) Pi*E - Sigma*V;
    
    for i = 1:(N-1)
        sR = dR( R, L, E, V);
        sL = dL( R, L, E, V);
        sE = dE( R, L, E, V);
        sV = dV( R, L, E, V);
        
        sR(2) = dR( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
        sL(2) = dL( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
        sE(2) = dE( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
        sV(2) = dV( R + (h/2)*sR(1), L + (h/2)*sL(1), E + (h/2)*sE(1),  V + (h/2)*sV(1));
        
        sR(3) = dR( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
        sL(3) = dL( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
        sE(3) = dE( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
        sV(3) = dE( R + (h/2)*sR(2), L + (h/2)*sL(2), E + (h/2)*sE(2),  V + (h/2)*sV(2));
        
        sR(4) = dR( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
        sL(4) = dL( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
        sE(4) = dE( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
        sV(4) = dE( R + h*sR(3), L + h*sL(3), E + h*sE(3),  V + h*sV(3));
        
        R = R + (h/6)*sum(rk4.*sR);
        L = L + (h/6)*sum(rk4.*sL);
        E = E + (h/6)*sum(rk4.*sE);
        V = V + (h/6)*sum(rk4.*sV);
        LC = 1000*(1-Tau)+R+L+E;
    end
    LC_TOT = LC;
end