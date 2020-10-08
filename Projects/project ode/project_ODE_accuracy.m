clear,clc

h = 1e-5;

t = 0:h:120;
N = length(t);

% Memory allocation
R = zeros(1,N);
L = zeros(1,N);
E = zeros(1,N);
V = zeros(1,N);
LC = zeros(1,N);

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
R(1) = 2e2;
L(1) = 0;
E(1) = 0;
V(1) = 4e-7;
LC(1) = 1000*(1-Tau)+R(1)+L(1)+E(1);

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
    sR(1) = dR( R(i), L(i), E(i), V(i));
    sL(1) = dL( R(i), L(i), E(i), V(i));
    sE(1) = dE( R(i), L(i), E(i), V(i));
    sV(1) = dV( R(i), L(i), E(i), V(i));

    sR(2) = dR( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
    sL(2) = dL( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
    sE(2) = dE( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
    sV(2) = dV( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));

    sR(3) = dR( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
    sL(3) = dL( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
    sE(3) = dE( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
    sV(3) = dE( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));

    sR(4) = dR( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
    sL(4) = dL( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
    sE(4) = dE( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
    sV(4) = dE( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));

    R(i+1) = R(i) + (h/6)*sum(rk4.*sR);       
    L(i+1) = L(i) + (h/6)*sum(rk4.*sL);       
    E(i+1) = E(i) + (h/6)*sum(rk4.*sE);
    V(i+1) = V(i) + (h/6)*sum(rk4.*sV);
    LC(i+1) = 1000*(1-Tau)+R(i)+L(i)+E(i);   
end

L1 = LC(end);
LC(end) = 0;

h = 1;
j = 1;

while(true)
    L2(j) = gridapprox(h);
    h = h/2;  
    
    if j > 1
        LFG(j-1) = abs(L2(j)-L2(j-1));
    end
    
    if j > 1 && LFG(j-1) < 1e-5
        break
    end
    
    j = j + 1;
end
LFG'
h
t1 = 0:h:120;
N = length(t1);

function Lend = gridapprox(h)
    t = 0:h:120;
    N = length(t);
    % Memory allocation
    R = zeros(1,N);
    L = zeros(1,N);
    E = zeros(1,N);
    V = zeros(1,N);
    LC = zeros(1,N);
    
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
R(1) = 2e2;
L(1) = 0;
E(1) = 0;
V(1) = 4e-7;
LC(1) = 1000*(1-Tau)+R(1)+L(1)+E(1);
        
    %%%%% Initial value conditions
    R(1) = 2e2;
    L(1) = 0;
    E(1) = 0;
    V(1) = 4e-7;
    LC(1) = 1000*(1-Tau)+R(1)+L(1)+E(1);
    
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
        sR(1) = dR( R(i), L(i), E(i), V(i));
        sL(1) = dL( R(i), L(i), E(i), V(i));
        sE(1) = dE( R(i), L(i), E(i), V(i));
        sV(1) = dV( R(i), L(i), E(i), V(i));
        
        sR(2) = dR( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
        sL(2) = dL( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
        sE(2) = dE( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
        sV(2) = dV( R(i) + (h/2)*sR(1), L(i) + (h/2)*sL(1), E(i) + (h/2)*sE(1),  V(i) + (h/2)*sV(1));
        
        sR(3) = dR( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
        sL(3) = dL( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
        sE(3) = dE( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
        sV(3) = dE( R(i) + (h/2)*sR(2), L(i) + (h/2)*sL(2), E(i) + (h/2)*sE(2),  V(i) + (h/2)*sV(2));
        
        sR(4) = dR( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
        sL(4) = dL( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
        sE(4) = dE( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
        sV(4) = dE( R(i) + h*sR(3), L(i) + h*sL(3), E(i) + h*sE(3),  V(i) + h*sV(3));
        
        R(i+1) = R(i) + (h/6)*sum(rk4.*sR);
        L(i+1) = L(i) + (h/6)*sum(rk4.*sL);
        E(i+1) = E(i) + (h/6)*sum(rk4.*sE);
        V(i+1) = V(i) + (h/6)*sum(rk4.*sV);
        LC(i+1) = 1000*(1-Tau)+R(i)+L(i)+E(i);
    end
    Lend = LC(end-1);
end