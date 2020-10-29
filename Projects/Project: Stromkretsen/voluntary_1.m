clear, clc, close all
format longe

%%%%%%%%%%%%%%%%% ~ TO DO ~ %%%%%%%%%%%%%%%%%%%%%
% relationship between U_0 and L_0
%%%%%%%%%%%%%%%%% ~ VARIABLES ~ %%%%%%%%%%%%%%%%%
L_0        = 1;
U_0        = 240;
I_0        = [zeros(1,length(U_0));U_0./L_0];
C          = 1e-6;
N          = 1e6;
t_period   = [0, 0.002];
options    = odeset('RelTol',1e-7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% ~ INITIALIZATION ~ %%%%%%%%%%%%%%%
I_ode_max  = [];
I_ode_max1 = [];
I_rk4_max  = [];
I_rk4_t    = [];
t_rk4      = [];
t_ode      = [];
I_ode      = [];
I_ode_t    = [];
results    = zeros(20,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% ~ System of first order ODE ~ %%%%%%%%
dI = @(t,y) [y(2); 2*y(1)/(1+y(1)^2)*y(2)^2-y(1)*((1+y(1)^2)/(L_0*C))];

%%%%%%%%%% ~ COMPUTATION OF U_MAX ~ %%%%%%%%%%%%%%
b = 2150;
a = 2140;
i = 1;
while(true)
%     computation of ODEs
    c = (a+b)/2;
    [~, Ia] = rk4(dI, t_period, N, [0; a]);
    [~, Ic] = rk4(dI, t_period, N, [0; c]);
    Ira = 10-max(Ia(1,:));
    Irc = 10-max(Ic(1,:));
    if Ira*Irc > 0
        a = c;
    else
        b = c;
    end
    U_max(i) = c;
    if abs((b-a)/2) < 0.5e-12
        break
    end
    i = i + 1;
end
%plotting at I_max = 10
[ta, Iu] = rk4(dI, [0 0.01], N, [0; U_max(end)]);
plot(ta,Iu(1,:))
grid on
title(sprintf('Current at U_0 = %0.i',U_0))
xlabel('Time [seconds]')
ylabel('Current [Ampere]')
ylim([-13 13])

%%%%%%%%%%%%%% ~ L_0 AND U_max  ~ %%%%%%%%%%%%%%
% calculating the relationship between L_0 and U_max
b          = 2400;
a          = 240;
for i = 1:20
    L_0 = 0.1 * i;
    b   = 2400;
    a   = 240;
    while(true)
        c       = (a+b)/2;
        [~, Ia] = rk4(dI, [0; 0.01], N, [0; a/L_0]);
        [~, Ic] = rk4(dI, [0; 0.01], N, [0; c/L_0]);
        Ira      = 10-max(Ia(1,:));
        Irc      = 10-max(Ic(1,:));        
        if Ira*Irc > 0
            a = c;
        else
            b = c;
        end
        if abs((b-a)/2) < 0.5e-9
            break
        end
    end
    results(i, 2) = c;
    results(i, 1) = L_0;
end


plot(results(:,1),results(:,2))
title(sprintf('relationship between U*_0 and L_0 = %0.i',U_0))
xlabel('L_0 [H]')
ylabel('Power [Volt]')
ylim([-13 13])
grid on


%%%%%%%%%%%%% ~ LOCAL FUNCTIONS ~ %%%%%%%%%%%%%%
function [t,y] = rk4(f, tspan, N, y0 )
h=(tspan(end)-tspan(1))/N;
t=tspan(1):h:tspan(end);
y=y0;
    for i=1:N    
        f1=f(t(i),y(:,i));    
        f2=f(t(i)+h/2,y(:,i)+h/2*f1);    
        f3=f(t(i)+h/2,y(:,i)+h/2*f2);    
        f4=f(t(i)+h,y(:,i)+h*f3);    
        y(:,i+1)=y(:,i)+h/6*(f1+2*f2+2*f3+f4);
    end
end
