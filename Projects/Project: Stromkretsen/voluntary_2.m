clear, clc, close all
format long

%%%%%%%%%%%%%%%%% ~ VARIABLES ~ %%%%%%%%%%%%%%%%%
L_0        = 1;
U_0        = [240,1200,2400];
I_0        = [zeros(1,length(U_0));U_0./L_0];
C          = 1e-6;
N          = 1e6;
t_period   = [0, 0.01];
options    = odeset('RelTol',1e-9);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% ~ System of first order ODE ~ %%%%%%%%
dUI = @(t,y) [y(2)*(1+y(1)^2)/L_0; -y(1)/C];

for i = 1:length(U_0)
    I0 = I_0(:,i)';
    %computation of ODEs
    [t_ode, I_ode] = ode45(dUI, t_period, I0, options);
    [t_rk4, I_rk4] = rk4(dUI, t_period, N, I_0(:,i));
    
    %plotting figure
    figure(i)
    tiledlayout(6,6);
    ax1 = nexttile([3 6]);
    plot(ax1,t_ode, I_ode(:,1), '-b', t_rk4, I_rk4(1,:), '-.r', t_ode, t_ode*0, 'k')
    grid on
    legend('ode45', 'rk4', 'Location','best')
    title(sprintf('Current at U_0 = %0.i V',U_0(i)))
    xlabel('Time [seconds]')
    ylabel('Current [Ampere]')
    
    figure(i)
    ax2 = nexttile([3 6]);
    plot(ax2,t_ode, I_ode(:,2), '-b', t_rk4, I_rk4(2,:), '-.r', t_ode, t_ode*0, 'k')
    grid on
    legend('ode45', 'rk4', 'Location','best')
    title(sprintf('Voltage at U_0 = %0.i V',U_0(i)))
    xlabel('Time [seconds]')
    ylabel('Voltage [V]')


end

%%%%%%%%%%%%% ~ I_MAX + T_PERIOD ~ %%%%%%%%%%%%%%%
N = 1e2;
j = 1;
warning('off')
while(true)
    for i = 1:length(U_0)
        %computation of ODEs
        [t_rk4, I_rk4] = rk4(dUI, t_period, N, I_0(:,i));
        %maximum amplitude of current, and its time (index)
        [I_rk4_max, I_rk4_t] = max(I_rk4(1,1:(end/5)));  
        %store value
        I_MAX(j,i) = I_rk4_max;
        T_PERIOD(j,i) = 4*I_rk4_t*t_period(2)/N;
    end
    
    j = j + 1;
    N = N * 2;
    
    if N >= 1e6
        break
    end
end
warning('on')

% absolute error estimation
for k = 1:length(I_MAX)-1
    for l=1:3 
        t_error(k,l) = abs(T_PERIOD(k+1,l)-T_PERIOD(k,l));
        I_error(k,l) = abs(I_MAX(k+1,l)-I_MAX(k,l));
    end
end


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