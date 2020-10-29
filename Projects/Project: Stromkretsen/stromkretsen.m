clear, clc, close all
format long

%%%%%%%%%%%%%%%%% ~ TO DO ~ %%%%%%%%%%%%%%%%%%%%%
% relationship between U_0 and L_0
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
dI = @(t,y) [y(2); 2*y(1)/(1+y(1)^2)*y(2)^2-y(1)*((1+y(1)^2)/(L_0*C))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(U_0)
    I0 = I_0(:,i)';
    %computation of ODEs
    [t_ode, I_ode] = ode45(dI, t_period, I0, options);
    [t_rk4, I_rk4] = rk4(dI, t_period, N, I_0(:,i));
    %computation of E(t) to prove its constant
    dEdt = C*(L_0./(1+I_rk4(1,:).^2).*I_rk4(2,:)).^2+L_0*log(1+I_rk4(1,:).^2);
    
    %plotting figure
    figure(i)
    tiledlayout(6,6);
    ax1 = nexttile([4 6]);
    plot(ax1,t_ode, I_ode(:,1), '-b', t_rk4, I_rk4(1,:), '-.r', t_ode, t_ode*0, 'k')
    grid on
    legend('ode45', 'rk4', 'Location','best')
    title(sprintf('Current at U_0 = %0.i',U_0(i)))
    xlabel('Time [seconds]')
    ylabel('Current [Ampere]')
    
    ax2 = nexttile([2 6]);
    plot(ax2, t_rk4, dEdt)
    ylabel('Energy [Wh]')
    title(sprintf('E at U_0 = %0.i',U_0(i)))
    xlabel('Time [seconds]')
    axis manual
    axis(ax2, [0 0.1 0 0.01])
    grid on
    xlim([0 0.001])
    ylim ([0 2*max(dEdt)])
    axis auto
end


%%%%%%%%%%%%%%%% ~ FOURIER ~ %%%%%%%%%%%%%%%%%%
I_fourier = 0;
N = 1e6;
% obtained from accuracy.m
period_index = [ ceil(N/(10^8/31190145)) , ceil(N/(10^8/26085965)),  ceil(N/(10^8/15470926)) ];
period = [ 0.006238029174805   0.005217193603516   0.003094185791016 ];

for i = 1:length(U_0)
    
    %computation of ODE
    [t, I] = rk4(dI, t_period, N, I_0(:,i));
    I = I';
    I = I(:,1);
    t = t';
    
    %calculation of omega
    w = 2 * pi / period(i);
    
    %conversion of h to match period
    t_p = t(1:period_index(i));
    TN = period_index(i);
    h = period(i) / TN;
    x = 1:TN;
    
    for k = 1:14
    sine_0 = sin(k*w.*t(x(1)));
    sine_n = sin(k*w.*t(x(end)));
    sine   = sin(k*w.*t(x(2:end-1))); 
    a(k) = (2/period(i)) * (h/2) * sum(I(x(1))*sine_0  ...
                                 +     I(end)*sine_n  ...
                                 +     I(x(2:end-1)).*sine*2);
    end
    ak = a';
    I_fourier = @(t) ak(1)*sin(1*w*t) + ak(2)*sin(2*w*t) + ak(3)*sin(3*w*t) + ak(4)*sin(4*w*t) + ak(5)*sin(5*w*t) + ak(6)*sin(6*w*t) + ak(7)*sin(7*w*t) + ak(8)*sin(8*w*t) + ak(9)*sin(9*w*t) + ak(10)*sin(10*w*t) + ak(11)*sin(11*w*t) + ak(12)*sin(12*w*t) + ak(13)*sin(13*w*t) + ak(14)*sin(14*w*t);
    I_fourier3 = @(t) ak(1)*sin(1*w*t) + ak(2)*sin(2*w*t) + ak(3)*sin(3*w*t);
    I_fourier_val = I_fourier(t);
    I_fourier_val3 = I_fourier3(t);

    %plot
    figure(i+3)   
    plot(t,I(:,1), '-.r',t,I_fourier_val,'b')
    legend('rk4','fourier')
    title(sprintf('Fourier vs Runge U_0 = %0.i',U_0(i)))
    xlabel('Time [seconds]')
    ylabel('Current [Ampere]')
    grid on
    figure(i+6)
    plot(t,I_fourier_val3,'m--',t,I_fourier_val,'b')
    legend('fourier with 3 a_k','fourier with 14 a_k')
    title(sprintf('Fourier U_0 = %0.i',U_0(i)))
    xlabel('Time [seconds]')
    ylabel('Current [Ampere]')
    grid on
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
