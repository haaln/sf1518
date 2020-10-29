clear, clc, close all
format shortg

%%%%%%%%%%%%%%%%% ~ VARIABLES ~ %%%%%%%%%%%%%%%%%
L_0        = 1;
U_0        = [240,1200,2400];
I_0        = [zeros(1,length(U_0));U_0./L_0];
C          = 1e-6;
t_period   = [0, 0.01];
options    = odeset('RelTol',1e-9);

N = 1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% ~ INITIALIZATION ~ %%%%%%%%%%%%%%%
I_ode_max   = [];
I_ode_max1  = [];
I_rk4_max   = [];
I_rk4_t     = [];
t_rk4       = [];
t_ode       = [];
I_ode       = [];
I_ode_t     = [];
I_MAX       = [];
T_PERIOD    = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% ~ System of first order ODE ~ %%%%%%%%
dI = @(t,y) [y(2); 2*y(1)/(1+y(1)^2)*y(2)^2-y(1)*((1+y(1)^2)/(L_0*C))];
    
%%%%%%%%% ~ CALCULATING THE PERIOD ~ %%%%%%%%%%
j = 1;
warning('off')
while(true)
    for i = 1:length(U_0)
        %computation of ODEs
        [t_rk4, I_rk4] = rk4(dI, t_period, N, I_0(:,i));
        %maximum amplitude of current, and its time (index)
        [I_rk4_max, I_rk4_t] = max(I_rk4(1,1:(end/5)));  
        %store value
        I_MAX(j,i) = I_rk4_max;
        T_PERIOD(j,i) = 4*I_rk4_t*t_period(2)/N
    end
    
    j = j + 1;
    N = N * 2;
    
    if N >= 1e8
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 2 : bisection method
% the absolute error for this method is 1e-5 (OR IS IT 0.5e-5?)
% a = 2;
% N = 1e5;
% b = round(N/(57/20));
% for i = 1:length(U_0)
%     %computation of ODEs
%     [t, I] = rk4(dI, t_period, N, I_0(:,i));
% 
%     while(true)
%         c = round( (b+a)/2 );
%         if I(1,c)*I(1,a) > 0
%             a = c;
%         else
%             b = c;
%         end
%	  if (b-a)<=2
%	      c = round( (b+a)/2 );
%	      break
%	  end
%         if abs(0.5*I(1,b)-0.5*I(1,a)) < 0.5e-5
%             c = round( (b+a)/2 );
%             break
%         end
% 
%     end
% 
% end
% t_index(i) = c;
% period(i) = 2*c*0.01/N;


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
