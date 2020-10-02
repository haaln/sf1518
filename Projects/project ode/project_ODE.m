clear, clc, close('all')

%%%%% CONSTANTS
Gamma = 1.36;
My = 1.36e-3;
Tau = 0.2;
Beta = 0.00027;
Rho = 0.1;
Alpha = 3.6e-2;
Sigma = 2;
Delta = 0.33;
Pi = 100;

%step length
h = 1e-3;
time = 0;
t = 0:h:120;
N = length(t);

%mem alloc
Y1 = zeros(1,N);
Y2 = zeros(1,N);
R = zeros(1,N);
L = zeros(1,N);
E = zeros(1,N);
V = zeros(1,N);
LC = zeros(1,N);

%%%%% STARTING VALUES
R(1) = 2e2;
L(1) = 0;
E(1) = 0;
V(1) = 100;
LC(1) = 1000*(1-Tau)+R(1)+L(1)+E(1);

% %euler-forward algorithm
% for i = 1:(N-1)
%     %y1 --start value
%     Y1 = [ R(i); L(i); E(i); V(i) ];
%     %f(t,y(t))
%     Yprime = [ Gamma*Tau - My*R(i) - Beta*R(i)*V(i); Rho*Beta*R(i)*V(i) - My*L(i) - Alpha*L(i) ; (1-Rho)*Beta*R(i)*V(i) + Alpha*L(i) - Delta*E(i) ; Pi*E(i) - Sigma*V(i)  ];
%     %euler one step
%     Y2 = Y1 + Yprime*h; 
%     %assign new start variables
%     R(i+1) = Y2(1);
%     L(i+1) = Y2(2);
%     E(i+1) = Y2(3);
%     V(i+1) = Y2(4);
%     LC(i+1) = 1000*(1-Tau)+R(i)+L(i)+E(i);
% end


% euler-backward algorithm
% for i = 1:(N-1)
%     Y2 = [ (R(i)+h*Tau*Gamma)/(1+h*My+(h*Beta)*V(i)); (L(i)+h*Beta*Rho*R(i)*V(i))/(1+My*h+h*Alpha); (E(i)+h*Beta*R(i)*V(i)-h*Beta*Rho*R(i)*V(i)+Alpha*L(i)*h)/(1+Delta*h); (V(i)+Pi*E(i)*h)/(1+Sigma*h)];
%     R(i+1) = Y2(1);
%     L(i+1) = Y2(2);
%     E(i+1) = Y2(3);
%     V(i+1) = Y2(4);
%     LC(i+1) = 1000*(1-Tau)+R(i)+L(i)+E(i);
% end

% % runge-kutta 4 algorithm
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



%plot1

tiledlayout(1,2);
colororder({'k','k'})
yyaxis left
ax1 = nexttile([2 1]);
plot(ax1,t,LC,'k')
text(-7,1150,'A','FontSize',15,'FontWeight','bold')
text(32,890,'CD4 lymphocytes','FontSize',13)
ylabel('CD4 lymphocytes','FontWeight','bold')
axis(ax1, [-10 120 0 1200])
yyaxis right
semilogy(ax1,t,V,'k')
text(20,10,'Cell-free virus','FontSize',13)
ylabel('Free virions V','color','k','FontWeight','bold')
axis(ax1, [-10 120 1e-1 1e4])
yticklabels({'0','0.1','10','100','1000','10,000'})
xlabel('Days from infection','FontWeight','bold')

%plot2
ax2 = nexttile([2 1]);
colororder({'k','k'})
yyaxis left
axis(ax2, [0 120 0 250])
plot(t,R,'k')
text(-7,240,'B','FontSize',15,'FontWeight','bold')
text(58,20,'R','FontAngle','italic','FontSize',13)
ax.YColor = 'k';
ylabel('R','color','k','FontWeight','bold','FontAngle','italic')
axis(ax2, [0 120 0 250])
yyaxis right
semilogy(t,L,'k--',t,E,'k-')
text(58,5,'L','FontAngle','italic','FontSize',13)
text(58,1,'E','FontAngle','italic','FontSize',13)
ylabel('L and E','color','k','FontWeight','bold')
axis(ax2, [-10 120 1e-1 1e2])
yticklabels({'0.1','1','10','100'})
xlabel('Days from infection','FontWeight','bold')

%graph size and position
x0=550;
y0=550;
width=850;
height=450;
set(gcf,'position',[x0,y0,width,height])

%xticks are not set in increments of 30, from 0 to 120
%ylabelticks for both semilogy graph are too small
%ylabel for plot2 "L and E" do not have cursive writing style for variables

%function
f = @(t,y) [Gamma*Tau - My*y(1) - Beta*y(1)*y(4); Rho*Beta*y(1)*y(4) - My*y(2) - Alpha*y(2) ; (1-Rho)*Beta*y(1)*y(4) + Alpha*y(2) - Delta*y(3) ; Pi*y(3) - Sigma*y(4)];
%ode45 explicit
[t1,xa1] = ode45(f,[0 120],[200 0 0 100]);
ode45_partitions = length(t1)
%ode23 implicit
[t2,xa2] = ode23(f,[0 120],[200 0 0 100]);
ode23_partitions = length(t2)

% error_ode45 = abs(L(end)-xa1(end,2))
% error_ode23 = abs(L(end)-xa2(end,2))
% you have not calculated the error to be below 1e-5

