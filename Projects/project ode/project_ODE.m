clear, clc, close('all')

%%%%% CONSTANTS
gamma = 1.36;
my = 1.36e-3;
tau = 0.2;
Beta = 0.00027;
rho = 0.1;
alpha = 3.6e-2;
sigma = 2;
delta = 0.33;
not_pi = 100;

%step length
h = 1e-5;
time = 0;
X = 0:h:120;

%mem alloc
Y1 = zeros(1,length(X));
Y2 = zeros(1,length(X));
R = zeros(1,length(X));
L = zeros(1,length(X));
E = zeros(1,length(X));
V = zeros(1,length(X));
LC = zeros(1,length(X));
%%%%% STARTING VALUES
R(1) = 2e2;
L(1) = 0;
E(1) = 0;
V(1) = 1e2;
LC(1) = 1000*(1-tau)+R(1)+L(1)+E(1);

%euler-forward algorithm
for i = 1:length(X)-1
    %y1 --start value
    Y1 = [ R(i); L(i); E(i); V(i) ];
    %f(t,y(t))
    Yprime = [ gamma*tau - my*R(i) - Beta*R(i)*V(i); rho*Beta*R(i)*V(i) - my*L(i) - alpha*L(i) ; (1-rho)*Beta*R(i)*V(i) + alpha*L(i) - delta*E(i) ; not_pi*E(i) - sigma*V(i)  ];
    %euler one step
    Y2 = Y1 + Yprime*h; 
    %assign new start variables
    R(i+1) = Y2(1);
    L(i+1) = Y2(2);
    E(i+1) = Y2(3);
    V(i+1) = Y2(4);
    LC(i+1) = 1000*(1-tau)+R(i)+L(i)+E(i);
end

%plot1
t = tiledlayout(1,2);
yyaxis left
ax1 = nexttile([2 1]);
plot(ax1,X,LC)
ylabel('CD4 lymphocytes')
axis(ax1, [0 120 0 1200])
% xticklabels({'0','30','60','90','120'})
yyaxis right
semilogy(ax1,X,V)
ylabel('Free virions')
axis(ax1, [0 120 1e-1 1e4])
yticklabels({'0','0.1','10','100','1000','10,000'})
xlabel('Days from infection')

%plot2
ax2 = nexttile([2 1]);
yyaxis left
axis(ax2, [0 120 0 250])
plot(X,R)
ylabel('R')
axis(ax2, [0 120 0 250])
yyaxis right
semilogy(X,L,X,E)
ylabel('L and E')
axis(ax2, [0 120 1e-1 1e2])
yticklabels({'0','0.1','10','100'})
xlabel('Days from infection')

od45 = @(t,y) [gamma*tau - my*y(1) - Beta*y(1)*y(4); rho*Beta*y(1)*y(4) - my*y(2) - alpha*y(2) ; (1-rho)*Beta*y(1)*y(4) + alpha*y(2) - delta*y(3) ; not_pi*y(3) - sigma*y(4)];
[t1,xa1] = ode45(od45,[0 120],[200 0 0 100]);
ode45_partitions = length(t1)

od23 = @(t,y) [gamma*tau - my*y(1) - Beta*y(1)*y(4); rho*Beta*y(1)*y(4) - my*y(2) - alpha*y(2) ; (1-rho)*Beta*y(1)*y(4) + alpha*y(2) - delta*y(3) ; not_pi*y(3) - sigma*y(4)];
[t2,xa2] = ode23(od45,[0 120],[200 0 0 100]);
ode23_partitions = length(t2)

% 0.5 gigabyte program  l m f a o
