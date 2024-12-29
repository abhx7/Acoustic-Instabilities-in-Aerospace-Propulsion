N =1000; %number of timesteps
%T = 40; %total time of experiment %fig4a,4b,4c
T = 80; %total time of experiment %fig5a,5b,5c,5d fig6a,6b,6c,6d
dt = (T/N);

%define system parameters
g = 1.4;
R = 287;
Lw= 3.6; %equivalent length of wire
lambda = 0.0328; %heat conductivity of air
Cv = 719; %specific heat of air
pbar = 101325; %mean pressure of air
rho = 1.225; %mean density of air
Tbar = pbar/(rho*R); %mean temperature of air
Tw = 1000; 
dT = abs(Tw - Tbar); %temperature difference
%dw = 0.1; %diameter of the wire
%S = 1; %cross section area of tube
c0 = 399.6; %speed of sound
u0 = 0.5; %mean velocity of air
M = u0/c0;
%K = (g-1)*((2*Lw*dT)/(S*c0*pbar*sqrt(3)))*sqrt(pi*lambda*Cv*rho*dw*u0/2);


%K = 0.1 %fig4a
%K = 0.6 %fig4b,4d
K = 4.5 %fig5a,5b,5c,5d
%K = 0.1 %fig6a,6b,6c,6d

%xf = 0.29; %fig4a,4b,4d
xf = 0.25; %fig5a,5b,5c,5d fig6a,6b,6c,6d

%tau = 0.2; %fig4a
%tau = 0.5; %fig4b,4d
tau = 0.45*pi; %fig5a,5b,5c,5d fig6a,6b,6c,6d

% Set initial values
u = zeros(1, N);
p = zeros(1, N);

J = 3; %number of Galerkin modes

% Preallocate arrays
y1 = zeros(J, N);
y2 = zeros(J, N);

for j = 1:J
    % Set initial values
    y1(j,1) = 0;
    %y1(1,1) = 0.15 %fig4a
    %y1(1,1) = 0.2; %fig4b,4d
    y1(1,1) = 0.18; %fig5a,5b,5c,5d fig6a,6b,6c,6d
    y2(j,1) = 0;
    u(1) = u(1) + y1(j,1)*cos(j*pi*xf);
    p(1) = p(1) + y2(j,1)*((-g*M)/(j*pi))*sin(j*pi*xf);
        
    n1 = round(tau/dt);
    % Fourth-order Runge-Kutta method
    for n = 1:n1     
        [k1, l1] = equations1(y1(j,n), y2(j,n),j);
        [k2, l2] = equations1(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j);
        [k3, l3] = equations1(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j);
        [k4, l4] = equations1(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u(n+1) = u(n+1) + y1(j,n+1)*cos(j*pi*xf);
        p(n+1) = p(n+1) + y2(j,n+1)*((-g*M)/(j*pi))*sin(j*pi*xf);
    end

    for n = n1+1:N-1
        [k1, l1] = equations2(y1(j,n), y2(j,n),j,K,u(n-n1),xf);
        [k2, l2] = equations2(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j,K,u(n-n1),xf);
        [k3, l3] = equations2(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j,K,u(n-n1),xf);
        [k4, l4] = equations2(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j,K,u(n-n1),xf);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u(n+1) = u(n+1) + y1(j,n+1)*cos(j*pi*xf);
        p(n+1) = p(n+1) + y2(j,n+1)*((-g*M)/(j*pi))*sin(j*pi*xf);
    end

end       

figure(6)
tiledlayout(2,2)
nexttile
plot(linspace(0,T,N), u, 'b','linewidth',2);
title('Variation of u''');
xlabel('t (s)');
ylabel('u'' ');
nexttile
plot(linspace(0,T,N), y1(1,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('\eta_1');
title('Variation of \eta_1 (t)');
nexttile
plot(linspace(0,T,N), y1(2,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('\eta_2');
title('Variation of \eta_2 (t)');
nexttile
plot(linspace(0,T,N), y1(3,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('\eta_3');
title('Variation of \eta_3 (t)');


figure(1);
plot(linspace(0,T,N), u, 'b','linewidth',2);
title('Variation of u''');
xlabel('t (s)');
ylabel('u'' ');
fontsize(gcf,scale=1.5)
grid on;

e=((0.5*(p.^2)+0.5*(g*M*u).^2)/((g*M)^2));
ee = envelope(e,100,"peak");
figure(5);
plot(linspace(0,T,N), ee, 'r','linewidth',2);
title('Energy Fluctuations Amplitude');
xlabel('t (s)');
ylabel('Energy/(\gamma M)^2 ');
legend('Non Linear')
fontsize(gcf,scale=1.5)
grid on;

figure(2);
plot(linspace(0,T,N), y1(1,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('\eta_1');
title('Variation of \eta_1 (t)');
fontsize(gcf,scale=1.5)
grid on;

figure(3);
plot(linspace(0,T,N), y1(2,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('Variation of \eta_2');
title('\eta_2 (t)');
fontsize(gcf,scale=1.5)
grid on;

figure(4);
plot(linspace(0,T,N), y1(3,:), 'g','linewidth',2);
xlabel('t (s)');
ylabel('\eta_3');
title('Variation of  \eta_3 (t)');
%ylim([-0.5 0.5])%fig6d
fontsize(gcf,scale=1.5)
grid on;


% Define the functions f1(x, y1, y2) and f2(x, y1, y2) 
% where y1 is eta and y2 is eta dot

% for t < tau
function [f1,f2] = equations1(y1, y2, j)
    kj = j*pi;
    wj = kj;
    w1 = pi;
    c1 = 0.1; c2 = 0.06;
    zetaj = (1/(2*pi))*(c1*(wj/w1) + c2*sqrt(w1/wj));
    
    f1 =  y2;
    f2 =  - (kj^2) * y1 -2*zetaj*wj*y2;
end

% for t > tau
function [f1,f2] = equations2(y1, y2, j, K, u,xf)
    kj = j*pi;
    wj = kj;
    w1 = pi;
    c1 = 0.1; c2 = 0.06;
    zetaj = (1/(2*pi))*(c1*(wj/w1) + c2*sqrt(w1/wj));

    f1 =  y2;
    f2 = - (kj^2) * y1 -2*zetaj*wj*y2 - ((2*j*pi*K))*(sqrt(abs((1/3) + u)) - sqrt(1/3))*sin(j*pi*xf) ;
end

