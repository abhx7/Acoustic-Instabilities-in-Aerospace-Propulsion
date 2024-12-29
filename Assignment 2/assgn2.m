% Initial conditions
y1_0 = 2000;
y2_0 = 0;

% Spatial domain
L = 4;
N = 41;
x = linspace(0, L, N);
dx = x(2) - x(1);

%define system parameters
n = 9;
L = 4;
gamma = 1.4;
R = 287;
Pbar = 101532;

p=[];u=[];%initialising arrays for plotting 
omega = zeros(1,4);
    

%looping for different linear temperature profiles
for j=1:5
    %different linear temperature profiles
    m = 50-j*50;
    T0 = 100 + j*200;

    pend=1000000000;
    for ij = 1:N
        w = 2*pi*(n*sqrt(gamma*R*(T0 + m*x(ij))))/(4*L);

        % Preallocate arrays
        y1 = zeros(1, N);
        y2 = zeros(1, N);
        uu = zeros(1, N);
    
        % Set initial values
        y1(1) = y1_0;
        y2(1) = y2_0;
        uu(1) = -((4*L*sqrt(R*T0))/(2*pi*n*Pbar*sqrt(1.4)))*y2(1);

        % Fourth-order Runge-Kutta method
        for i = 1:N-1
            [k1, l1] = equations(x(i), y1(i), y2(i),m,T0,n,L,gamma,R,w);
            [k2, l2] = equations(x(i) + 0.5*dx, y1(i) + 0.5*dx*k1, y2(i) + 0.5*dx*l1,m,T0,n,L,gamma,R,w);
            [k3, l3] = equations(x(i) + 0.5*dx, y1(i) + 0.5*dx*k2, y2(i) + 0.5*dx*l2,m,T0,n,L,gamma,R,w);
            [k4, l4] = equations(x(i) + dx, y1(i) + dx*k3, y2(i) + dx*l3,m,T0,n,L,gamma,R,w);
        
            y1(i+1) = y1(i) + (dx/6) * (k1 + 2*k2 + 2*k3 + k4);
            y2(i+1) = y2(i) + (dx/6) * (l1 + 2*l2 + 2*l3 + l4);
    
            uu(i+1) = -(((R*(T0 + m*i*dx)))/(w*Pbar))*y2(i+1);
        end
        if abs(y1(41))<pend
            Y1 = y1;
            Y2 = y2;
            U = uu;
            pend = abs(Y1(41));
            omega(j) = w;
        end
    end    
    p = [p Y1];
    u = [u U];
end

omega/(2*pi)

% Plotting solutions for system of equations
figure(1);
axis equal
plot(x, Y1, 'b', x, Y2, 'r');
title('Solution of Coupled ODEs');
xlabel('x');
ylabel('y1, y2');
legend('y1', 'y2');
grid on;

%plot of pressure amplitudes for different linear axial temperature profiles
figure2=figure(2);
plot(x, abs(p(1:41)), 'k',x, abs(p(42:82)), 'b',x, abs(p(83:123)), 'r',x, abs(p(124:164)), 'g',x, abs(p(165:205)), 'c','linewidth',2);
title('Pressure Amplitude');
xlabel('x');
ylabel('|p|');
legend('m = 0, T_0 = 300','m = -50, T_0 = 500','m = -100, T_0 = 700','m = -150, T_0 = 900','m = -200, T_0 = 1100');
grid on;
saveas(figure2,'pvsx5.png')

%plot of velocity amplitudes for different linear axial temperature profiles
figure3=figure(3);
plot(x, abs(u(1:41)), 'k',x, abs(u(42:82)), 'b',x, abs(u(83:123)), 'r',x, abs(u(124:164)), 'g',x, abs(u(165:205)), 'c','linewidth',2);
title('Velocity Amplitude');
xlabel('x');
ylabel('|v|');
legend('m = -0, T_0 = 300','m = -50, T_0 = 500','m = -100, T_0 = 700','m = -150, T_0 = 900','m = -200, T_0 = 1100');
grid on;
saveas(figure3,'vvsx5.png')

%plot of comparison of pressure and velocity amplitudes
figure4=figure(4);
yyaxis left
plot(x, abs(p(1:41)),'linewidth',2);
xlabel('Distance (m)');
ylabel('Pressure Amplitude (Pa)');
yyaxis right
plot(x, abs(u(1:41)), 'linewidth',2);
ylabel('Velocity Amplitude (m/s)');
legend('Pressure', 'Velocity');
saveas(figure4,'pv5.png') 

% Define the functions f1(x, y1, y2) and f2(x, y1, y2)
function [dy1_dx, dy2_dx] = equations(x, y1, y2,m,T0,n,L,gamma,R,w)
    %define temperature profile
    Tbar = T0 + m*x;
    
    % Define the system of equations
    dy1_dx = y2;
    dy2_dx = -((w^2)/(gamma*R*Tbar)) * y1 - (m/Tbar) * y2;
end