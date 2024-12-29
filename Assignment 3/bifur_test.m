N = 1000; %number of timesteps
T = 40; %total time of experiment
dt = (T/N);
xf = 0.3; %m

%define system parameters
g = 1.4;
c0 = 399.6; %speed of sound
u0 = 0.5; %mean velocity of air
M = u0/c0;

J = 3; %number of Galerkin modes

% Set initial values
ks=0:0.05:2;%non-dimensional heater power
n1=length(ks)

%tau=0.2:0.1:0.8;%time-lag
%n3=length(tau);
tau = 0.2;

% Preallocate arrays
u_rms1 = zeros(1,n1);

u1 = zeros(n1, N);
p1 = zeros(n1, N);

y1 = zeros(J, N);
y2 = zeros(J, N);

y1(1,1)=0.18;

%forward
for k=1:n1
    % Set initial values
    for j = 1:J  
        u1(k,1) = u1(k,1) + y1(j,1)*cos(j*pi*xf);
    end

    n2 = round(tau/dt);
    % Fourth-order Runge-Kutta method
    for n = 1:n2
      for j = 1:J  
        [k1, l1] = equations1(y1(j,n), y2(j,n),j);
        [k2, l2] = equations1(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j);
        [k3, l3] = equations1(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j);
        [k4, l4] = equations1(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u1(k,n+1) = u1(k,n+1) + y1(j,n+1)*cos(j*pi*xf);
      end
    end

    for n = n2+1:N-1
      for j=1:J
        [k1, l1] = equations2(y1(j,n), y2(j,n),j,ks(k),u1(k,n-n2),xf);
        [k2, l2] = equations2(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j,ks(k),u1(k,n-n2),xf);
        [k3, l3] = equations2(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j,ks(k),u1(k,n-n2),xf);
        [k4, l4] = equations2(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j,ks(k),u1(k,n-n2),xf);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u1(k,n+1) = u1(k,n+1) + y1(j,n+1)*cos(j*pi*xf);
      end
    end

    %figure(k);
    %plot(0:N-1, u1(k,:), 'k','linewidth',2);
    %xlabel('t');
    %ylabel('|u''|');
    %grid on;

    %getting initial values for next
    y1(:,1) = y1(:,N);
    y2(:,1) = y2(:,N);
 
    u_rms1(k) = rms(u1(k,500:N));
end

% Preallocate arrays
u_rms2 = zeros(1,n1);

u2 = zeros(n1, N);
p2 = zeros(n1, N);

%backward 
for k=1:n1
    % Set initial values
    %y1(:,1) = y1(:,N);
    %y2(:,1) = y2(:,N);
     for j = 1:J  
          u2(n1-k+1,1) = u2(n1-k+1,1) + y1(j,1)*cos(j*pi*xf);
    end
    
    n2 = round(tau/dt);
    % Fourth-order Runge-Kutta method
    for n = 1:n2     
      for j=1:J
        [k1, l1] = equations1(y1(j,n), y2(j,n),j);
        [k2, l2] = equations1(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j);
        [k3, l3] = equations1(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j);
        [k4, l4] = equations1(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u2(n1-k+1,n+1) = u2(n1-k+1,n+1) + y1(j,n+1)*cos(j*pi*xf);
      end
    end

    for n = n2+1:N-1
      for j=1:J
        [k1, l1] = equations2(y1(j,n), y2(j,n),j,ks(n1-k+1),u2(n1-k+1,n-n2),xf);
        [k2, l2] = equations2(y1(j,n) + 0.5*dt*k1, y2(j,n) + 0.5*dt*l1,j,ks(n1-k+1),u2(n1-k+1,n-n2),xf);
        [k3, l3] = equations2(y1(j,n) + 0.5*dt*k2, y2(j,n) + 0.5*dt*l2,j,ks(n1-k+1),u2(n1-k+1,n-n2),xf);
        [k4, l4] = equations2(y1(j,n) + dt*k3, y2(j,n) + dt*l3,j,ks(n1-k+1),u2(n1-k+1,n-n2),xf);
    
        y1(j,n+1) = y1(j,n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        y2(j,n+1) = y2(j,n) + (dt/6) * (l1 + 2*l2 + 2*l3 + l4);

        u2(n1-k+1,n+1) = u2(n1-k+1,n+1) + y1(j,n+1)*cos(j*pi*xf);
      end
    end

    %figure(k);
    %plot(0:N-1, u2(k,:), 'k','linewidth',2);
    %xlabel('t');
    %ylabel('|u''|');
    %grid on;
 
    u_rms2(n1-k+1) = rms(u2(n1-k+1,500:N));
end



figure(1);
plot(ks,u_rms1,'b*','linewidth',2);
hold on;
plot(ks,u_rms2,'rx','linewidth',2);
xlabel('K');
ylabel('|u_1|');
legend('Forward','Backward')
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
    f2 = -2*zetaj*wj*y2 - (kj^2) * y1;
end

% for t > tau
function [f1,f2] = equations2(y1, y2, j, K, u,xf)
    kj = j*pi;
    wj = kj;
    w1 = pi;
    c1 = 0.1; c2 = 0.06;
    zetaj = (1/(2*pi))*(c1*(wj/w1) + c2*sqrt(w1/wj));

    f1 =  y2;
    f2 = -2*zetaj*wj*y2 - (kj^2) * y1 - j*pi*K*(abs(sqrt((1/3) + u)) - sqrt(1/3))*sin(j*pi*xf) ;
end

