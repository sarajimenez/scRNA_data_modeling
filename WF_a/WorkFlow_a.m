%{
VSRP projec: Learning Generative Causal Models from Sparse Temporal Observations during Cellular Reprogramming
King Abdullah University of Science and Technology (KAUST)

Developer: Sara Jimenez Correa (Living Systems Laboratory)
Advisor: Jesper Tegner 
%}

clc; clear; close all; 

tic
%% Main

% The following variables are used in the main script and target function (solutions_a.m)

global xobs
global t0 
global xpre
global x0 

%% Fixed points of the "true system" (generated data) --> xa'=0 and xb'=0 simultaneously

syms xa xb

Xss = myfun(0,[xa;xb]); % Steady state
eq_a = Xss(1) == 0;
eq_b = Xss(2) == 0;    

[xae,xbe] = solve(eq_a,eq_b,xa,xb);

xae = double(xae);
xbe = double(xbe);

%% Computing the vector field for the "true system" (generated data)

% Initial conditions equally distributed 
xa = linspace(0,3,20);
xb = linspace(0,3,20);

[xa1,xb1] = meshgrid(xa,xb);

% Pre-location
u = zeros(size(xa1));
v = zeros(size(xb1));

for i = 1:numel(xa1)
    Xprime = myfun(0,[xa1(i);xb1(i)]); % Solve at time zero 
    u(i) = Xprime(1);
    v(i) = Xprime(2);    
end

figure(1)
quiver(xb1,xa1,v,u,'r'),xlabel('xb'),ylabel('xa'),title('Vector field of the true system'),axis tight equal;

%% Plotting solutions on the vector field of the "true system" (generated data)

hold on
for xa0 = [0 3.0]
    for xb0 = [0 0.05 0.075 0.1 0.3 0.5 0.75 1.0 1.5 2.0 3.0]
        [t, S] = ode45(@myfun,[0,10],[xa0,xb0]); 
        plot(S(:,2),S(:,1),'b')
    end
end

hold on
for xb0 = [0 3.0]
    for xa0 = [0 0.05 0.075 0.1 0.3 0.5 0.75 1.0 1.5 2.0 3.0]
        [t, S] = ode45(@myfun,[0,10],[xa0,xb0]); 
        plot(S(:,2),S(:,1),'b')
    end
end
    
%% Generate data over time

% Initial conditions

t0 = [0:0.1:10];
xa0 = 0; 
xb0 = 0;

[t x] = ode45(@myfun,t0,[xa0 xb0]);

% Normalization --> values between 0 and 1

m = max(max(x(:,1)),max(x(:,2))); % Maximum value 

xa = x(:,1)./m;
xb = x(:,2)./m;

xobs = [xa xb];

% Initial conditions for the fitting --> normalized values
x0 = [xobs(1,1) xobs(1,2)];

%% Optimization set-up particle swarm

fun = @solutions_a; % "model data base" 

% Parameter search space
ub = [0.9,1,1,0.9,1,1,1,1,1,1,1,1];
lb = [0,0,0.1,0,0,0.1,0.1,0.1,-1,-1,-1,-1];

options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');

rng default  % For reproducibility
nvars = 12; % Number of parameters to estimate 
[param] = particleswarm(fun,nvars,lb,ub,options);
        
figure(2)
subplot(1,2,1), plot(t,xobs(:,1),'r',t,xpre(:,1),'k'),title('x_A'),legend('Observed','Predicted'),xlabel('time'),ylabel('Expression'); hold on;
subplot(1,2,2), plot(t,xobs(:,2),'r',t,xpre(:,2),'k'),title('x_B');legend('Observed','Predicted'),xlabel('time'),ylabel('Expression'); hold on;


%% Fixed points of the solution system --> x1'=0 and x2'=0 simultaneously 

syms x1 x2

Xss_s = sigmoidal_s(0,[x1;x2],param); % Steady state of the solution system 
eq_1 = Xss_s(1) == 0;
eq_2 = Xss_s(2) == 0;    

[x1e,x2e] = solve(eq_1,eq_2,x1,x2);

x1e = double(x1e);
x2e = double(x2e);

%% Computing the vector field for the "decoded system"

% Initial conditions equally distributed 
x1 = linspace(0,3,20);
x2 = linspace(0,3,20);

[x1_s,x2_s] = meshgrid(x1,x2);

% Pre-location
u_s = zeros(size(x1_s));
v_s = zeros(size(x2_s));

for i = 1:numel(x1_s)
    Xprime_s = sigmoidal_s(0,[x1_s(i);x2_s(i)],param); % Solve at time zero 
    u_s(i) = Xprime_s(1);
    v_s(i) = Xprime_s(2);    
end

figure(4)
quiver(x2_s,x1_s,v_s,u_s,'r'),xlabel('x2'),ylabel('x1'),title('Vector field of the decoded system'),axis tight equal;

hold on
plot(x2e,x1e,'.k','MarkerSize',20) % Fixed point 

%% Plotting solutions on the vector field of the "decoded system"

hold on
for x10 = [0 3.0]
    for x20 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,10],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

hold on
for x20 = [0 3.0]
    for x10 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,10],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

toc
