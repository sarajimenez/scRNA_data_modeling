%{
VSRP projec: Learning Generative Causal Models from Sparse Temporal Observations during Cellular Reprogramming
King Abdullah University of Science and Technology (KAUST)

Developer: Sara Jimenez Correa (Living Systems Laboratory)
Advisor: Jesper Tegner 
%}

clc; clear; close all; 

tic
%% Main

% The following variables are used in the main script and target function (solutions_b.m)

global t0 
global x0 
global xpre

% Initial conditions of the ODE's
t0 = [0:0.1:5];
x0 = [0  1];

%% Optimization set-up particle swarm

fun = @solutions_b; % "model data base" 

% Parameter search space
ub = [0.9,2,2,0.9,2,2,1,1,1,1,1,1];
lb = [0,0,0.1,0,0,0.1,0.1,0.1,-1,-1,-1,-1];

options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');

rng default  % For reproducibility
nvars = 12; % Number of parameters to estimate 
[param] = particleswarm(fun,nvars,lb,ub,options);
        
figure(2)
subplot(1,2,1), plot(t0,xpre(:,1),'k'),title('x_A'),legend('Predicted'),xlabel('time'),ylabel('Expression'); hold on;
subplot(1,2,2), plot(t0,xpre(:,2),'k'),title('x_B');legend('Predicted'),xlabel('time'),ylabel('Expression'); hold on;


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
x1 = linspace(0,2,20);
x2 = linspace(0,2,20);

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
for x10 = [0 2.0]
    for x20 = [0 0.3 0.5 1.0 1.5 2.0]
        [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

hold on
for x20 = [0 2.0]
    for x10 = [0 0.3 0.5 1.0 1.5 2.0]
        [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

toc
