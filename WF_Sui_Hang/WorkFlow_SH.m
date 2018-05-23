%{
VSRP projec: Learning Generative Causal Models from Sparse Temporal Observations during Cellular Reprogramming
King Abdullah University of Science and Technology (KAUST)

Developer: Sara Jimenez Correa (Living Systems Laboratory)
Advisor: Jesper Tegner 
%}

clc; clear; close all;

tic
%% Main

% Optimization set-up particle swarm

fun = @solutions_SH; % "model data base" 

% Parameter search space
ub = [2,2,2,2,1,6];
lb = [0,0,0,0,0,2];

options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');
% options = optimoptions('particleswarm','SwarmSize',100,'Display','iter');

rng default  % For reproducibility
nvars = 6; % Number of parameters to estimate 
[param, exitflag, fval, output] = particleswarm(fun,nvars,lb,ub,options);

  
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

figure(1)
quiver(x2_s,x1_s,v_s,u_s,'r'),xlabel('x2'),ylabel('x1'),title('Vector field of the decoded system'),axis tight equal;

hold on
plot(x2e,x1e,'.k','MarkerSize',20) % Fixed point 

%% Plotting solutions on the vector field of the "decoded system"

hold on
for x10 = [0 3.0]
    for x20 = [0 0.05 0.075 0.1 0.3 0.5 0.75 1.0 1.5 2.0 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,15],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

hold on
for x20 = [0 3.0]
    for x10 = [0 0.05 0.075 0.1 0.3 0.5 0.75 1.0 1.5 2.0 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,15],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

toc
