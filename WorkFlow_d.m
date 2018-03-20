% VSRP projec: Learning Generative Causal Models from Sparse Temporal Observations during Cellular Reprogramming
% King Abdullah University of Science and Technology (KAUST)

% Developer: Sara Jimenez Correa (Living Systems Laboratory)
% Advisor: Jesper Tegner

clc; clear; close all; 

tic
%% Main

% The following variables are used in the main script and in the "model
% data base" (solutions)

global xobs
global t0 
global xpre
global x0 

%% Fixed points of the "true system" --> xa'=0 and xb'=0 simultaneously

syms xa xb

Xss = myfun(0,[xa;xb]); % Steady state
eq_a = Xss(1) == 0;
eq_b = Xss(2) == 0;    

[xae,xbe]=solve(eq_a,eq_b,xa,xb);

xae=double(xae);
xbe=double(xbe);

%% Computing the vector field for the "true system"

% Initial conditions equally distributed 
xa=linspace(0,3,20);
xb=linspace(0,3,20);

[xa1,xb1]=meshgrid(xa,xb);

% Pre-location
u=zeros(size(xa1));
v=zeros(size(xb1));

for i = 1:numel(xa1)
    Xprime = myfun(0,[xa1(i);xb1(i)]); % Solve at time zero 
    u(i) = Xprime(1);
    v(i) = Xprime(2);    
end

figure(1)
quiver(xb1,xa1,v,u,'r'),xlabel('xb'),ylabel('xa'),title('Vector field of the true system'),axis tight equal;

%% Plotting solutions on the vector field of the "true system"

hold on
for xa0 = [0.5 0.5 0 3.0]
    for xb0 = [3.0 0 0.5 0]
        [t, S] = ode45(@myfun,[0,50],[xa0,xb0]); 
        plot(S(:,2),S(:,1),'b')
    end
end
    
%% Generate data over time

% Initial conditions

t0=[0:0.1:10];
xa0=2.5;
xb0=0.5;
xi=[xa0 xb0];

[t x]=ode45(@myfun,t0,xi);

% Normalization --> values between 0 and 1

xa=x(:,1);
xa=x(:,1)./max(xa);

xb=x(:,2);
xb=x(:,2)./max(xb);

xobs=[xa xb];

% Initial conditions for the fitting --> normalized values
x0=[xobs(1,1) xobs(1,2)];

%% Optimization set-up particle swarm

fun=@solutions; % "model data base" 

% Parameter search space
ub=[1,1,1,1,1,1,1,1,1,1,1,1];
lb=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];

options=optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');

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

[x1e,x2e]=solve(eq_1,eq_2,x1,x2);

x1e=double(x1e);
x2e=double(x2e);

%% Computing the vector field for the "decoded system"

% Initial conditions equally distributed 
x1=linspace(0,3,20);
x2=linspace(0,3,20);

[x1_s,x2_s]=meshgrid(x1,x2);

% Pre-location
u_s=zeros(size(x1_s));
v_s=zeros(size(x2_s));

t=0;
for i = 1:numel(x1_s)
    Xprime_s = sigmoidal_s(t,[x1_s(i);x2_s(i)],param);
    u_s(i) = Xprime_s(1);
    v_s(i) = Xprime_s(2);    
end

figure(3)
quiver(x2_s,x1_s,v_s,u_s,'r'),xlabel('x2'),ylabel('x1'),title('Vector field of the decoded system'),axis tight equal;

hold on
plot(x2e,x1e,'.k','MarkerSize',20) % Fixed point 

%% Plotting solutions on the vector field of the "decoded system"

hold on
for x10 = [0.5 0.5 0 3.0]
    for x20 = [3.0 0 0.5 0]
        [t, S] = ode45(@sigmoidal_s,[0,20],[x10,x20],[],param); 
        plot(S(:,2),S(:,1),'b')
    end
end

toc
%% ODEs for generating data

% Model for generating data
% Following: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29?39. https://doi.org/10.1016/j.bpj.2010.03.058

function dx = myfun(t,x)

% Parameters of the model for generating data

qa=1;
qb=1;
ka=1;
kb=1;
S=0.5;
n=4;

pa=1.2; % High autoactivation of xa
pb=1.2; % High expression of xb

% pa=0.6; % Low autoactivation of xa
% pb=0.6; % Low autoactivation of xb

xa=x(1);
xb=x(2);

dxa=(pa.*xa.^n)./(S^n+xa.^n)+(qa*S^n)./(S^n+xb.^n)-ka.*xa;
dxb=(pb.*xb.^n)./(S^n+xb.^n)+(qb*S^n)./(S^n+xa.^n)-kb.*xb;

dx=[dxa dxb]';

end

%% Model for decoding data 

% Model of solutions: S-Systems
% Following: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133

function r=solutions(param)

global xobs
global t0 
global xpre
global x0 

[t, x]=ode45(@sigmoidal,t0,x0); % The function change according to the model that we want to test

xpre=x;

function dx = sigmoidal(t,x) % 12 parameters 

    x1=x(1);
    x2=x(2);
    
    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1=param(1)+(1-param(1))./(1+exp((-4*param(2)/param(3))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2=param(4)+(1-param(4))./(1+exp((-4*param(5)/param(6))*(x2-param(6))/(1-param(4))));
    % Input
    xI=0;
    
    % W11=param(9), W12=param(10)
    dx1=-x1./param(7)+param(9).*s1+param(10).*s2+xI;
    
    % W21=param(11), W22=param(12)    
    dx2=-x2./param(8)+param(11).*s1+param(12).*s2+xI;
    
    dx=[dx1 dx2]';
    
end

% Cost function ||xobs-xpre||
r=sqrt(sum((xpre(:,1)-xobs(:,1)).^2))+sqrt(sum((xpre(:,2)-xobs(:,2)).^2));

end

function dx = sigmoidal_s(t,x,param) % 12 parameters 

    x1=x(1);
    x2=x(2);
    
    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1=param(1)+(1-param(1))./(1+exp((-4*param(2)/param(3))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2=param(4)+(1-param(4))./(1+exp((-4*param(5)/param(6))*(x2-param(6))/(1-param(4))));
    % Input
    xI=0;
    
    % W11=param(9), W12=param(10)
    dx1=-x1./param(7)+param(9).*s1+param(10).*s2+xI;
    
    % W21=param(11), W22=param(12)    
    dx2=-x2./param(8)+param(11).*s1+param(12).*s2+xI;
    
    dx=[dx1 dx2]';
    
end