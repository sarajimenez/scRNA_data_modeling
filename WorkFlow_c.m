% VSRP projec: Learning Generative Causal Models from Sparse Temporal Observations during Cellular Reprogramming
% King Abdullah University of Science and Technology (KAUST)

% Developer: Sara Jimenez Correa (Living Systems Laboratory)
% Advisor: Jesper Tegner

clc; clear; close all; 

tic
%% Main

% The following variables are used in the main script and in the "model
% data base" script (solutions)

global xobs
global t0 
global xpre
global x0 

%% Fixed points --> xa'=0 and xb'=0 simultaneously (these are equilibrium points)

syms xa xb

% Parameters of the model for generating data "True system"

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

eq_a=(pa.*xa.^n)./(S^n+xa.^n)+(qa*S^n)./(S^n+xb.^n)-ka.*xa==0;
eq_b=(pb.*xb.^n)./(S^n+xb.^n)+(qb*S^n)./(S^n+xa.^n)-kb.*xb==0;

[xae,xbe]=solve(eq_a,eq_b,xa,xb);

xae=double(xae);
xbe=double(xbe);


%% Computing the vector field for the "True system"

xa=linspace(0,3,20);
xb=linspace(0,3,20);

[xa1,xb1]=meshgrid(xa,xb);

u=zeros(size(xa1));
v=zeros(size(xb1));

t=0;
for i = 1:numel(xa1)
    Xprime = myfun(t,[xa1(i);xb1(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);    
end

figure(1)
quiver(xb1,xa1,v,u,'r'); 
xlabel('xb')
ylabel('xa')
axis tight equal;

%% Plotting solutions on the vector field of the "true system"

hold on
for xa0 = [0.5 0.5 0 2.5]
    for xb0 = [2.5 0 0.5 0]
        [t, S] = ode45(@myfun,[0,25],[xa0,xb0]); 
        plot(S(:,2),S(:,1),'b')
    end
end
    

%% Generate data over time

% Initial conditions

t0=[0:0.1:20];
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

x0=[xobs(1,1) xobs(1,2)];


%% Optimization set-up particle swarm

fun=@solutions; % "model data base" script

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

toc
%% ODEs for generating data

% Model for generating data
% Following: Stefan Semrau and Alexander van Oudenaarden, Annu. Rev. Cell Dev. Biol. 2015. 31:317-345

function dx = myfun(t,x)

% Parameters of the model for generating data "True system"

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
% Following: Daniels, B. C., & Nemenman, I. (2015). Automated adaptive inference of phenomenological dynamical models. Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133

function r=solutions(param)

global xobs
global t0 
global xpre
global x0 

[t, x]=ode45(@sigmoidal,t0,x0); % The function change according to the model that we want to test

xpre=x;

function dx = SSystems(t,x) % 10 parameters

    x1=x(1);
    x2=x(2);
    
    xI=0;
    
    dx1=param(1).*xI.^(param(2)).*x1.^(param(3)).*x2.^(param(4))-param(5).*xI.^(param(6)).*x1.^(param(7)).*x2.^(param(8));
    dx2=x1.^(param(9)).*x2.^(param(10))-1;
    
    dx=[dx1 dx2]';
    
end

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