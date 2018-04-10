%% Model for decoding data 

% Model of solutions: Sigmoidal 
% Following: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133

function cf=solutions(param)

global xobs
global t0 
global xpre
global x0 

[t, x]=ode45(@sigmoidal,t0,x0); % The function change according to the model that we want to test

xpre=x;

% residuals ||xobs-xpre|| (20% of importance)
r=sqrt(sum((xpre(:,1)-xobs(:,1)).^2))+sqrt(sum((xpre(:,2)-xobs(:,2)).^2));

% Non-linear constraint related with the stability of the fixed points (60%
% of importance)
% [F1, F2, F3] = nlc(param);

% Constraint related with the feasibility of the problem (20% of importance)
if param(7) || param(8) || param(3) || param(6) == 0
    F4=1;
else
    F4=0;
end

if param(1) || param(4) == 1
    F5=1;
else
    F5=0;
end

% Cost function
% cf=0.2*r+0.2*F1+0.2*F2+0.2*F3+0.1*F4+0.1*F5;
cf=1/3*r+1/3*F4+1/3*F5;

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

end
