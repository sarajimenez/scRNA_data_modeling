%{ 
Model for decoding data 
Model of solutions: Sigmoidal 
Following: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133 
%}

function cf = solutions_d(param)

% global t0 
% global x0 
% global xpre
% 
% The initial conditions can be changed 
% [t, x]=ode45(@sigmoidal_s,t0,x0,param); % The function change according to the model that we want to test
% 
% xpre = x;

% Non-linear constraint related with the stability of the fixed points 
[F1, F2, F3, G] = nlc_d(param);

% Constraint related with the relation of the parameters

H1 = exp(param(9)*param(10));
H2 = exp(param(11)*param(12));
H4 = exp(-param(10)*param(11));
H5 = exp(-param(9)*param(12));
H6 = exp(-(param(2)/(2*param(3))-1));
H7 = exp(-(param(5)/(2*param(6))-1));
    
% Cost function
cf = 1*F1(1,1) + 1*F1(2,1) + 1*F2(1,1) + 1*F2(2,1) + 1*F3(1,1) + 1*F3(2,1) + 1*G(1,1) + 1*G(2,1) + 1*H1 + 1*H2 + 1*H4 + 1*H5 + 1*H6 + 1*H7;

function [F1, F2, F3, G] = nlc_d(param)

x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions

% sigmoids 

% Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
s1 = param(1) + 1/(1+exp(-param(3)*(x_1-param(2)/param(3))));
% Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
s2 = param(4) + 1/(1+exp(-param(6)*(x_2-param(5)/param(6))));

% Input
xI = 0;
    
% W11=param(9), W12=param(10)
f = -x_1./param(7)+param(9).*s1+param(10).*s2+xI;
    
% W21=param(11), W22=param(12)    
g = -x_2./param(8)+param(11).*s1+param(12).*s2+xI;

Xj = [f; g];

assignin('base','Xj',Xj) % Save the jacobian in the workspace

% Jacobian of the function wrt x_1 and x_2
J = jacobian([Xj(1);Xj(2)],[x_1,x_2]);

assignin('base','J',J) % Save the jacobian in the workspace

% Pre-location
ev = zeros(2,2);
D = zeros(2,1);
T = zeros(2,1);
SRT = zeros(2,1);
F1 = zeros(2,1);
F2 = zeros(2,1);
F3 = zeros(2,1);
G = zeros(2,1);
eq1_e = zeros(2,1);
eq2_e = zeros(2,1);

% Nullclines
 
% W11=param(9), W12=param(10)
eq1 = param(7).*(param(9).*s1+param(10).*s2) - x_1;
% W21=param(11), W22=param(12)    
eq2 = param(8).*(param(11).*s1+param(12).*s2) - x_2;

s = [0.0 1.0; 1.0 0.0]; % Desired states (each row is a state)

for j = 1:2
       
    % Evaluate the jacobian in the "desired states"
    J_e = subs(J,{x_1,x_2},{s(j,1),s(j,2)}); % Here we modify the states  
    J_e = double(J_e);
    
    % Eigen values 
    ev(:,j) = eig(J_e);
    
    % The determinant of the jacobian evaluated at the desired state must be positive 
    D(j) = ev(1,j).*ev(2,j);
    
    F1(j) = exp(-D(j));
    
    % The trace of the jacobian evaluated at the desired state must be negative
    T(j) = ev(1,j)+ev(2,j);
    
    F2(j) = exp(T(j));
    
    % The eigenvalues must be real then the term inside the square root must be positive
    SRT(j) = T(j)^2-4.*D(j);
    
    F3(j) = exp(-SRT(j));
    
    % Intersections of the nullclines
    eq1_e(j) = subs(eq1,{x_1,x_2},{s(j,1),s(j,2)}); eq1_e(j) = double(eq1_e(j));
    
    eq2_e(j) = subs(eq2,{x_1,x_2},{s(j,1),s(j,2)}); eq2_e(j) = double(eq2_e(j));
    
    G(j) = abs(eq1_e(j) - eq2_e(j));
        
end

assignin('base','ev',ev) % Save the eigenvalues in the workspace
assignin('base','eq1_e',eq1_e) % Save in the workspace
assignin('base','eq2_e',eq2_e) % Save in the workspace


end

end
