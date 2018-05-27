function [F1, F2, F3, G] = nlc_e(param)

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
    
% W11=param(7), W12=param(8)
f = -x_1./1.0+param(7).*s1+param(8).*s2+xI;
    
% W21=param(9), W22=param(10)    
g = -x_2./1.0+param(9).*s1+param(10).*s2+xI;

% Jacobian of the function wrt x_1 and x_2
J = jacobian([f;g],[x_1,x_2]);

assignin('base','J',J) % Save the jacobian in the workspace

% Pre-location
ev = zeros(2,2);
D = zeros(2,1);
T = zeros(2,1);
SRT = zeros(2,1);
F1 = zeros(2,1);
F2 = zeros(2,1);
F3 = zeros(2,1);
% G = zeros(2,2);
% Xe = zeros(2,2);

s = [0 1; 1 0]; % Desired states (each row is a state)
assignin('base','s',s) % Save in the workspace

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
    
end

    % Nullclines
 
    eq1 = f == 0; 
    eq2 = g == 0;
    
    % Intersections of the nullclines
    [x1e, x2e] = solve(eq1,eq2,x_1,x_2);
    
    Xe(:,1) = double(x1e);
    Xe(:,2) = double(x2e);
    
    G = abs(Xe - s);
        

assignin('base','ev',ev) % Save the eigenvalues in the workspace
assignin('base','Xe',Xe) % Save in the workspace

end