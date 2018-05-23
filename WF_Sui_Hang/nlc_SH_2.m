function [F1, F2, F3, G] = nlc_SH_2(param)

qa = 1;
qb = 1;
ka = 1;
kb = 1;
S = 0.5;
n = 4;
pa = param(1);
pb = param(2);
    
x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions

f = (pa.*x_1.^n)./(S^n+x_1.^n)+(qa*S^n)./(S^n+x_2.^n)-ka.*x_1;
        
g = (pb.*x_2.^n)./(S^n+x_2.^n)+(qb*S^n)./(S^n+x_1.^n)-kb.*x_1;

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
%G = zeros(2,1);

s = [0.01 1.6; 1.6 0.01]; % Desired states (each row is a state)

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
    
    % Nullclines
    
    eq1 = f == 0; 
    eq2 = g == 0;
    
    % Intersections of the nullclines
    [x1e, x2e] = solve(eq1,eq2,x_1,x_2);
    
    x1e=double(x1e);
    x2e=double(x2e);
    
    Xe(j,1) = double(x1e);
    Xe(j,2) = double(x2e);
    
    G(j,:) = abs(Xe(j,:) - s(j,:));
        
end

assignin('base','ev',ev) % Save the eigenvalues in the workspace
assignin('base','Xe',Xe) % Save in the workspace


end