%{ 
Model for decoding data 
Model of solutions: Sigmoidal 
Following: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29?39. https://doi.org/10.1016/j.bpj.2010.03.058
%}

function cf = solutions_SH(param)

qa = param(1);
qb = param(2);
ka = param(3);
kb = param(4);
S = param(5);
n = param(6);
pa = param(7);
pb = param(8);

% Non-linear constraint related with the stability of the fixed points 
[F1, F2, F3, G] = nlc_SH(param);
    
% Cost function
cf = 1*F1(1,1) + 1*F1(2,1) + 1*F2(1,1) + 1*F2(2,1) + 1*F3(1,1) + 1*F3(2,1) + 1*G(1,1) + 1*G(2,1);

function [F1, F2, F3, G] = nlc_SH(param)

qa = param(1);
qb = param(2);
ka = param(3);
kb = param(4);
S = param(5);
n = param(6);
pa = param(7);
pb = param(8);
    
x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions

f = (pa.*x_1.^n)./(S^n+x_1.^n)+(qa*S^n)./(S^n+x_2.^n)-ka.*x_1;
        
g = (pb.*x_2.^n)./(S^n+x_2.^n)+(qb*S^n)./(S^n+x_1.^n)-kb.*x_1;

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
 
eq1 = ((pa.*x_1.^n)./(S^n+x_1.^n)+(qa*S^n)./(S^n+x_2.^n))./ka - x_1;
  
eq2 = ((pb.*x_2.^n)./(S^n+x_2.^n)+(qb*S^n)./(S^n+x_1.^n))./kb - x_2;

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
